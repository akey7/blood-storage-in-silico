module UfbaSampler

using Distributed
using COBREXA, HiGHS, JuMP, MathOptInterface
using Base.Iterators
import ConstraintTrees as C
import SBMLFBCModels as S
import AbstractFBCModels as A
import AbstractFBCModels: stoichiometry
import AbstractFBCModels.CanonicalModel: Model, Reaction, Metabolite, Gene, Coupling
using CSV
using DataFrames
using DataFramesMeta
using ThreadsX
using OrderedCollections

include("FbaModelBuilder.jl")
using .FbaModelBuilder

export sample_fluxes,
    ufba_all_additives_all_times,
    load_metabolite_bounds,
    query_metabolite_bounds,
    histograms_for_reaction_in_additive,
    plot_all_histograms,
    fba,
    add_sinks_for_unmatched_metabolites!,
    find_metabolite_matches,
    is_metabolite_in_exchange,
    case_3_constraint_tree!,
    list_objectives_in_model,
    optimize_case_3,
    display_jump_results,
    analyze_case_3,
    make_ufba_models_for_additives_and_times,
    execute_all_ufba_jobs,
    map_reaction_ids_to_reaction_strings,
    extract_case3_sinks,
    init_workers!,
    execute_ufba_job,
    count_n_all_zero_fluxes,
    does_manual_prune_list_match_sink_name,
    load_flux_bounds_overrides,
    sbml_add_constant_to_selfclosing_parameters!,
    extract_added_case3_sink_ids,
    find_metabolites_with_exchanges


"""
    init_workers!(; project=Base.active_project())

Activate `project` and load required packages on all current workers. This is necessary
for multiprocessing in uFBA sampling.

Call this after `addprocs(...)` (or any time you add more workers).
"""
function init_workers!(; project::AbstractString = Base.active_project())
    for p in workers()
        Distributed.remotecall_eval(
            Main,
            p,
            quote
                import Pkg
                Pkg.activate($project)
                using COBREXA, HiGHS, JuMP, MathOptInterface
            end,
        )
    end
    return nothing
end

"""
    map_reaction_ids_to_reaction_strings(model::A.AbstractFBCModel)

Maps reaction_ids in the given model to human-readable reaction strings specifying reactants and products with an arrow pointing in the direction specified by the bounds of the reaction.

# Arguments
1. `model::A.AbstractFBCModel`: The model to create the reaction strings from.

# Returns
`Dict{String,String}`

Returns a dicitonary mapping reaction ids in the model to a human-readable reaction string.
"""
function map_reaction_ids_to_reaction_strings(model::A.AbstractFBCModel)
    result = OrderedDict()
    for rxn_id in sort(string.(keys(model.reactions)))
        stoi = model.reactions[rxn_id].stoichiometry
        rxn = model.reactions[rxn_id]
        lhs = replace(
            join(
                [!isapprox(v, -1.0) ? "$(abs(v)) $k" : k for (k, v) in stoi if v < 0],
                " + ",
            ),
            "M_" => "",
        )
        rhs = replace(
            join(
                [!isapprox(v, 1.0) ? "$(abs(v)) $k" : k for (k, v) in stoi if v > 0],
                " + ",
            ),
            "M_" => "",
        )
        if rxn.lower_bound < 0.0 && isapprox(rxn.upper_bound, 0.0)
            result[rxn_id] = "$lhs <-- $rhs ($(rxn.lower_bound), $(rxn.upper_bound))"
        elseif isapprox(rxn.lower_bound, 0.0) && rxn.upper_bound > 0.0
            result[rxn_id] = "$lhs --> $rhs ($(rxn.lower_bound), $(rxn.upper_bound))"
        else
            result[rxn_id] = "$lhs <-> $rhs ($(rxn.lower_bound), $(rxn.upper_bound))"
        end
    end
    return result
end

"""
    load_metabolite_bounds()

Loads the rates of metabolite oncentration change from the `output/concentration_rates.csv` file. This file is produced by the `AbsoluteQuant` module from absolute (or approximately absolute) metabolomics quantifcation data over time.

Downstream handling of this DataFrame expects to find the following columns in the csv: additive, metabolite, final_time, intercept, rate, lb, ub.

# Returns
`DataFrame`

Returns the loaded DataFrame.
"""
function load_metabolite_bounds()
    metabolite_bounds_filename = joinpath("output", "concentration_rates.csv")
    metabolite_bounds_df = CSV.read(metabolite_bounds_filename, DataFrame)
    return metabolite_bounds_df
end

"""
    load_flux_bounds_overrides()

Loads `input/flux_bounds_overrides.csv`. This file contains rate bounds for fluxes that **override** the specifications in the RBC-GEM. If the `.csv` file does not exist, then `nothing` is returned.

# Returns
`Union{DataFrame,Nothing}`

Returns the flux bounds overrides DataFrame or `nothing` if the source file does not exist.
"""
function load_flux_bounds_overrides()
    flux_bounds_filename = joinpath("input", "flux_bounds_overrides.csv")
    return isfile(flux_bounds_filename) ? CSV.read(flux_bounds_filename, DataFrame) :
           nothing
end

"""
    sbml_add_constant_to_selfclosing_parameters!(infile::AbstractString; outfile::AbstractString = infile, default_constant::AbstractString = "true")

This is a patch because COBREXA is writing corrupt SBML files. This opens the file and fixes the problem.

# Arguments
1. `infile::AbstractString`: Filename to patch.
2. `outfile::AbstractString = infile`: Out file to write
3. `default_constant::AbstractString = "true"`: Constant to patch with.
"""
function sbml_add_constant_to_selfclosing_parameters!(
    infile::AbstractString;
    outfile::AbstractString = infile,
    default_constant::AbstractString = "true",
)
    s = read(infile, String)
    s2 = replace(
        s,
        Regex(raw"<parameter\b(?![^>]*\bconstant=)([^>]*)\s*/>") =>
            SubstitutionString("<parameter\\1 constant=\"$default_constant\"/>"),
    )
    write(outfile, s2)
end

"""
    save_ufba_model_sbml(model::A.AbstractFBCModel, additive::AbstractString, final_time::Int64)

Save the given uFBA model to the filesystem for later retrieval. Models are saved in SBML format in the `output/ufba_models` folder.

# Arguments:
1. `model::A.AbstractFBCModel`: uFBA model to save.
2. `additive::AbstractString`: Additive the uFBA model is in.
3. `final_time::Int64`: Final time of the uFBA model.
"""
function save_ufba_model_sbml(
    model::A.AbstractFBCModel,
    additive::AbstractString,
    final_time::Int64,
)
    filename = joinpath("output", "ufba_models", "uFBA $(additive)_$(final_time).xml")
    sbml_fbc = convert(S.SBMLFBCModel, model)
    save_model(sbml_fbc, filename)
    sbml_add_constant_to_selfclosing_parameters!(filename)
    println("Wrote $filename")
end

"""
    query_metabolite_bounds(metabolite_bounds_df, additive, metabolite, final_time)

Find the rate of concentration chage for the metabolite in the given additive at the given final time. Returns `nothing` if not found.

# Arguments
1. `metabolite_bounds_df`: DataFrame as loaded by [`load_metabolite_bounds`](@ref BloodStorageInSilico.UfbaSampler.load_metabolite_bounds).
2. `additive`: String of the additive as specified in the DataFrame.
3. `metabolite`: Metabolite id.
4. `final_time`: The final time point of the interval.

# Returns
`Tuple{Float64,Float64}`

Using the 95% confidence interval of rate in the original DataFrame, a tuple with the lower and upper bounds of this interval.
"""
function query_metabolite_bounds(metabolite_bounds_df, additive, metabolite, final_time)
    query_df = @rsubset(
        metabolite_bounds_df,
        :additive == additive,
        :metabolite == metabolite,
        :final_time == final_time
    )
    if nrow(query_df) > 0
        lower_bound = query_df[1, :lb]
        upper_bound = query_df[1, :ub]
        return (lower_bound, upper_bound)
    else
        return nothing
    end
end

"""
    fba(model::A.AbstractFBCModel; n_chains::Int64 = 10)

Standard flux balance analysis of the given `model`. Returns samples of fluxes upon success, `nothing` for infeasible solutions.

# Arguments
1. `model::A.AbstractFBCModel`: The model to optimize.

2. `n_chains::Int64 = 10`: Number of chains to sample. Must be a keyword and defaults to 10.

# Returns
`Union{Nothing,DataFrame}`

Returns samples of fluxes upon success, `nothing` for infeasible solutions.
"""
function fba(model::A.AbstractFBCModel; n_chains::Int64 = 10)
    @info "Standard FBA sampling, N chains $n_chains"
    println("> Simple optimization attempt")
    solution = flux_balance_analysis(model; optimizer = HiGHS.Optimizer)
    if isnothing(solution)
        println("Simple optimization failed")
        return nothing, nothing
    else
        println("Simple optimization succeeded!")
        display(solution.fluxes)
        println("> Flux sampling")
        samples_df = sample_fluxes(model; n_chains = n_chains)
        return solution, samples_df
    end
end

"""
    is_metabolite_in_exchange(model::A.AbstractFBCModel, metabolite::AbstractString)

Determines whether a metabolite is in an exchange by detecting a substring in the id of the reaction in which the metabolite is found.

# Arguments
1. `model::A.AbstractFBCModel`: The model with the reactions to check.

2. `metabolite::AbstractString`: Metabolite id to search for.

# Returns
`Bool`

`true` if the metabolite is in an exchange, `false` otherwise.
"""
function is_metabolite_in_exchange(model::A.AbstractFBCModel, metabolite::AbstractString)
    exchange_substring = "EX_$(metabolite[1:end-2])"
    for rxn in keys(model.reactions)
        if contains(rxn, exchange_substring)
            return true
        end
    end
    return false
end

"""
    find_metabolite_matches(model::A.AbstractFBCModel, metabolite_bounds_df::DataFrame, additive::AbstractString, final_time::Int64)

Creates a DataFrame of the metabolites found, not found, or in exchange for each additive and time point in the flux balance constratint tree for the given model. This is useful for determining which metabolites have been measured and are available at each time point for each additive. In other words, this is a data quality check function.

# Arguments
1. `model::A.AbstractFBCModel`: The model to get the `ConstraintTree` from.

2. `metabolite_bounds_df::DataFrame`: The metabolite bounds DataFrame to search.

3. `additive::AbstractString`: Additive being searched.

4. `final_time::Int64`: Final time being searched.

# Returns
`DataFrame`

Returns a `DataFrame` of metabolite measurement availability.
"""
function find_metabolite_matches(
    model::A.AbstractFBCModel,
    metabolite_bounds_df::DataFrame,
    additive::AbstractString,
    final_time::Int64,
)
    @info "Matching metabolites, additive: $additive, final_time: $final_time"
    ct = flux_balance_constraints(model)
    status_rows = []
    found_count = 0
    not_found_count = 0
    in_exchange_count = 0
    for k ∈ keys(ct.flux_stoichiometry)
        short_metabolite_id = string(k)[3:end]
        bounds = query_metabolite_bounds(
            metabolite_bounds_df,
            additive,
            short_metabolite_id,
            final_time,
        )
        if isnothing(bounds)
            status_row = (
                additive = additive,
                metabolite = short_metabolite_id,
                status = "not found",
            )
            push!(status_rows, status_row)
            not_found_count += 1
            # ct.flux_stoichiometry[k].bound = C.Between(-1000.0, 1000.0)
        elseif is_metabolite_in_exchange(model, short_metabolite_id)
            status_row = (
                additive = additive,
                metabolite = short_metabolite_id,
                status = "in exchange",
            )
            in_exchange_count += 1
            lb, ub = bounds
            if isapprox(lb, 0.0) && isapprox(ub, 0.0)
                @warn "$additive $short_metabolite_id is fixed at 0.0"
            end
        else
            status_row =
                (additive = additive, metabolite = short_metabolite_id, status = "found")
            push!(status_rows, status_row)
            found_count += 1
            lb, ub = bounds
            if isapprox(lb, 0.0) && isapprox(ub, 0.0)
                @warn "$additive $short_metabolite_id is fixed at 0.0"
            end
        end
    end
    metabolite_status_df = DataFrame(status_rows)
    println(
        "Found $found_count, in exchange $in_exchange_count, not found $not_found_count",
    )
    return metabolite_status_df
end

"""
    does_manual_prune_list_match_sink_name(sink_name::String, sink_opt_outs::Union{Vector{String},Nothing} = nothing)

Determines if the given sink name contains any of the substrings in the given sink opt-outs list.

# Arguments
1. `sink_name::String`: The name of the sink.
2. `sink_opt_outs::Union{Vector{String},Nothing} = nothing`: If specified, contains a list of substrings that are matched against the given sink name.

# Returns
`Bool`

Returns `true` if one of the provided substrings matches the given sink name. Returns `false` if the substring list is not provided or none of the substrings are found
"""
function does_manual_prune_list_match_sink_name(
    sink_name::String,
    sink_opt_outs::Union{Vector{String},Nothing} = nothing,
)
    if isnothing(sink_opt_outs)
        return false
    else
        for sink_opt_out in sink_opt_outs
            if contains(sink_name, sink_opt_out)
                return true
            end
        end
        return false
    end
end

function find_metabolites_with_exchanges(model::A.AbstractFBCModel)
    exchange_ids = [
        reaction_id for (reaction_id, _) in model.reactions if contains(reaction_id, "R_EX")
    ]
    metabolite_ids = [replace(exchange_id, "R_EX_" => "") for exchange_id in exchange_ids]
    return metabolite_ids
end

"""
    add_sinks_for_unmatched_metabolites!(model::A.AbstractFBCModel, metabolite_status_df::DataFrame, additive::AbstractString, prune_zero_sinks::Union{Vector{String},Nothing}; sink_opt_outs::Union{Vector{String},Nothing} = nothing)

Add sinks for unmeasured (umatched) metabolites in the model. This is part of the uFBA process. This method mutates the given model in place.

# Arguments
1. `model::A.AbstractFBCModel`: Model to add sinks to. **This model is mutated in place.**
2. `metabolite_status_df::DataFrame`: Metabolite measurement availability DataFrame.
3. `additive::AbstractString`: Additive for measurement search.
4. `prune_zero_sinks::Union{Vector{String},Nothing}`: If specified, the provided list of zero flux sinks are not added (pruned) to the model. If `nothing`, no sinks are pruned.
5. `sink_opt_outs::Union{Vector{String},Nothing} = nothing`: If a `Vector{String}`, sinks with specified substrings are ensured to not be added to the model. For example, placing `2pg_c` in this list will ensure that NO sink for `2pg_c` will be added. This parameter provides another way to manually opt-out of sinks, rather than simply relying on the automated zero-flux pruning process. If this parameter is `nothing`, no manual pruning is performed in this way.

# Returns
`Vector{String}`

Returns a vector of strings with the reaction ids of all sinks added to the model.
"""
function add_sinks_for_unmatched_metabolites!(
    model::A.AbstractFBCModel,
    metabolite_status_df::DataFrame,
    additive::AbstractString,
    prune_zero_sinks::Union{Vector{String},Nothing};
    sink_opt_outs::Union{Vector{String},Nothing} = nothing,
)
    if isnothing(prune_zero_sinks)
        @info "Add sinks for unmatched metabolites, DO NOT prune sinks automatically"
    else
        @info "Add sinks for unmatched metabolites, automatic pruning of $(length(prune_zero_sinks))"
    end
    if isnothing(sink_opt_outs)
        @info "Add sinks for unmatched metabolites, DO NOT prune sinks manually"
    else
        @info "Add sinks for unmatched metabolites, manual pruning of $(length(prune_zero_sinks))"
    end
    metabolites_with_exchanges = find_metabolites_with_exchanges(model)
    prune_zero_sinks_2 = isnothing(prune_zero_sinks) ? [] : prune_zero_sinks
    not_found_df = @chain metabolite_status_df begin
        @rsubset(:status == "not found", :additive == additive)
        @select(:metabolite)
    end
    added_sink_ids = []
    for metabolite in sort(unique(not_found_df.metabolite))
        sink_up_name = "R_UNKNOWN_SK_UP_$metabolite"
        if !(
            does_manual_prune_list_match_sink_name(sink_up_name, sink_opt_outs) ||
            sink_up_name in prune_zero_sinks_2
        )
            sink_up = Reaction(
                name = sink_up_name,
                stoichiometry = Dict("M_$(metabolite)" => -1.0),
                lower_bound = -1000.0,
                upper_bound = 0.0,
            )
            model.reactions[sink_up_name] = sink_up
            push!(added_sink_ids, sink_up_name)
        else
            # println("Skipping zero flux sink $sink_up_name")
        end
        sink_down_name = "R_UNKNOWN_SK_DOWN_$metabolite"
        if !(
            does_manual_prune_list_match_sink_name(sink_down_name, sink_opt_outs) ||
            sink_down_name in prune_zero_sinks_2
        )
            sink_down = Reaction(
                name = sink_down_name,
                stoichiometry = Dict("M_$(metabolite)" => -1.0),
                lower_bound = 0.0,
                upper_bound = 1000.0,
            )
            model.reactions[sink_down_name] = sink_down
            push!(added_sink_ids, sink_down_name)
        else
            # println("Skipping zero flux sink $sink_down_name")
        end
    end
    return added_sink_ids
end

"""
    list_objectives_in_model(model::A.AbstractFBCModel)

Lists all objectives in `model` to stdout.

# Argument
1. `model::A.AbstractFBCModel`: Model in which to search for ojectives.
"""
function list_objectives_in_model(model::A.AbstractFBCModel)
    obj_coeffs = A.AbstractFBCModels.objective(model)
    rxn_ids = keys(model.reactions)
    for (id, coeff) in zip(rxn_ids, obj_coeffs)
        if !isapprox(coeff, 0.0)
            println("Reaction $id has objective coefficient: $coeff")
        end
    end
end

@doc raw"""
    case_3_constraint_tree!(model::A.AbstractFBCModel, metabolite_status_df::DataFrame, additive::AbstractString)

Sets objective in the model's `ConstraintTree` to prune fluxes according to Case 3 in the Bordbar paper.

``\min \sum_{i=1}^{m} \lvert \Delta x_i \rvert + \sum_{j=1}^{n} \lvert v_j \rvert``

Where ``|\Delta x_i|`` denotes magnitude of the rate of change of the unmeasured metabolites and ``|v_j|`` is the magnitude of the reaction fluxes in the network.

# Arguments
1. `model::A.AbstractFBCModel`: The model in which **the `ConstraintTree` will be mutated**

2. `metabolite_status_df::DataFrame`: DataFrame from [`find_metabolite_matches`](@ref BloodStorageInSilico.UfbaSampler.find_metabolite_matches) to find unmeasured metabolites.

3. `additive::AbstractString`: Additive to search for metabolite measurement availability.

# Returns
`ConstraintTree`

The mutated `ConstraintTree` modified with the objective for Case 3.
"""
function case_3_constraint_tree!(
    model::A.AbstractFBCModel,
    metabolite_status_df::DataFrame,
    additive::AbstractString,
)
    ct = flux_balance_constraints(model)
    flux_ids = collect(keys(ct.fluxes))
    unfound_metabolites =
        @rsubset(metabolite_status_df, :status == "not found", :additive == additive)
    unfound_metabolite_ids = Symbol.(unique(unfound_metabolites.metabolite))

    abs_flux_vars = C.variables(keys = flux_ids, bounds = C.Between(0.0, Inf))
    ct_flux_pos = C.zip(abs_flux_vars, ct.fluxes) do abs_flux_var, flux
        C.Constraint(abs_flux_var.value + flux.value, C.Between(0.0, Inf))
    end
    ct_flux_neg = C.zip(abs_flux_vars, ct.fluxes) do abs_flux_var, flux
        C.Constraint(abs_flux_var.value - flux.value, C.Between(0.0, Inf))
    end

    abs_metabolite_vars =
        C.variables(keys = unfound_metabolite_ids, bounds = C.Between(0.0, Inf))
    ct_metabolite_pos =
        C.zip(abs_metabolite_vars, ct.flux_stoichiometry) do abs_metabolite_var, stoi
            C.Constraint(abs_metabolite_var.value + stoi.value, C.Between(0.0, Inf))
        end
    ct_metabolite_neg =
        C.zip(abs_metabolite_vars, ct.flux_stoichiometry) do abs_metabolite_var, stoi
            C.Constraint(abs_metabolite_var.value - stoi.value, C.Between(0.0, Inf))
        end

    ct =
        ct +
        C.ConstraintTree(:abs_flux_vars => abs_flux_vars) +
        C.ConstraintTree(:flux_pos => ct_flux_pos) +
        C.ConstraintTree(:flux_neg => ct_flux_neg) +
        C.ConstraintTree(:abs_metabolite_vars => abs_metabolite_vars) +
        C.ConstraintTree(:metabolite_pos => ct_metabolite_pos) +
        C.ConstraintTree(:metabolite_neg => ct_metabolite_neg)

    ct.objective = C.Constraint(
        C.sum(a.value for (rxn_id, a) in abs_flux_vars; init = 0.0) +
        C.sum(a.value for (metabolite_id, a) in abs_metabolite_vars; init = 0.0),
    )

    return ct
end

"""
    optimize_case_3(ct::C.ConstraintTree, objective::C.LinearValue)

Create a JuMP model with the given Case 3 `ConstraintTree` and optimize it to find zero flux reactions to prune.

# Arguments
1. `ct::C.ConstraintTree`: `ConstraintTree` with Case 3 objective.

2. `objective::C.LinearValue`: Objective to optimize the constraint tree for. This can be the objective for the `ConstraintTree` passed as the first argument.

# Returns
`C.Tree{Float64}`

`C.Tree{Float64}` with the optimization results substituted in. These results can be used to prune a model.
"""
function optimize_case_3(ct::C.ConstraintTree, objective::C.LinearValue)
    @info "Optimizing case 3"

    num_vars = C.variable_count(ct)
    model = JuMP.Model(HiGHS.Optimizer)
    JuMP.@variable(model, x[1:num_vars])
    JuMP.@objective(model, JuMP.MIN_SENSE, C.substitute(objective, x))

    C.traverse(ct) do c
        b = c.bound
        if b isa C.EqualTo
            JuMP.@constraint(model, C.substitute(c.value, x) == b.equal_to)
        elseif b isa C.Between
            val = C.substitute(c.value, x)
            isinf(b.lower) || JuMP.@constraint(model, val >= b.lower)
            isinf(b.upper) || JuMP.@constraint(model, val <= b.upper)
        end
    end

    JuMP.set_silent(model)
    JuMP.optimize!(model)
    if is_solved_and_feasible(model)
        println("Case 3 optimization success!")
        result_ct = deepcopy(ct)
        var_values = JuMP.value.(model[:x])
        solution_tree = C.substitute_values(result_ct, var_values)
        return solution_tree
    else
        println("OH NO CASE 3 OPTIMIZATION FAILED!")
        return nothing
    end
end

"""
    function analyze_case_3(case_3_optimize_result::C.Tree{Float64})

Analyze the results of the Case 3 optimization to make lists of of sinks added for unmeasured metabolites that have zero flux and non-zero flux. Also gathers these results into a DataFrame for easier manual inspection.

# Argument
1. `case_3_optimize_result::C.Tree{Float64}`: Case 3 optimization result.

# Returns
`Tuple{Vector{String},Vector{String},DataFrame}`

Tuple of reaction ids for zero flux Case 3 sinks, non-zero flux Case 3 sinks, and a status DataFrame for manual inspection.
"""
function analyze_case_3(case_3_optimize_result::C.Tree{Float64})
    zero_case3_sinks = [
        k for (k, v) in case_3_optimize_result.fluxes if
        isapprox(v, 0.0) && contains(string(k), "R_UNKNOWN_SK")
    ]
    nonzero_case3_sinks = [
        k for (k, v) in case_3_optimize_result.fluxes if
        !isapprox(v, 0.0) && contains(string(k), "R_UNKNOWN_SK")
    ]
    sink_status_rows = []
    for zero_case3_sink in zero_case3_sinks
        row = (sink = zero_case3_sink, is_non_zero = false)
        push!(sink_status_rows, row)
    end
    for nonzero_case3_sink in nonzero_case3_sinks
        row = (sink = nonzero_case3_sink, is_non_zero = true)
        push!(sink_status_rows, row)
    end
    sink_status_df = DataFrame(sink_status_rows)
    return zero_case3_sinks, nonzero_case3_sinks, sink_status_df
end

"""
    count_n_all_zero_fluxes(samples_df)

Counts the number of fluxes which have every sample at zero flux and return the ids of reactions that are blocked. This is to assist in finding potentially broken reactions in uFBA jobs.

# Arguments
1. `samples_df`: DataFrame result of an apparently successful sampling run.

# Returns
`Tuple{Int64,Vector{String}}`

Returns a tuple with the count of the fluxes which have all samples at zero as the first element and a vector of blocked reaction_ids as the second element.
"""
function count_n_all_zero_fluxes(samples_df)
    n_all_zero_fluxes = 0
    blocked_reaction_ids = []
    for col_name in names(samples_df)
        col = samples_df[!, col_name]
        n_zeros = sum(isapprox.(col, 0.0, atol = 1.0e-10))
        if n_zeros == length(col)
            n_all_zero_fluxes += 1
            push!(blocked_reaction_ids, col_name)
        end
    end
    return n_all_zero_fluxes, blocked_reaction_ids
end

"""
    execute_ufba_job(job, n_chains = 10)

Execute a uFBA job specified by the first argument with the given number of chains.

# Arguments
1. `job`: `NamedTuple` with the following keys: `additive` to specify the additive solution, `final_time` to specify the time point of the simulation, `pruned_model` to specify the model to optimize, `metabolite_bounds_df` rates of chage of metabolites in a DataFrame.

2. `n_chains`: Number of chains to sample. Defaults to 10.

# Returns
`Tuple{Union{Nothing,DataFrame},Union{Int64,Missing}}`

Returns a tuple of two items. First, a DataFrame of sampled fluxes if successful, or `nothing` is the optimization failed. Second, an Int64 of the number of fluxes that have zeros for all samples or missing if the sampling failed.
"""
function execute_ufba_job(job, n_chains = 10)
    additive = job.additive
    final_time = job.final_time
    pruned_model = job.pruned_model
    metabolite_bounds_df = job.metabolite_bounds_df
    @info "execute_ufba_job: additive: $additive, final_time: $final_time"
    ct = flux_balance_constraints(pruned_model)
    for k in keys(ct.flux_stoichiometry)
        short_metabolite_id = string(k)[3:end]
        bounds = query_metabolite_bounds(
            metabolite_bounds_df,
            additive,
            short_metabolite_id,
            final_time,
        )
        if isnothing(bounds)
            ct.flux_stoichiometry[k].bound = C.EqualTo(0.0)
        else
            lb, ub = bounds
            ct.flux_stoichiometry[k].bound = C.Between(lb, ub)
        end
    end
    objective_flux = flux_balance_analysis(pruned_model; optimizer = HiGHS.Optimizer)
    if isnothing(objective_flux)
        println("OH NO uFBA SIMPLE OPTIMIZATION FAILED!")
        return nothing, missing
    else
        println("Simple optimization succeeded! Sampling fluxes...")
        samples_df = sample_fluxes(pruned_model; n_chains = n_chains)
        n_all_zero_fluxes, blocked_reaction_ids = count_n_all_zero_fluxes(samples_df)
        samples_df[!, :additive] .= additive
        samples_df[!, :final_time] .= final_time
        return samples_df, n_all_zero_fluxes, blocked_reaction_ids
    end
end

"""
    execute_all_ufba_jobs(jobs, n_chains = 10)

Executes and aggregates results from all uFBA jobs specified.

# Arguments
1. `jobs`: Vector of all jobs to be executed.
2. `n_chains`: The number of sampling chains for each job. Defaults to 10.

# Returns
`Tuple{DataFrame,DataFrame,DataFrame,DataFrame}`

A tuple of the following three DataFrames: All sampling results, statuses of each attempted sampling job, counts of statuses across all sampling jobs, and per-model blocked reaction ids.
"""
function execute_all_ufba_jobs(jobs, n_chains = 10)
    all_sampling_dfs_1 = map(jobs) do job
        execute_ufba_job(job, n_chains)
    end
    all_results = [
        (sdf, n_all_zero_fluxes, blocked_reaction_ids) for
        (sdf, n_all_zero_fluxes, blocked_reaction_ids) in all_sampling_dfs_1
    ]
    status_rows = vcat(
        eachrow([
            (
                additive = job.additive,
                final_time = job.final_time,
                status = isnothing(sdf) ? "fail" : "ok",
                n_all_zero_fluxes = n_all_zero_fluxes,
            ) for (job, (sdf, n_all_zero_fluxes, _)) in zip(jobs, all_results)
        ])...,
    )
    blocked_reaction_ids_rows = []
    for (job, (_, _, blocked_reaction_ids)) in zip(jobs, all_results)
        for blocked_reaction_id in blocked_reaction_ids
            blocked_reaction_ids_row = (
                additive = job.additive,
                final_time = job.final_time,
                blocked_reaction_id = blocked_reaction_id,
            )
            push!(blocked_reaction_ids_rows, blocked_reaction_ids_row)
        end
    end
    sampling_df = vcat([sdf for (sdf, _) in all_sampling_dfs_1 if !isnothing(sdf)]...)
    status_df = DataFrame(status_rows)
    status_counts_df = @chain status_df begin
        @groupby(:status)
        combine(nrow => :Count)
    end
    blocked_reaction_ids_df = DataFrame(blocked_reaction_ids_rows)
    return sampling_df, status_df, status_counts_df, blocked_reaction_ids_df
end

"""
    make_ufba_models_for_additives_and_times(metabolite_bounds_df::DataFrame, n_models::Int64)

Create all models that represent each combination of additive and final time point.

# Arguments
1. `metabolite_bounds_df::DataFrame`: The bounds of rates of concentration change for the metabolites.
2. `n_models::Int64`: Number of models to generate. If `-1`, all possible models are created.
3. `exchanges::Union{Nothing,Vector{String}} = nothing`: Passed to `create_fba_model`. If specified, a list of exchanges to add to all uFBA models. If not specified, no exchanges are added to uFBA models.

# Returns
`Vector{NamedTuple}`

Returns a vector of `NamedTuple` with specifications for jobs for each model. Each `NamedTuple` has the following properties:

1. `additive`: Additive
2. `final_time`: Final time point
3. `full_model`: The full model created before pruning
4. `pruned_model`: The model after pruning.
5. `sink_status_df`: DataFrame of status of sinks
6. `metabolite_bounds_df`: Metabolite rate DataFrame used to create the model
7. `zero_case3_sinks`: Sinks that have zero flux that were pruned out
8. `nonzero_case3_sinks`: Sinks that have non-zero flux
9. `added_sink_ids`: Sinks that were added to the model according to the call to [`add_sinks_for_unmatched_metabolites!`](@ref BloodStorageInSilico.UfbaSampler.add_sinks_for_unmatched_metabolites!). More direct than inferring from zero_case3_sinks and non_zero_case3_sinks.
"""
function make_ufba_models_for_additives_and_times(
    metabolite_bounds_df::DataFrame,
    n_models::Int64;
    exchanges::Union{Nothing,Vector{String}} = nothing,
    flux_bounds_overrides_df::Union{Nothing,DataFrame} = nothing,
)
    base_rbc_gem = load_base_rbc_gem()
    final_times = sort(unique(metabolite_bounds_df.final_time))
    additives = sort(unique(metabolite_bounds_df.additive))
    pairs =
        n_models == -1 ? collect(product(additives, final_times)) :
        collect(product(additives, final_times))[1:n_models]
    n_pairs = length(pairs)
    result = map(enumerate(pairs)) do p
        (i, (additive, final_time)) = p
        @info "make_ufba_models_for_additives_and_times: $i of $n_pairs"
        full_model, _ = create_fba_model(
            base_rbc_gem;
            exchanges = exchanges,
            flux_bounds_overrides_df = flux_bounds_overrides_df,
        )
        metabolite_status_df = find_metabolite_matches(
            full_model,
            metabolite_bounds_df,
            additive,
            final_time,
        )
        add_sinks_for_unmatched_metabolites!(
            full_model,
            metabolite_status_df,
            additive,
            nothing,
        )
        ct = case_3_constraint_tree!(full_model, metabolite_status_df, additive)
        case_3_optimize_result_ct = optimize_case_3(ct, ct.objective.value)
        if isnothing(case_3_optimize_result_ct)
            @error "Failed to optimize case 3 for additive: $additive, final_time: $final_time"
        end
        zero_case3_sinks, nonzero_case3_sinks, sink_status_df =
            analyze_case_3(case_3_optimize_result_ct)
        sink_status_df[!, :additive] .= additive
        sink_status_df[!, :final_time] .= final_time
        pruned_model, _ = create_fba_model(
            base_rbc_gem;
            exchanges = exchanges,
            flux_bounds_overrides_df = flux_bounds_overrides_df,
        )
        added_sink_ids = add_sinks_for_unmatched_metabolites!(
            pruned_model,
            metabolite_status_df,
            additive,
            string.(zero_case3_sinks),
        )
        save_ufba_model_sbml(pruned_model, additive, final_time)
        (
            additive = additive,
            final_time = final_time,
            full_model = deepcopy(full_model),
            pruned_model = deepcopy(pruned_model),
            sink_status_df = sink_status_df,
            metabolite_bounds_df = deepcopy(metabolite_bounds_df),
            zero_case3_sinks = zero_case3_sinks,
            nonzero_case3_sinks = nonzero_case3_sinks,
            added_sink_ids = added_sink_ids,
        )
    end
    return result
end

"""
    extract_case3_sinks(ufba_jobs)

Extracts the status of the sinks for unmeasured metabolites for all jobs given and gathers the result into a DataFrame.

# Arguments
1. `ufba_jobs`: The finished ufba_jobs. Each job is a `NamedTuple` with `additive`, `final_time`, and `zero_case3_sinks` properties.

# Returns
`Tuple{DataFrame,DataFrame}`

Returns two DataFrames:
1. Status of unmeasured metabolite sinks for each uFBA job.
2. Aggregated status of unmeasured metabolites sinks for each uFBA job.
"""
function extract_case3_sinks(ufba_jobs)
    status_rows = []
    for ufba_job in ufba_jobs
        for zero_case3_sink in ufba_job.zero_case3_sinks
            row = (
                additive = ufba_job.additive,
                final_time = ufba_job.final_time,
                sink = zero_case3_sink,
                status = "zero",
            )
            push!(status_rows, row)
        end
        for nonzero_case3_sink in ufba_job.nonzero_case3_sinks
            row = (
                additive = ufba_job.additive,
                final_time = ufba_job.final_time,
                sink = nonzero_case3_sink,
                status = "nonzero",
            )
            push!(status_rows, row)
        end
    end
    status_df = DataFrame(status_rows)
    status_aggregated_df = @chain status_df begin
        @groupby(:additive, :final_time, :status)
        combine(nrow => :count)
    end
    return status_df, status_aggregated_df
end

"""
    extract_added_case3_sink_ids(jobs)

Extract and return a DataFrame of the sinks added to each uFBA model from the finished uFBA jobs.

# Arguments
1. `jobs`: The result of the call to [`make_ufba_models_for_additives_and_times`](@ref BloodStorageInSilico.UfbaSampler.make_ufba_models_for_additives_and_times)

# Returns
`DataFrame`

Returns a DataFrame with the following columns:
1. `additive`: The additive
2. `final_time`: Final time of the model
3. `metabolite_id`: The metabolite the sink is for
4. `direction`: up or down depending on the direction of the sink.
5. `added_sink_id`: The reaction id of the corresponding sink.

The DataFrame is sorted by additive, final time. metabolite id, and direction.
"""
function extract_added_case3_sink_ids(jobs)
    rows = []
    for job in jobs
        additive = job.additive
        final_time = job.final_time
        for added_sink_id in job.added_sink_ids
            direction = contains(added_sink_id, "UP") ? "up" : "down"
            metabolite_id =
                replace(added_sink_id, "R_UNKNOWN_SK_UP_" => "", "R_UNKNOWN_SK_DOWN_" => "")
            row = (
                additive = additive,
                final_time = final_time,
                metabolite_id = metabolite_id,
                direction = direction,
                added_sink_id = added_sink_id,
            )
            push!(rows, row)
        end
    end
    unsorted_df = DataFrame(rows)
    sorted_df = @orderby(unsorted_df, :additive, :final_time, :metabolite_id, :direction)
    return sorted_df
end

"""
    sample_fluxes(model; n_chains::Int64, tolerance::Float64)

Sample the allowable flux space of the `model`. Use the `julia -p X...` -p command line option to set the number of workers for this operation.

# Arguments
1. `model`: Model to be sampled.
2. `n_chains::Int64`: The number of chains to calculate, with each chain producing ~126 samples. Defaults to 10 chains.
3. `tolerance::Float64`: The tolerance bounds on the objective.

# Returns
`DataFrame`
1. Returns a `DataFrame` with each reaction as a column and each row a flux sample.
"""
function sample_fluxes(model; n_chains::Int64 = 10, tolerance::Float64 = 0.99)
    println("N Chains: $n_chains")
    s = flux_sample(
        model,
        optimizer = HiGHS.Optimizer,
        objective_bound = relative_tolerance_bound(tolerance),
        n_chains = n_chains,
        workers = workers(),
        collect_iterations = [10],
    )
    s_dict = Dict()
    for reaction_id ∈ keys(s)
        s_dict[reaction_id] = s[reaction_id]
    end
    DataFrame(s_dict)
end

end

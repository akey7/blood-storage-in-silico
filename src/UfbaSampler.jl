module UfbaSampler

using Distributed
using COBREXA, HiGHS, JuMP, MathOptInterface
using Base.Iterators
import ConstraintTrees as C
import SBMLFBCModels as S
import AbstractFBCModels as A
using CSV
using DataFrames
using DataFramesMeta
using ThreadsX
using OrderedCollections
using Chain
using ProgressMeter

include("FbaModelBuilder.jl")
using .FbaModelBuilder
include("PruningOptimizations.jl")
using .PruningOptimizations
include("MetaboliteBounds.jl")
using .MetaboliteBounds

export sample_fluxes,
    ufba_all_additives_all_times,
    histograms_for_reaction_in_additive,
    plot_all_histograms,
    fba,
    is_metabolite_in_exchange,
    list_objectives_in_model,
    display_jump_results,
    make_ufba_models_for_additives_and_times,
    execute_all_ufba_jobs,
    map_reaction_ids_to_reaction_strings,
    extract_sinks,
    init_workers!,
    execute_ufba_job,
    count_n_all_zero_fluxes,
    load_flux_bounds_overrides,
    sbml_add_constant_to_selfclosing_parameters!,
    extract_added_sink_ids,
    decompose_sink_id,
    optimize_constriant_tree

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
                include("src/UfbaSampler.jl")
                using .UfbaSampler
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
`Tuple{Dict{String,String},DataFrame}`

Returns a tuple with two elements:
1. A dictionary mapping reaction ids in the model to a human-readable reaction string and
2. A DataFrame with `:reaction_id` and `:reaction_string` columns.
"""
function map_reaction_ids_to_reaction_strings(model::A.AbstractFBCModel)
    result_dict = OrderedDict()
    reaction_ids = []
    reaction_strings = []
    for rxn_id in sort(string.(keys(model.reactions)))
        push!(reaction_ids, rxn_id)
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
            rxn_string = "$lhs <-- $rhs ($(rxn.lower_bound), $(rxn.upper_bound))"
            result_dict[rxn_id] = rxn_string
            push!(reaction_strings, rxn_string)
        elseif isapprox(rxn.lower_bound, 0.0) && rxn.upper_bound > 0.0
            rxn_string = "$lhs --> $rhs ($(rxn.lower_bound), $(rxn.upper_bound))"
            result_dict[rxn_id] = rxn_string
            push!(reaction_strings, rxn_string)
        else
            rxn_string = "$lhs <-> $rhs ($(rxn.lower_bound), $(rxn.upper_bound))"
            result_dict[rxn_id] = rxn_string
            push!(reaction_strings, rxn_string)
        end
    end
    result_df = DataFrame(reaction_id = reaction_ids, reaction_string = reaction_strings)
    return result_dict, result_df
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
    # println("Wrote $filename")
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
        workers_config = workers()
        samples_df = sample_fluxes(model, workers_config; n_chains = n_chains)
        return solution, samples_df
    end
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
`Tuple{Union{Nothing,DataFrame},Union{Int64,Missing}, Vector}`

Returns a tuple of three items: 
1. First, a DataFrame of sampled fluxes if successful, or `nothing` is the optimization failed. 
2. Second, an Int64 of the number of fluxes that have zeros for all samples or missing if the sampling failed.
3. A vector of blocked reaction ids

If the pruning optimization with [`optimize_case_1`](@ref BloodStorageInSilico.UfbaSampler.PruningOptimizations.optimize_case_1) fails, returns `nothing, missing, []`.
"""
function execute_ufba_job(job, n_chains = 10)
    additive = job.additive
    final_time = job.final_time
    pruned_with_metabolite_bounds_ct = job.pruned_with_metabolite_bounds_ct
    if !isnothing(pruned_with_metabolite_bounds_ct)
        # @info "execute_ufba_job(): additive: $additive, final_time: $final_time"
        # objective_flux =
        #     optimized_values(pruned_with_metabolite_bounds_ct; optimizer = HiGHS.Optimizer)
        objective_value = pruned_with_metabolite_bounds_ct.objective.value
        optimization_status, _ =
            optimize_constriant_tree(pruned_with_metabolite_bounds_ct, objective_value)
        if optimization_status == :fail
            # println("OH NO uFBA SIMPLE OPTIMIZATION FAILED!")
            return nothing, missing, missing
        else
            # println("Simple optimization succeeded! Sampling fluxes...")
            workers_config = workers()
            samples_df = sample_fluxes(
                pruned_with_metabolite_bounds_ct,
                workers_config;
                n_chains = n_chains,
            )
            n_all_zero_fluxes, blocked_reaction_ids = count_n_all_zero_fluxes(samples_df)
            samples_df[!, :additive] .= additive
            samples_df[!, :final_time] .= final_time
            return samples_df, n_all_zero_fluxes, blocked_reaction_ids
        end
    else
        # @error "execute_ufba_job(): optimize_case_1() failed for additive: $additive, final_time: $final_time, skipping"
        return nothing, missing, []
    end
end

"""
    execute_all_ufba_jobs(jobs, rxn_ids_to_strings_df; n_chains = 10)

Executes and aggregates results from all uFBA jobs specified. For jobs returned as the failure case from [`execute_ufba_job`](@ref BloodStorageInSilico.UfbaSampler.execute_ufba_job), creates a row in the statuses of each sampling job DataFrame with `n_all_zero_fluxes` as a `missing` value. Displays a nice green status bar as it goes.

# Arguments
1. `jobs`: Vector of all jobs to be executed.
2. `rxn_ids_to_strings_df`: The DataFrame mapping reactions ids to strings made by [`map_reaction_ids_to_reaction_strings`](@ref BloodStorageInSilico.UfbaSampler.map_reaction_ids_to_reaction_strings)
3. `n_chains = 10`: The number of sampling chains for each job. Defaults to 10.

# Returns
`Tuple{DataFrame,DataFrame,DataFrame,DataFrame}`

A tuple of the following four DataFrames: 
1. All sampling results,
2. Statuses of each attempted sampling job,
3. Counts of statuses across all sampling jobs, and
4. Per-model blocked reaction ids with reaction strings joined in.
"""
function execute_all_ufba_jobs(jobs, rxn_ids_to_strings_df; n_chains = 10)
    n_jobs = length(jobs)
    prog = Progress(n_jobs, "Optimizing and sampling uFBA jobs")
    job_results = map(jobs) do job
        job_result = execute_ufba_job(job, n_chains)
        next!(prog)
        return job_result
    end
    all_results = [
        (sdf, n_all_zero_fluxes, blocked_reaction_ids) for
        (sdf, n_all_zero_fluxes, blocked_reaction_ids) in job_results
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
        if !ismissing(blocked_reaction_ids)
            for blocked_reaction_id in blocked_reaction_ids
                blocked_reaction_ids_row = (
                    additive = job.additive,
                    final_time = job.final_time,
                    blocked_reaction_id = blocked_reaction_id,
                )
                push!(blocked_reaction_ids_rows, blocked_reaction_ids_row)
            end
        end
    end
    sampling_df = vcat([sdf for (sdf, _) in job_results if !isnothing(sdf)]...)
    status_df = DataFrame(status_rows)
    status_counts_df = @chain status_df begin
        @groupby(:status)
        combine(nrow => :Count)
    end
    blocked_reaction_ids_df = DataFrame(blocked_reaction_ids_rows)
    joined_blocked_reaction_ids_df = innerjoin(
        blocked_reaction_ids_df,
        rxn_ids_to_strings_df,
        on = :blocked_reaction_id => :reaction_id,
    )
    return sampling_df, status_df, status_counts_df, joined_blocked_reaction_ids_df
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
5. `metabolite_bounds_df`: Metabolite rate DataFrame used to create the model
6. `zero_sinks`: Sinks that have zero flux that were pruned out
7. `nonzero_sinks`: Sinks that have non-zero flux
8. `added_sink_ids`: Sinks that were added to the model according to the call to [`add_sinks_for_unmatched_metabolites!`](@ref BloodStorageInSilico.UfbaSampler.MetaboliteBounds.add_sinks_for_unmatched_metabolites!). More direct than inferring from zero_sinks and non_zero_sinks.
9. `pruning_method`: The pruning method, currently hardcoded to `:case1`
10. `pruned_with_metabolite_bounds_ct`: A ConstraintTree with metabolite bounds and the pruned set of sinks added, ready for optimziation.
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
    prog = Progress(n_pairs, "Preparing uFBA models")
    result = map(enumerate(pairs)) do p
        (i, (additive, final_time)) = p
        additive_string = String(additive)
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
        first_sink_specifications = (
            metabolite_status_df = metabolite_status_df,
            additive = additive,
            prune_zero_sinks = nothing,
            sink_opt_outs = nothing,
        )
        add_sinks_for_unmatched_metabolites!(full_model, first_sink_specifications)
        case1_ct = flux_balance_constraints(full_model)
        prune_optimize_status, prune_optimize_result =
            optimize_case_1(case1_ct; write_lp_path = nothing)
        if prune_optimize_status == :ok
            case_1_analysis = analyze_case_1_pruning_optimization(prune_optimize_result)
            prune_zero_sinks = string.(case_1_analysis.prune)
            nonzero_sinks = string.(case_1_analysis.keep)
            pruned_model, _ = create_fba_model(
                base_rbc_gem;
                exchanges = exchanges,
                flux_bounds_overrides_df = flux_bounds_overrides_df,
            )
            second_sink_specifications = (
                metabolite_status_df = metabolite_status_df,
                additive = additive,
                prune_zero_sinks = prune_zero_sinks,
                sink_opt_outs = nothing,
            )
            added_sink_ids = add_sinks_for_unmatched_metabolites!(
                pruned_model,
                second_sink_specifications,
            )

            # This SBML will have sinks (if added) but not metabolite bounds.
            # For the graph analysis that is not important at this time.
            save_ufba_model_sbml(pruned_model, additive, final_time)

            pruned_with_metabolite_bounds_ct = flux_balance_constraints(pruned_model)
            add_metabolite_bounds_to_constraint_tree!(
                pruned_with_metabolite_bounds_ct,
                metabolite_bounds_df,
                additive_string,
                final_time,
            )
            next!(prog)
            return (
                additive = additive,
                final_time = final_time,
                full_model = deepcopy(full_model),
                pruned_model = deepcopy(pruned_model),
                metabolite_bounds_df = deepcopy(metabolite_bounds_df),
                zero_sinks = prune_zero_sinks,
                nonzero_sinks = nonzero_sinks,
                added_sink_ids = added_sink_ids,
                pruning_method = :case1,
                pruned_with_metabolite_bounds_ct = pruned_with_metabolite_bounds_ct,
                pruning_breaks_df = nothing,
            )
        else
            pruning_breaks_df = DataFrame(
                additive = additive,
                final_time = final_time,
                broken_case_1_constraint = prune_optimize_result,
            )
            next!(prog)
            return (
                additive = additive,
                final_time = final_time,
                full_model = deepcopy(full_model),
                pruned_model = deepcopy(pruned_model),
                metabolite_bounds_df = deepcopy(metabolite_bounds_df),
                zero_sinks = nothing,
                nonzero_sinks = nothing,
                added_sink_ids = nothing,
                pruning_method = :case1,
                pruned_with_metabolite_bounds_ct = nothing,
                pruning_breaks_df = pruning_breaks_df,
            )
        end
    end
    return result
end

"""
    decompose_sink_id(sink_id)

Extract the direction and metabolite id from a given sink id/name.

# Arguments
1. `sink_id`: The id of the sink.

# Returns
`Tuple{String,String}`

Returns a tuple of metabolite id and direction.
"""
function decompose_sink_id(sink_id)
    sink_str = String(sink_id)
    metabolite_id = @chain sink_str begin
        replace("R_UNKNOWN_SK_DOWN_" => "")
        replace("R_UNKNOWN_SK_UP_" => "")
    end
    direction = occursin(sink_str, "UP") ? "up" : "down"
    return metabolite_id, direction
end

"""
    extract_sinks(ufba_jobs)

Extracts the status of the sinks for unmeasured metabolites for all jobs given and gathers the result into a DataFrame.

# Arguments
1. `ufba_jobs`: The finished ufba_jobs. Each job is a `NamedTuple` with `additive`, `final_time`, and `zero_sinks` properties.

# Returns
`DataFrame`

Returns two DataFrames:
1. Status of unmeasured metabolite sinks for each uFBA job.
"""
function extract_sinks(ufba_jobs)
    status_rows = []
    for ufba_job in ufba_jobs
        pruning_method = ufba_job.pruning_method
        for zero_sink in ufba_job.zero_sinks
            metabolite_id, direction = decompose_sink_id(zero_sink)
            row = (
                pruning_method = pruning_method,
                additive = ufba_job.additive,
                final_time = ufba_job.final_time,
                status = "zero",
                metabolite_id = metabolite_id,
                direction = direction,
                sink = zero_sink,
            )
            push!(status_rows, row)
        end
        for nonzero_sink in ufba_job.nonzero_sinks
            metabolite_id, direction = decompose_sink_id(nonzero_sink)
            row = (
                pruning_method = pruning_method,
                additive = ufba_job.additive,
                final_time = ufba_job.final_time,
                status = "nonzero",
                metabolite_id = metabolite_id,
                direction = direction,
                sink = nonzero_sink,
            )
            push!(status_rows, row)
        end
    end
    status_df = DataFrame(status_rows)
    sorted_df =
        @orderby(status_df, :additive, :final_time, :status, :metabolite_id, :direction)
    return sorted_df
end

"""
    extract_added_sink_ids(jobs)

Extract and return a DataFrame of the sinks added to each uFBA model from the finished uFBA jobs.

# Arguments
1. `jobs`: The result of the call to [`make_ufba_models_for_additives_and_times`](@ref BloodStorageInSilico.UfbaSampler.make_ufba_models_for_additives_and_times)

# Returns
`DataFrame`

Returns a DataFrame with the following columns:
1. `pruning_method`: The pruning method (right now, always `:case1`)
2. `additive`: The additive
3. `final_time`: Final time of the model
4. `metabolite_id`: The metabolite the sink is for
5. `direction`: up or down depending on the direction of the sink.
6. `added_sink_id`: The reaction id of the corresponding sink.

The DataFrame is sorted by additive, final time. metabolite id, and direction.
"""
function extract_added_sink_ids(jobs)
    rows = []
    for job in jobs
        additive = job.additive
        final_time = job.final_time
        pruning_method = job.pruning_method
        for added_sink_id in job.added_sink_ids
            direction = contains(added_sink_id, "UP") ? "up" : "down"
            metabolite_id =
                replace(added_sink_id, "R_UNKNOWN_SK_UP_" => "", "R_UNKNOWN_SK_DOWN_" => "")
            row = (
                pruning_method = pruning_method,
                additive = additive,
                final_time = final_time,
                metabolite_id = metabolite_id,
                direction = direction,
                added_sink_id = added_sink_id,
            )
            push!(rows, row)
        end
    end
    if length(rows) > 0
        unsorted_df = DataFrame(rows)
        sorted_df =
            @orderby(unsorted_df, :additive, :final_time, :metabolite_id, :direction)
        return sorted_df
    else
        empty_df = DataFrame(
            pruning_method = [],
            additive = [],
            final_time = [],
            metabolite_id = [],
            direction = [],
            added_sink_id = [],
        )
        return empty_df
    end
end

function substitute_jump(val::C.LinearValue, vars)
    e = JuMP.AffExpr() # unfortunately @expression(model, 0) is not type stable and gives an Int
    for (i, w) in zip(val.idxs, val.weights)
        if i == 0
            JuMP.add_to_expression!(e, w)
        else
            JuMP.add_to_expression!(e, w, vars[i])
        end
    end
    return e
end

function constraint_jump!(jump_model, expr, b::C.EqualTo; base_name::String)
    JuMP.@constraint(jump_model, expr == b.equal_to, base_name = "$(base_name)_eq")
end

function constraint_jump!(jump_model, expr, b::C.Between; base_name::String)
    isinf(b.lower) ||
        JuMP.@constraint(jump_model, expr >= b.lower, base_name = "$(base_name)_lb")
    isinf(b.upper) ||
        JuMP.@constraint(jump_model, expr <= b.upper, base_name = "$(base_name)_ub")
end

function optimize_constriant_tree(
    ct::C.ConstraintTree,
    objective_value::Union{Nothing,C.Value} = nothing;
    silent::Bool = true,
)
    # Adding functionality to optimization_model() in COBREXA.jl
    ct_paths = []
    jump_model = JuMP.Model(HiGHS.Optimizer)
    JuMP.@variable(jump_model, x[1:C.variable_count(ct)])
    isnothing(objective_value) ||
        JuMP.@objective(jump_model, JuMP.MAX_SENSE, substitute_jump(objective_value, x))
    C.itraverse(ct) do path, con
        ct_path = join(path, ".")
        if ct_path != "objective"
            push!(ct_paths, ct_path)
            isnothing(con.bound) || constraint_jump!(
                jump_model,
                substitute_jump(con.value, x),
                con.bound;
                base_name = ct_path,
            )
        end
    end
    silent && JuMP.set_silent(jump_model)
    JuMP.optimize!(jump_model)
    status = JuMP.termination_status(jump_model)
    if status in [JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL]
        solved_values = JuMP.value.(x)
        solution_tree = C.substitute_values(ct, solved_values)
        return :ok, solution_tree
    elseif status == JuMP.MOI.INFEASIBLE
        conflicted_constraints = optimization_failure_analysis(jump_model)
        return :fail, conflicted_constraints
    else
        return :fail, ["No further information is available."]
    end
end

"""
    sample_fluxes(constraints::C.ConstraintTree, workers_config; n_chains::Int64, tolerance::Float64)

Sample the allowable flux space of the ConstraintTree `ct`. Uses default ACHR sampling method.

Use the `julia -p X...` -p command line option to set the number of workers for this operation.

# Arguments
1. `constraints::C.ConstraintTree`: ConstraintTree to be sampled.
2. `workers_config`: The output of `Distributed.workers()` used to call this function. This sets up the workers.
3. `n_chains::Int64`: The number of chains to calculate. Defaults to 10 chains.
4. `tolerance::Float64`: The tolerance bounds on the objective.

# Returns
`DataFrame`

1. Returns a `DataFrame` with each reaction as a column and each row a flux sample.
"""
function sample_fluxes(
    constraints,
    workers_config;
    n_chains::Int64 = 10,
    tolerance::Float64 = 0.99,
)
    optimizer = HiGHS.Optimizer
    objective = constraints.objective.value
    settings = []
    method = sample_chain_achr
    seed = UInt64(123)
    collect_iterations = [32]

    objective_flux = optimized_values(
        constraints;
        objective = objective,
        output = constraints.objective,
        optimizer = optimizer,
        settings = settings,
    )
    isnothing(objective_flux) && return nothing
    constraints *= :objective_bound^C.Constraint(objective, objective_flux)
    warmup = vcat(
        (
            transpose(v) for (_, vs) in constraints_variability(
                constraints,
                constraints.fluxes;
                optimizer = optimizer,
                settings = settings,
                output = (_, om) -> JuMP.value.(om[:x]),
                output_type = Vector{Float64},
                workers = workers_config,
            ) for v in vs
        )...,
    )

    # I could use kwargs... in the following call but am not using that at
    # this time.

    samples = sample_constraints(
        method,
        constraints;
        seed,
        output = constraints.fluxes,
        start_variables = warmup,
        n_chains,
        collect_iterations,
        workers = workers_config,
    )

    samples_dict = Dict()
    for reaction_id in keys(samples)
        samples_dict[reaction_id] = samples[reaction_id]
    end

    samples_df = DataFrame(samples_dict)
    return samples_df
end

end

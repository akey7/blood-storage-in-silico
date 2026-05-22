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
using Statistics

include("FbaModelBuilder.jl")
using .FbaModelBuilder
include("PruningOptimizations.jl")
using .PruningOptimizations
include("MetaboliteBounds.jl")
using .MetaboliteBounds

export sample_fluxes,
    ufba_all_additives_all_times,
    histograms_for_reaction_in_additive,
    is_metabolite_in_exchange,
    list_objectives_in_model,
    display_jump_results,
    make_ufba_models_for_additives_and_times,
    execute_all_ufba_jobs,
    extract_pruning_overview,
    map_reaction_ids_to_reaction_strings,
    init_workers!,
    execute_ufba_job,
    count_n_all_zero_fluxes,
    load_flux_bounds_overrides,
    sbml_add_constant_to_selfclosing_parameters!,
    decompose_sink_id,
    optimize_constraint_tree,
    extract_broken_constraints,
    extract_unmeasured_relaxations,
    load_reaction_names_and_subsystems,
    load_subsystem_category_map,
    extract_constraint_bounds

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
    load_reaction_names_and_subsystems()

Reads the `input/Reaction Id to Subsystem and Name Map.csv` to find extra information about reactions.

# Returns
`Union{DataFrame,Nothing}`

Returns `nothing` (which will cause errors later) if the file is not found, or the DataFrame contained in that file.
"""
function load_reaction_names_and_subsystems()
    filename = joinpath("input", "Reaction Id to Subsystem and Name Map.csv")
    return isfile(filename) ? CSV.read(filename, DataFrame) : nothing
end

"""
    load_subsystem_category_map()

Reads the subsystem to category map DataFrame. Throws an error and stops if the map file is not found.

# Returns
`DataFrame`

Returns the mapping DataFrame.
"""
function load_subsystem_category_map()
    filename = joinpath("input", "Subsystem Category Map.csv")
    if !isfile(filename)
        error(
            "Can't load $(filename), which is necessary for mapping reaction ids to categories.",
        )
    end
    return CSV.read(filename, DataFrame)
end

function remap_category(original_category)
    if original_category == "Other" ||
       original_category == "Pseudoreactions" ||
       original_category == "Transport reactions"
        return "Other, Transport, & Pseudoreactions"
    else
        return original_category
    end
end

"""
    map_reaction_ids_to_reaction_strings(
        model::A.AbstractFBCModel,
        reaction_names_and_subsystems_df::DataFrame,
        subsystem_category_map_df::DataFrame,
    )

Maps reaction_ids in the given model to human-readable reaction strings specifying reactants and products with an arrow pointing in the direction specified by the bounds of the reaction.

# Arguments
1. `model::A.AbstractFBCModel`: The model to create the reaction strings from.
2. `reaction_names_and_subsystems_df::DataFrame`: DataFrame with `rxn_id`, `subsystem`, and `reaction_name` columns.
3. `subsystem_category_map_df::DataFrame`: DataFrame with `name` (name of the reaction subsystem) and `category` columns.

# Returns
`Tuple{Dict{String,Dict{Symbol,String}},DataFrame}`

Returns a tuple with two elements:
1. A dictionary mapping reaction ids in the model to a human-readable reaction strings, subsystems, and reaction names
2. A DataFrame with `reaction_id`, `reaction_string`, `name`, and `subsystem`, and `category` columns.
"""
function map_reaction_ids_to_reaction_strings(
    model::A.AbstractFBCModel,
    reaction_names_and_subsystems_df::DataFrame,
    subsystem_category_map_df::DataFrame,
)
    subsystem_category_map_dict = Dict()
    for row in eachrow(subsystem_category_map_df)
        subsystem = row.name
        category = row.category
        subsystem_category_map_dict[subsystem] = category
    end
    result_dict = OrderedDict()
    reaction_ids = []
    reaction_strings = []
    reaction_subsystems = []
    reaction_categories = []
    for rxn_id in sort(string.(keys(model.reactions)))
        rxn_df = @rsubset(reaction_names_and_subsystems_df, :rxn_id == String(rxn_id))
        rxn_name = nrow(rxn_df) > 0 ? rxn_df[1, :reaction_name] : "Sink or unknown name"
        rxn_subsystem = nrow(rxn_df) > 0 ? rxn_df[1, :subsystem] : "Other"
        rxn_category = get(subsystem_category_map_dict, rxn_subsystem, "Other")
        push!(reaction_ids, rxn_id)
        push!(reaction_subsystems, rxn_subsystem)
        push!(reaction_categories, rxn_category)
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
            result_dict[rxn_id] = Dict(
                :rxn_string => rxn_string,
                :name => rxn_name,
                :subsystem => rxn_subsystem,
                :category => remap_category(rxn_category),
            )
            push!(reaction_strings, rxn_string)
        elseif isapprox(rxn.lower_bound, 0.0) && rxn.upper_bound > 0.0
            rxn_string = "$lhs --> $rhs ($(rxn.lower_bound), $(rxn.upper_bound))"
            result_dict[rxn_id] = Dict(
                :rxn_string => rxn_string,
                :name => rxn_name,
                :subsystem => rxn_subsystem,
                :category => remap_category(rxn_category),
            )
            push!(reaction_strings, rxn_string)
        else
            rxn_string = "$lhs <-> $rhs ($(rxn.lower_bound), $(rxn.upper_bound))"
            result_dict[rxn_id] = Dict(
                :rxn_string => rxn_string,
                :name => rxn_name,
                :subsystem => rxn_subsystem,
                :category => remap_category(rxn_category),
            )
            push!(reaction_strings, rxn_string)
        end
    end
    rxn_ids_strings_df = DataFrame(
        reaction_id = reaction_ids,
        reaction_string = reaction_strings,
        reaction_subsystem = reaction_subsystems,
        reaction_category = remap_category.(reaction_categories),
    )
    result_df = @chain rxn_ids_strings_df begin
        leftjoin(reaction_names_and_subsystems_df, on = :reaction_id => :rxn_id)
        @orderby(:reaction_id)
    end
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
    prune_status = job.prune_status
    prune_method = job.prune_method
    if prune_status == :ok
        # @info "execute_ufba_job(): additive: $additive, final_time: $final_time"
        objective_value = pruned_with_metabolite_bounds_ct.objective.value
        fba_status, fba_breaks =
            optimize_constraint_tree(pruned_with_metabolite_bounds_ct, objective_value)
        if fba_status == :fail
            # println("OH NO uFBA SIMPLE OPTIMIZATION FAILED!")
            result = (
                samples_df = nothing,
                n_all_zero_fluxes = missing,
                blocked_reaction_ids = missing,
                prune_status = prune_status,
                fba_status = fba_status,
                fba_breaks = fba_breaks,
                job_status = :prune_ok_fba_fail,
            )
            return result
        else
            # println("Simple optimization succeeded! Sampling fluxes...")
            workers_config = workers()
            samples_df, sinks_df = sample_fluxes(
                pruned_with_metabolite_bounds_ct,
                workers_config;
                n_chains = n_chains,
            )
            n_all_zero_fluxes, blocked_reaction_ids = count_n_all_zero_fluxes(samples_df)
            samples_df[!, :additive] .= additive
            samples_df[!, :final_time] .= final_time
            if !isnothing(sinks_df)
                sinks_df[!, :prune_method] .= prune_method
                sinks_df[!, :additive] .= additive
                sinks_df[!, :final_time] .= final_time
                @rtransform!(
                    sinks_df,
                    :metabolite_id = replace(string(:sink_id), "R_REVSK_" => "")
                )
                @select!(
                    sinks_df,
                    :prune_method,
                    :additive,
                    :final_time,
                    :sink_id,
                    :metabolite_id,
                    :median_flux
                )
            end
            result = (
                samples_df = samples_df,
                sinks_df = sinks_df,
                n_all_zero_fluxes = n_all_zero_fluxes,
                blocked_reaction_ids = blocked_reaction_ids,
                prune_status = prune_status,
                fba_status = fba_status,
                fba_breaks = fba_breaks,
                job_status = :ok,
            )
            return result
        end
    else
        # @error "execute_ufba_job(): optimize_case_1() failed for additive: $additive, final_time: $final_time, skipping"
        result = (
            samples_df = nothing,
            sinks_df = nothing,
            n_all_zero_fluxes = missing,
            blocked_reaction_ids = missing,
            prune_status = prune_status,
            fba_status = missing,
            fba_breaks = nothing,
            job_status = :prune_fail_fba_fail,
        )
        return result
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
`NamedTuple`

Returns a named tuple with the following four DataFrames: 
1. `sampling_df`: All sampling results,
2. `status_df`: Statuses of each attempted sampling job,
3. `status_counts_df`: Counts of statuses across all sampling jobs, and
4. `joined_blocked_reaction_ids_df`: Per-model blocked reaction ids with reaction strings joined in.
5. `prune_breaks_df`: DataFrame of constraints conflicted during pruning, consolidated into one DataFrame.
6. `fba_breaks_df`: DataFrame of constraints conflicted during initial FBA optimization, consolidated into one DataFrame.
"""
function execute_all_ufba_jobs(jobs, rxn_ids_to_strings_df; n_chains = 10)
    n_jobs = length(jobs)
    prog = Progress(n_jobs, "Optimizing and sampling uFBA jobs")
    job_results = map(jobs) do job
        job_result = execute_ufba_job(job, n_chains)
        next!(prog)
        return job_result
    end
    status_rows = []
    blocked_reaction_ids_rows = []
    for (job, job_result) in zip(jobs, job_results)
        status_row = (
            prune_method = job.prune_method,
            additive = job.additive,
            final_time = job.final_time,
            job_status = job_result.job_status,
            n_all_zero_fluxes = job_result.n_all_zero_fluxes,
        )
        push!(status_rows, status_row)
        blocked_reaction_ids = job_result.blocked_reaction_ids
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
    status_df = @chain status_rows begin
        DataFrame()
        @orderby(:additive, :final_time)
    end
    blocked_reaction_ids_df = DataFrame(blocked_reaction_ids_rows)
    sampling_dfs = [
        job_result.samples_df for
        job_result in job_results if !isnothing(job_result.samples_df)
    ]
    sinks_dfs = [
        job_result.sinks_df for job_result in job_results if !isnothing(job_result.sinks_df)
    ]
    sampling_df = vcat(sampling_dfs...)
    sinks_df = length(sinks_dfs) > 0 ? vcat(sinks_dfs...) : nothing
    status_counts_df = @chain status_df begin
        @groupby(:job_status)
        combine(nrow => :Count)
        @orderby(:Count)
    end
    broken_constraints = extract_broken_constraints(jobs, job_results)
    joined_blocked_reaction_ids_df =
        join_blocked_reaction_ids(blocked_reaction_ids_df, rxn_ids_to_strings_df)
    result = (
        sampling_df = sampling_df,
        sinks_df = sinks_df,
        status_df = status_df,
        status_counts_df = status_counts_df,
        joined_blocked_reaction_ids_df = joined_blocked_reaction_ids_df,
        prune_breaks_df = broken_constraints.prune_breaks_df,
        fba_breaks_df = broken_constraints.fba_breaks_df,
    )
    return result
end

"""
    join_blocked_reaction_ids(blocked_reaction_ids_df, rxn_ids_to_strings_df)

Joins the blocked reaction ids to their reaction strings. If there are no blocked reactions, returns an empty DataFrame.

# Arguments
1. `blocked_reaction_ids_df`: DataFrame of blocked reactions ids
2. `rxn_ids_to_strings_df`: DataFrame mapping blocked reaction ids to strings.

# Returns
`DataFrame`

Blocked reaction ids joined to their reaction strings.
"""
function join_blocked_reaction_ids(blocked_reaction_ids_df, rxn_ids_to_strings_df)
    if nrow(blocked_reaction_ids_df) > 0
        result = @chain blocked_reaction_ids_df begin
            innerjoin(rxn_ids_to_strings_df, on = :blocked_reaction_id => :reaction_id)
            @orderby(:additive, :final_time, :blocked_reaction_id)
        end
        return result
    else
        return DataFrame(additive = [], final_time = [], blocked_reaction_id = [])
    end
end

"""
    function make_ufba_models_for_additives_and_times(metabolite_bounds_df::DataFrame, n_models::Int64; exchanges::Union{Nothing,Vector{String}} = nothing; flux_bounds_overrides_df::Union{Nothing,DataFrame} = nothing, metabolites_to_ignore::Vector{String} = nothing, prune_method::Symbol = :case3, relax_quantile::Float64 = 0.1, sink_opt_ins::Vector{String})

Create all models that represent each combination of additive and final time point.

# Arguments
1. `metabolite_bounds_df::DataFrame`: The bounds of rates of concentration change for the metabolites.
2. `n_models::Int64`: Number of models to generate. If `-1`, all possible models are created.
3. `exchanges::Union{Nothing,Vector{String}} = nothing`: Passed to `create_fba_model`. If specified, a list of exchanges to add to all uFBA models. If not specified, no exchanges are added to uFBA models.
4. `flux_bounds_overrides_df::Union{Nothing,DataFrame} = nothing`: If specified, a DataFrame of per-reaction flux bounds overrides.
5. `metabolites_to_ignore::Vector{String} = nothing`: If specified, these metabolite bounds are ignored.
6. `prune_method::Symbol = :case3`: Prune method to use. Can be either `:case1` or `:case3`.
7. `relax_strategy::Symbol = :q`: Strategy to find realxation amount. Either `:q` or `:frac_minimum` as noted in [`suggested_unmeasured_metabolite_bounds`](@ref BloodStorageInSilico.UfbaSampler.MetaboliteBounds.suggested_unmeasured_metabolite_bounds).
8. `relax_quantile::Float64 = 0.1`: Relaxation quantile to use. See [`suggested_unmeasured_metabolite_bounds`](@ref BloodStorageInSilico.UfbaSampler.MetaboliteBounds.suggested_unmeasured_metabolite_bounds) for more information.
9. `sink_opt_ins::Vector{String}`: Vector of sinks to create no matter what the pruning results are.
10. `frac_minimum::Float64 = 0.1`: Fraction of minimum measurement for relaxation of unmeasured metabolites as found in [`suggested_unmeasured_metabolite_bounds`](@ref BloodStorageInSilico.UfbaSampler.MetaboliteBounds.suggested_unmeasured_metabolite_bounds) for more information.

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
9. `prune_method`: The pruning method.
10. `pruned_with_metabolite_bounds_ct`: A ConstraintTree with metabolite bounds and the pruned set of sinks added, ready for optimziation.
11. `prune_optimize_status`: Either `:ok` (for a successful prune optimization) or `:fail` for a failed prune optimization.
12. `prune_breaks_df`: If pruning was a `:fail` as indicated by `prune_optimize_status`, this field is populated with a DataFrame reporting the broken constraints. If the pruning was `:ok`, this field is `nothing`.
13. `pruned_default_lb`: Default lower bound for unmeasured metabolites in pruned model.
14. `pruned_default_ub`: Default upper bound for unmeasured metabolites in pruned model.
15. `pruned_unmeasured_metabolites`: Metabolites that were not measured.
16. `pruned_measured_metabolites`: Metabolites that were measured.
17. `relax_quantile`: The quantile of absolute value sof bounds on measured metabolites that was used for unmeasured metabolites.
18. `frac_minimum`: Fraction of minimum measurement used for metabolite relaxation
"""
function make_ufba_models_for_additives_and_times(
    metabolite_bounds_df::DataFrame,
    n_models::Int64;
    exchanges::Union{Nothing,Vector{String}} = nothing,
    flux_bounds_overrides_df::Union{Nothing,DataFrame} = nothing,
    metabolites_to_ignore::Vector{String} = nothing,
    prune_method::Symbol = :case3,
    relax_strategy::Symbol = :q,
    relax_quantile::Float64 = 0.1,
    sink_opt_ins::Vector{String} = nothing,
    frac_minimum::Float64 = 0.1,
)
    base_rbc_gem = load_base_rbc_gem()
    final_times = sort(unique(metabolite_bounds_df.final_time))
    additives = sort(unique(metabolite_bounds_df.additive))
    sink_opt_ins_2 = isnothing(sink_opt_ins) ? String[] : sink_opt_ins
    pairs =
        n_models == -1 ? collect(product(additives, final_times)) :
        collect(product(additives, final_times))[1:n_models]
    n_pairs = length(pairs)
    prog = Progress(n_pairs, "Preparing uFBA models")
    result = map(enumerate(pairs)) do p
        (i, (additive, final_time)) = p
        additive_string = String(additive)
        fba_model_result = create_fba_model(
            base_rbc_gem;
            exchanges = exchanges,
            flux_bounds_overrides_df = flux_bounds_overrides_df,
        )
        full_model = fba_model_result.model
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
            sink_opt_ins = sink_opt_ins,
            metabolites_to_ignore = metabolites_to_ignore,
        )
        add_sinks_for_unmatched_metabolites!(full_model, first_sink_specifications)
        pre_prune_with_metabolite_bounds_ct = flux_balance_constraints(full_model)
        measured_unmeasured = add_metabolite_bounds_to_constraint_tree!(
            pre_prune_with_metabolite_bounds_ct,
            metabolite_bounds_df,
            additive_string,
            final_time;
            metabolites_to_ignore = metabolites_to_ignore,
            relax_strategy = relax_strategy,
            relax_quantile = relax_quantile,
            frac_minimum = frac_minimum,
        )
        unmeasured_metabolite_ids = measured_unmeasured.unmeasured_metabolites
        prune_status, prune_result = optimize_for_pruning(
            prune_method,
            pre_prune_with_metabolite_bounds_ct,
            unmeasured_metabolite_ids,
        )
        if prune_status == :ok
            prune_analysis = analyze_pruning_optimization(prune_result)
            prune_zero_sinks = string.(prune_analysis.prune)
            nonzero_sinks = string.(prune_analysis.keep)
            pruned_fba_model_result = create_fba_model(
                base_rbc_gem;
                exchanges = exchanges,
                flux_bounds_overrides_df = flux_bounds_overrides_df,
            )
            pruned_model = pruned_fba_model_result.model
            second_sink_specifications = (
                metabolite_status_df = metabolite_status_df,
                additive = additive,
                prune_zero_sinks = prune_zero_sinks,
                sink_opt_ins = sink_opt_ins,
                metabolites_to_ignore = metabolites_to_ignore,
            )
            added_sink_ids = add_sinks_for_unmatched_metabolites!(
                pruned_model,
                second_sink_specifications,
            )

            # This SBML will have sinks (if added) but not metabolite bounds.
            # For the graph analysis that is not important at this time.
            save_ufba_model_sbml(pruned_model, additive, final_time)

            pruned_with_metabolite_bounds_ct = flux_balance_constraints(pruned_model)
            pruned_metabolite_bounds_result = add_metabolite_bounds_to_constraint_tree!(
                pruned_with_metabolite_bounds_ct,
                metabolite_bounds_df,
                additive_string,
                final_time;
                metabolites_to_ignore = metabolites_to_ignore,
                relax_strategy = relax_strategy,
                relax_quantile = relax_quantile,
                frac_minimum = frac_minimum,
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
                prune_method = prune_method,
                pruned_with_metabolite_bounds_ct = pruned_with_metabolite_bounds_ct,
                pre_prune_with_metabolite_bounds_ct = pre_prune_with_metabolite_bounds_ct,
                prune_status = prune_status,
                prune_breaks_df = nothing,
                pruned_default_lb = pruned_metabolite_bounds_result.default_lb,
                pruned_default_ub = pruned_metabolite_bounds_result.default_ub,
                pruned_unmeasured_metabolites = pruned_metabolite_bounds_result.unmeasured_metabolites,
                pruned_measured_metabolites = pruned_metabolite_bounds_result.measured_metabolites,
                relax_quantile = relax_strategy == :q ? relax_quantile : missing,
                frac_minimum = relax_strategy == :frac_minimum ? frac_minimum : missing,
            )
        else
            prune_breaks_df = DataFrame(
                prune_method = prune_method,
                additive = additive,
                final_time = final_time,
                broken_constraint = prune_result,
            )
            next!(prog)
            return (
                additive = additive,
                final_time = final_time,
                full_model = deepcopy(full_model),
                pruned_model = nothing,
                metabolite_bounds_df = deepcopy(metabolite_bounds_df),
                zero_sinks = nothing,
                nonzero_sinks = nothing,
                added_sink_ids = nothing,
                prune_method = prune_method,
                pruned_with_metabolite_bounds_ct = nothing,
                pre_prune_with_metabolite_bounds_ct = pre_prune_with_metabolite_bounds_ct,
                prune_status = prune_status,
                prune_breaks_df = prune_breaks_df,
                pruned_default_lb = nothing,
                pruned_default_ub = nothing,
                pruned_unmeasured_metabolites = nothing,
                pruned_measured_metabolites = nothing,
                relax_quantile = relax_strategy == :q ? relax_quantile : missing,
                frac_minimum = relax_strategy == :frac_minimum ? frac_minimum : missing,
            )
        end
    end
    return result
end

"""
    decompose_sink_id(sink_id)

Extract the metabolite id from a given sink id.

# Arguments
1. `sink_id`: The id of the sink.

# Returns
`String`

Returns the metabolite id associated with the sink id.
"""
function decompose_sink_id(sink_id)
    sink_str = String(sink_id)
    metabolite_id = replace(sink_str, "R_REVSK_" => "")
    return metabolite_id
end

"""
    extract_pruning_overview(ufba_jobs)

Extracts the status of sink pruning for metabolites for all jobs given and gathers the result into a DataFrame.

# Arguments
1. `ufba_jobs`: The finished ufba_jobs. Each job is a `NamedTuple` with `additive`, `final_time`, `nonzero_sinks`, and `zero_sinks` properties.

# Returns
`DataFrame`

Status of metabolite sinks for each uFBA job.
"""
function extract_pruning_overview(ufba_jobs)
    status_rows = []
    for ufba_job in ufba_jobs
        prune_method = ufba_job.prune_method
        zero_sinks = ufba_job.zero_sinks
        nonzero_sinks = ufba_job.nonzero_sinks
        if !isnothing(zero_sinks)
            for zero_sink in zero_sinks
                metabolite_id = decompose_sink_id(zero_sink)
                row = (
                    prune_method = prune_method,
                    additive = ufba_job.additive,
                    final_time = ufba_job.final_time,
                    status = "zero",
                    metabolite_id = metabolite_id,
                    sink = zero_sink,
                )
                push!(status_rows, row)
            end
        end
        if !isnothing(nonzero_sinks)
            for nonzero_sink in nonzero_sinks
                metabolite_id = decompose_sink_id(nonzero_sink)
                row = (
                    prune_method = prune_method,
                    additive = ufba_job.additive,
                    final_time = ufba_job.final_time,
                    status = "nonzero",
                    metabolite_id = metabolite_id,
                    sink = nonzero_sink,
                )
                push!(status_rows, row)
            end
        end
    end
    if isempty(status_rows)
        return DataFrame(
            prune_method = [],
            additive = [],
            final_time = [],
            status = [],
            metabolite_id = [],
        )
    else
        status_df = DataFrame(status_rows)
        sorted_df = @orderby(
            status_df,
            :prune_method,
            :additive,
            :final_time,
            :status,
            :metabolite_id
        )
        return sorted_df
    end
end

"""
    extract_unmeasured_relaxations(ufba_jobs)

Extracts the relaxation bounds used for unmeasured metabolites in all models into a DataFrame.

# Arguments
1. `ufba_jobs`: Original uFBA jobs created by [`make_ufba_models_for_additives_and_times`](@ref BloodStorageInSilico.UfbaSampler.make_ufba_models_for_additives_and_times)

# Returns
`DataFrame`

Returns a DataFrame with the following columns:
1. `additive`
2. `final_time`
3. `pruned_default_lb`: Default lower bound. `missing` if the model failed to prune.
4. `pruned_default_ub`: Default upper bound. `missing` if the model failed to prune.
"""
function extract_unmeasured_relaxations(ufba_jobs)
    status_rows = []
    for ufba_job in ufba_jobs
        additive = ufba_job.additive
        final_time = ufba_job.final_time
        pruned_default_lb = ufba_job.pruned_default_lb
        pruned_default_ub = ufba_job.pruned_default_ub
        if !isnothing(pruned_default_lb) && !isnothing(pruned_default_ub)
            status_row = (
                additive = additive,
                final_time = final_time,
                pruned_default_lb = pruned_default_lb,
                pruned_default_ub = pruned_default_ub,
            )
            push!(status_rows, status_row)
        else
            status_row = (
                additive = additive,
                final_time = final_time,
                pruned_default_lb = missing,
                pruned_default_ub = missing,
            )
            push!(status_rows, status_row)
        end
    end
    result_df = @chain status_rows begin
        DataFrame()
        @orderby(:additive, :final_time)
    end
    return result_df
end

function ct_to_rows!(
    metabolite_rows,
    flux_rows,
    ct,
    additive,
    final_time,
    prune_status,
    measured_metabolites,
)
    C.itraverse(ct) do path, con
        path_str = string.(path)
        branch = first(path_str)
        leaf = last(path_str)
        constraint_path = join(string.(path), ".")
        b = con.bound
        isnothing(b) && return
        bound_equal_to = b isa C.EqualTo ? b.equal_to : missing
        bound_ub = b isa C.Between ? b.upper : missing
        bound_lb = b isa C.Between ? b.lower : missing
        constraint_type = b isa C.EqualTo ? "equality" : "between"
        if branch == "fluxes"
            flux_row = (
                additive = additive,
                final_time = final_time,
                prune_status = prune_status,
                reaction_id = leaf,
                constraint_path = constraint_path,
                constraint_type = constraint_type,
                lb = bound_lb,
                ub = bound_ub,
                equal_to = bound_equal_to,
            )
            push!(flux_rows, flux_row)
        else
            is_measured = leaf in measured_metabolites ? "measured" : "unmeasured"
            metabolite_row = (
                additive = additive,
                final_time = final_time,
                prune_status = prune_status,
                metabolite_id = leaf,
                is_measured = is_measured,
                constraint_path = constraint_path,
                constraint_type = constraint_type,
                lb = bound_lb,
                ub = bound_ub,
                equal_to = bound_equal_to,
            )
            push!(metabolite_rows, metabolite_row)
        end
    end
end

function extract_constraint_bounds(ufba_jobs)
    metabolite_rows = []
    flux_rows = []
    for ufba_job in ufba_jobs
        additive = ufba_job.additive
        final_time = ufba_job.final_time
        measured_metabolites = string.(ufba_job.pruned_measured_metabolites)
        pruned_with_metabolite_bounds_ct = ufba_job.pruned_with_metabolite_bounds_ct
        prune_status = ufba_job.prune_status
        pre_prune_with_metabolite_bounds_ct = ufba_job.pre_prune_with_metabolite_bounds_ct
        ct =
            prune_status == :ok ? pruned_with_metabolite_bounds_ct :
            pre_prune_with_metabolite_bounds_ct
        ct_to_rows!(
            metabolite_rows,
            flux_rows,
            ct,
            additive,
            final_time,
            prune_status,
            measured_metabolites,
        )
    end
    metabolites_df = @chain metabolite_rows begin
        DataFrame()
        @orderby(:additive, :final_time, :is_measured, :metabolite_id)
    end
    fluxes_df = @chain flux_rows begin
        DataFrame()
        @orderby(:additive, :final_time, :reaction_id)
    end
    result = (metabolites_df = metabolites_df, fluxes_df = fluxes_df)
    return result
end

"""
    extract_broken_constraints(jobs, job_results)

Called by [`execute_all_ufba_jobs`](@ref BloodStorageInSilico.UfbaSampler.execute_all_ufba_jobs) to find all broken pruning and simple FBA optimization constraints during execution of all uFBA jobs.

TODO: Revisit in the future if I need to capture relaxation values that are attempted in a failed prune job.

# Arguments
1. `jobs`: Original uFBA jobs created by [`make_ufba_models_for_additives_and_times`](@ref BloodStorageInSilico.UfbaSampler.make_ufba_models_for_additives_and_times)
2. `job_results`: Results of executed uFBA jobs from a local variable in [`execute_all_ufba_jobs`](@ref BloodStorageInSilico.UfbaSampler.execute_all_ufba_jobs)

# Returns
`NamedTuple`

Returns a named tuple with the following fields
1. `prune_breaks_df`: DataFrame with constraints broken in the pruning optimizations. Has columns `additive`, `final_time`, `prune_break`.
2. `fba_breaks_df`: DataFrame with constraints broken in the simple FBA optimizations. Has columns `additive`, `final_time`, `fba_break`.
"""
function extract_broken_constraints(jobs, job_results)
    prune_breaks_dfs = []
    fba_breaks_dfs = []
    for (job, job_result) in zip(jobs, job_results)
        additive = job.additive
        final_time = job.final_time
        prune_status = job.prune_status
        fba_status = job_result.fba_status
        if prune_status == :fail
            single_prune_breaks_df = job.prune_breaks_df
            push!(prune_breaks_dfs, single_prune_breaks_df)
        end
        if !ismissing(fba_status) && fba_status == :fail
            fba_breaks = job_result.fba_breaks
            if !isnothing(fba_breaks)
                single_fba_break_df = DataFrame(
                    additive = additive,
                    final_time = final_time,
                    fba_break = fba_breaks,
                )
                push!(fba_breaks_dfs, single_fba_break_df)
            end
        end
    end
    unsorted_prune_breaks_df =
        length(prune_breaks_dfs) > 0 ? vcat(prune_breaks_dfs...) :
        DataFrame(additive = [], final_time = [], broken_constraint = [])
    prune_breaks_df =
        @orderby(unsorted_prune_breaks_df, :additive, :final_time, :broken_constraint)
    unsorted_fba_breaks_df =
        length(fba_breaks_dfs) > 0 ? vcat(fba_breaks_dfs...) :
        DataFrame(additive = [], final_time = [], fba_break = [])
    fba_breaks_df = @orderby(unsorted_fba_breaks_df, :additive, :final_time, :fba_break)
    result = (prune_breaks_df = prune_breaks_df, fba_breaks_df = fba_breaks_df)
    return result
end

"""
    substitute_jump(val::C.LinearValue, vars)

Called by [`optimize_constraint_tree`](@ref BloodStorageInSilico.UfbaSampler.optimize_constraint_tree) to create JuMP models.

Copied from `substitute_jump` in COBREXA.jl. Used to assemble a `C.LinearValue` into a `JuMP.AffExpr` to create JuMP constraints and objective from `ConstraintTrees`

# Arguments
1. `val::C.LinearValue`: The `LinearValue` from which to construct the expression.
2. `vars`: The JuMP variable(s) used to create the expression.

# Returns
`JuMP.AffExpr`

The `AffrExpr` for incorporation in the JuMP model.
"""
function substitute_jump(val::C.LinearValue, vars)
    e = JuMP.AffExpr()
    for (i, w) in zip(val.idxs, val.weights)
        if i == 0
            JuMP.add_to_expression!(e, w)
        else
            JuMP.add_to_expression!(e, w, vars[i])
        end
    end
    return e
end

"""
    constraint_jump!(jump_model, expr, b::C.EqualTo; base_name::String)

Called by [`optimize_constraint_tree`](@ref BloodStorageInSilico.UfbaSampler.optimize_constraint_tree) to create JuMP models.

Mutates the given JuMP model to add a constriant from an expression and given `C.EqualTo` bound. Also accepts a string used to name the objective (the suffix `_eq` is added to this string).

# Arguments
1. `jump_model`: The JuMP model to mutate.
2. `expr`: Expression of the constraint.
3. `b::C.EqualTo`: The bound of the constraint.
4. `base_name::String`: The base name of the constraint. `_eq` is appended to the base name.
"""
function constraint_jump!(jump_model, expr, b::C.EqualTo; base_name::String)
    JuMP.@constraint(jump_model, expr == b.equal_to, base_name = "$(base_name)_eq")
end

"""
    constraint_jump!(jump_model, expr, b::C.EqualTo; base_name::String)

Called by [`optimize_constraint_tree`](@ref BloodStorageInSilico.UfbaSampler.optimize_constraint_tree) to create JuMP models.

Mutates the given JuMP model to add a constriant from an expression and given `C.Between` bound. Also accepts a string used to name the objective (the suffixes `_lb` and `_ub` are added to this string).

# Arguments
1. `jump_model`: The JuMP model to mutate.
2. `expr`: Expression of the constraint.
3. `b::C.Between`: The bound of the constraint.
4. `base_name::String`: The base name of the constraint. `_lb` or `_ub` is appended to the base name depending on the side of the constraint.
"""
function constraint_jump!(jump_model, expr, b::C.Between; base_name::String)
    isinf(b.lower) ||
        JuMP.@constraint(jump_model, expr >= b.lower, base_name = "$(base_name)_lb")
    isinf(b.upper) ||
        JuMP.@constraint(jump_model, expr <= b.upper, base_name = "$(base_name)_ub")
end

"""
    optimize_constraint_tree(ct::C.ConstraintTree, objective_value::Union{Nothing,C.Value} = nothing; silent::Bool = true)

Creates a JuMP model from the given `ConstraintTree` and optimizes it. Draws inspiration from `optimization_model()` in COBREXA.jl and adds functionality to name constraints and debug broken models using those constraint names upon optimization failure.

# Arguments
1. `ct::C.ConstraintTree`: Tree containing the constraints.
2. `objective_value::Union{Nothing,C.Value} = nothing`: Objective value to use in the optimization. If in doubt, using the `ConstraintTree` to be optimized, pass `ct.objective.value`
3. `silent::Bool = true`: If `true`, silences JuMP during optimization to clean up script outputs.

# Returns
`Tuple{Symbol,Union{C.ConstraintTree,Vector{String}}}`

Returns a tuple with two elements.
1. `:ok` or `:fail`: The status of the optimization.
2. `C.ConstraintTree` of the optimized values OR `Vector{String}` of conflicted constraints that broke the optimization
"""
function optimize_constraint_tree(
    ct::C.ConstraintTree,
    objective_value::Union{Nothing,C.Value} = nothing;
    silent::Bool = true,
)
    jump_model = JuMP.Model(HiGHS.Optimizer)
    JuMP.@variable(jump_model, x[1:C.variable_count(ct)])
    isnothing(objective_value) ||
        JuMP.@objective(jump_model, JuMP.MAX_SENSE, substitute_jump(objective_value, x))
    C.itraverse(ct) do path, con
        ct_path = join(path, ".")
        if ct_path != "objective"
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
`Tuple{DataFrame,Union{Nothing,DataFrame}}`

Returns a tuple with two elements:
1. A `DataFrame` with each non-sink reaction as a column and each row a flux sample.
2. `nothing` or a DataFrame with sink reaction ids in one column and median flux in another column.
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
    # this time as I find kwargs to make the code a confusing mess.

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
        if !occursin("R_REVSK_", string(reaction_id))
            samples_dict[reaction_id] = samples[reaction_id]
        end
    end
    samples_df = DataFrame(samples_dict)

    sinks_rows = [
        (sink_id = sink_id, median_flux = median(samples[sink_id])) for
        sink_id in keys(samples) if occursin("R_REVSK_", string(sink_id))
    ]
    if length(sinks_rows) > 0
        sinks_df = DataFrame(sinks_rows)
        return samples_df, sinks_df
    else
        return samples_df, nothing
    end
end

end

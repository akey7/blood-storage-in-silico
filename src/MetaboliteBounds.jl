module MetaboliteBounds

using CSV
using DataFrames
using DataFramesMeta
import AbstractFBCModels as A
import ConstraintTrees as C
using COBREXA
import AbstractFBCModels: stoichiometry
import AbstractFBCModels.CanonicalModel: Model, Reaction, Metabolite, Gene, Coupling
using Statistics

export load_metabolite_bounds,
    query_metabolite_bounds,
    find_metabolite_matches,
    add_metabolite_bounds_to_constraint_tree!,
    print_metabolite_bounds_on_constraint_tree,
    add_sinks_for_unmatched_metabolites!,
    find_metabolites_with_exchanges,
    does_sink_id_match_list,
    load_metabolite_measurement_opt_outs,
    load_sink_opt_ins

function load_sink_opt_ins()
    sink_opt_ins_filename = joinpath("input", "sink_opt_ins.csv")
    opt_ins_df =
        isfile(sink_opt_ins_filename) ? CSV.read(sink_opt_ins_filename, DataFrame) :
        DataFrame(metabolite_id_of_enabled_sink = [])
    result = String.(sort(unique(opt_ins_df.metabolite_id_of_enabled_sink)))
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
    load_metabolite_measurement_opt_outs()

Load the metabolite opt-out list from `input/metabolite_measurement_opt_outs.csv`. This list is used to ignore absolute quant estimations in model construction.

# Returns
`Vector{String}`

Vector of strings of metabolite ids (without the leading `M_`) to opt out of.
"""
function load_metabolite_measurement_opt_outs()
    filename = joinpath("input", "metabolite_measurement_opt_outs.csv")
    if isfile(filename)
        df = CSV.read(filename, DataFrame)
        return String.(sort(unique(df.disabled_metabolite_id)))
    else
        return String[]
    end
end

"""
    query_metabolite_bounds(metabolite_bounds_df, additive, metabolite, final_time)

Find the rate of concentration chage for the metabolite in the given additive at the given final time. Returns `nothing` if not found.

# Arguments
1. `metabolite_bounds_df`: DataFrame as loaded by [`load_metabolite_bounds`](@ref BloodStorageInSilico.UfbaSampler.MetaboliteBounds.load_metabolite_bounds).
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
    does_manual_prune_list_match_sink_name(sink_id::String, matching_sink_ids::Union{Vector{String},Nothing} = nothing)

Determines if the given sink name contains any of the substrings in the given sink opt-outs list.

# Arguments
1. `sink_id::String`: The name of the sink.
2. `matching_sink_ids::Union{Vector{String},Nothing} = nothing`: If specified, contains a list of substrings that are matched against the given sink name.

# Returns
`Bool`

Returns `true` if one of the provided substrings matches the given sink name. Returns `false` if the substring list is not provided or none of the substrings are found
"""
function does_sink_id_match_list(
    sink_id::String,
    matching_sink_ids::Union{Vector{String},Nothing} = nothing,
)
    if isnothing(matching_sink_ids)
        return false
    else
        for matching_sink_id in matching_sink_ids
            if contains(sink_id, matching_sink_id)
                return true
            end
        end
        return false
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
    # @info "Matching metabolites, additive: $additive, final_time: $final_time"
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
                lb = missing,
                ub = missing,
            )
            push!(status_rows, status_row)
            not_found_count += 1
        elseif is_metabolite_in_exchange(model, short_metabolite_id)
            lb, ub = bounds
            status_row = (
                additive = additive,
                metabolite = short_metabolite_id,
                status = "in exchange",
                lb = lb,
                ub = ub,
            )
            in_exchange_count += 1
            if isapprox(lb, 0.0) && isapprox(ub, 0.0)
                @warn "$additive $short_metabolite_id is fixed at 0.0"
            end
        else
            lb, ub = bounds
            status_row = (
                additive = additive,
                metabolite = short_metabolite_id,
                status = "found",
                lb = lb,
                ub = ub,
            )
            push!(status_rows, status_row)
            found_count += 1
            if isapprox(lb, 0.0) && isapprox(ub, 0.0)
                @warn "$additive $short_metabolite_id is fixed at 0.0"
            end
        end
    end
    metabolite_status_df = DataFrame(status_rows)
    # println(
    #     "Found $found_count, in exchange $in_exchange_count, not found $not_found_count",
    # )
    return metabolite_status_df
end

"""
    suggested_unmeasured_metabolite_bounds(metabolite_bounds_df::DataFrame, additive::String, final_time::Int64; p::Float64 = 0.5)

Suggest upper and lower bounds for unmeasured metabolites for appropriate model relaxation. It does this by looking at the absolute values of the lower and upper bounds and finding the given percentile within that vector.

# Arguments
1. `metabolite_bounds_df::DataFrame`: Bounds of measured metabolites.
2. `additive::String`: Additive to search within the metabolite bounds.
3. `final_time::Int64`: Final time to search within the metabolite bounds.
4. `strategy::Symbol = :q`: If `:q`, looks for the value of the quantile noted in `p` parameter. If `:tenth_minimum`, 0.1x of the minimum absolute value of measure metabolite abundances.
5. `p::Float64 = 0.1`: Percentile of the measured absolute values to base the bounds off of. Defaults to searching for the median.

# Returns
`Tuple{Float64,Float64}`

Suggested lower and upper bounds for unmeasured metabolites.
"""
function suggested_unmeasured_metabolite_bounds(
    metabolite_bounds_df::DataFrame,
    additive::String,
    final_time::Int64;
    p::Float64 = 0.1,
    strategy::Symbol = :q,
)
    selection_df = @rsubset(
        metabolite_bounds_df,
        :additive == additive,
        :final_time == final_time,
        !isapprox(:lb, 0.0),
        !isapprox(:ub, 0.0)
    )
    if strategy == :tenth_minimum
        min_lb_abs = minimum(abs.(selection_df.lb))
        min_ub_abs = minimum(abs.(selection_df.ub))
        min_abs = minimum([min_lb_abs, min_ub_abs])
        overall = 0.1 * min_abs
        return -overall, overall
    elseif strategy == :q
        abs_bounds = []
        for lb in selection_df.lb
            push!(abs_bounds, abs(lb))
        end
        for ub in selection_df.ub
            push!(abs_bounds, abs(ub))
        end
        quantile_value = quantile(abs_bounds, p)
        return -quantile_value, quantile_value
    else
        throw(ArgumentError("Unknown strategy $strategy"))
    end
end

"""
    add_metabolite_bounds_to_constraint_tree!(ct::C.ConstraintTree, metabolite_bounds_df::DataFrame, additive::String, final_time::Int64; metabolites_to_ignore::Union{Vector{String},Nothing} = nothing, relax_quantile::Float64 = 0.5)

Adds dx/dt metabolite rate of change bounds to the given ConstraintTree. The constraint tree should come from `flux_balance_constraints()`. The bounds are created by replacing `C.EqualTo(0.0)` constraints on the `:flux_stoichiometry` branch with `C.Between(lb, ub)` constraints. Unmeasured metabolites have upper and lower bounds set to percentile measurement suggested by [`suggested_unmeasured_metabolite_bounds`](@ref BloodStorageInSilico.UfbaSampler.MetaboliteBounds.suggested_unmeasured_metabolite_bounds)

Mutates the given ConstraintTree in place.

# Arguments
1. `ct::C.ConstraintTree`: ConstraintTree to modify
2. `metabolite_bounds_df::DataFrame`: DataFrame with the upper and lower bounds of metabolite concentration dx/dt.
3. `additive::String`: Additive to find in the bounds DataFrame
4. `final_time::Int64`: Final time to find in the DataFrame.
5. `metabolites_to_ignore::Union{Vector{String},Nothing} = nothing`: If `nothing`, incorporates constraints for all metabolites in the DataFrame. If specified, ignores the metabolites specified (omit the leading `M_` in this list).
6. `relax_strategy::Symbol = :q`: Strategy to find realxation amount. Either `:q` or `:tenth_minimum` as noted in [`suggested_unmeasured_metabolite_bounds`](@ref BloodStorageInSilico.UfbaSampler.MetaboliteBounds.suggested_unmeasured_metabolite_bounds).
7. `relax_quantile::Float64 = 0.5`: The percentile of the metabolite measurements to set upper and lower bounds of unmeasured to. If unspecified, defaults to 0.5.

# Returns
`NamedTuple`

Returns a named tuple with the following fields:
1. `unmeasured_metabolites`: Vector of symbols of metabolites that did not have measurements that were incorporated into the constraint tree.
2. `measured_metabolites`: Vector of symbols of metabolites that have measurements that were incorporated into the constraint tree.
3. `default_ub`: The default upper bound of unmeasured metabolites.
4. `default_lb`: The default lower bound of unmeasured metabolites.
"""
function add_metabolite_bounds_to_constraint_tree!(
    ct::C.ConstraintTree,
    metabolite_bounds_df::DataFrame,
    additive::String,
    final_time::Int64;
    metabolites_to_ignore::Union{Vector{String},Nothing} = nothing,
    relax_strategy::Symbol = :q,
    relax_quantile::Float64 = 0.1,
)
    metabolites_to_ignore_2 = !isnothing(metabolites_to_ignore) ? metabolites_to_ignore : []
    default_lb, default_ub = suggested_unmeasured_metabolite_bounds(
        metabolite_bounds_df,
        additive,
        final_time;
        strategy = relax_strategy,
        p = relax_quantile,
    )
    unmeasured_metabolites = Symbol[]
    measured_metabolites = Symbol[]
    # n_metabolites_to_ignore_2 = length(metabolites_to_ignore_2)
    # @info "add_metabolite_bounds_to_constraint_tree!(): Ignoring $n_metabolites_to_ignore_2 metabolites"
    for k in keys(ct.flux_stoichiometry)
        short_metabolite_id = string(k)[3:end]
        if short_metabolite_id ∉ metabolites_to_ignore_2
            bounds = query_metabolite_bounds(
                metabolite_bounds_df,
                additive,
                short_metabolite_id,
                final_time,
            )
            if isnothing(bounds)
                # ct.flux_stoichiometry[k].bound = C.EqualTo(0.0)
                ct.flux_stoichiometry[k].bound = C.Between(default_lb, default_ub)
                push!(unmeasured_metabolites, k)
            else
                lb, ub = bounds
                ct.flux_stoichiometry[k].bound = C.Between(lb, ub)
                push!(measured_metabolites, k)
            end
        else
            # println("Skipping bounds for metabolite id $short_metabolite_id")
            continue
        end
    end
    result = (
        unmeasured_metabolites = unmeasured_metabolites,
        measured_metabolites = measured_metabolites,
        default_lb = default_lb,
        default_ub = default_ub,
    )
    return result
end

"""
    find_metabolites_with_exchanges(model::A.AbstractFBCModel)

Finds extracellular metabolites with exchanges in the provided model and returns a list of the metabolite ids found. Used by [`add_sinks_for_unmatched_metabolites!`](@ref BloodStorageInSilico.UfbaSampler.MetaboliteBounds.add_sinks_for_unmatched_metabolites!).

# Arguments
1. `model::A.AbstractFBCModel`: The model which has the metabolites and exchanges of interest.

# Returns
`Vector{String}`

Returns a list of metabolites with exchanges.
"""
function find_metabolites_with_exchanges(model::A.AbstractFBCModel)
    exchange_ids = [
        reaction_id for (reaction_id, _) in model.reactions if contains(reaction_id, "R_EX")
    ]
    metabolite_ids = [replace(exchange_id, "R_EX_" => "") for exchange_id in exchange_ids]
    return metabolite_ids
end

function should_sink_id_be_included(
    sink_id;
    unfound_metabolite_ids,
    sink_opt_ins,
    prune_zero_sinks,
)
    if does_sink_id_match_list(sink_id, sink_opt_ins)
        return true
    elseif does_sink_id_match_list(sink_id, prune_zero_sinks)
        return false
    elseif does_sink_id_match_list(sink_id, unfound_metabolite_ids)
        return true
    else
        return false
    end
end

"""
    add_sinks_for_unmatched_metabolites!(model::A.AbstractFBCModel, NamedTuple)

Add sinks for unmeasured (umatched) metabolites in the model UNLESS those metabolites are already part of an exchange. Exchanges take precedence, see [`find_metabolites_with_exchanges`](@ref BloodStorageInSilico.UfbaSampler.MetaboliteBounds.find_metabolites_with_exchanges) for details. This method mutates the given model in place.

# Arguments
1. `model::A.AbstractFBCModel`: Model to add sinks to. **This model is mutated in place.**
2. `NamedTuple`: A named tuple with additional data to use while adding sinks.

The named tuple needs the following elements
1. `metabolite_status_df`: The metabolite status DataFrame that specifies which metabolites have measurements and therefore do not need sinks.
2. `additive`: The additive to search for measurements in.
3. `prune_zero_sinks`: The vector of sinks to remove as determined by analyzing the Case 1 optimization. If `nothing`, no sinks are removed from this process.
4. `sink_opt_ins`: The manually defind vector of metabolite ids that must have sinks.

The order of precedence for opt-ins, opt-out, and pruning is:
1. `sink_opt_ins`: If a sink is opted-in, this takes first priority.
2. `prune_zero_sinks`: If a sink is pruned algorithmically, this takes last priority.

# Returns
`Vector{String}`

Returns a vector of strings with the reaction ids of all sinks finally added to the model after processing the sink specifications.
"""
function add_sinks_for_unmatched_metabolites!(
    model::A.AbstractFBCModel,
    sink_specifications::NamedTuple,
)
    metabolite_status_df = sink_specifications.metabolite_status_df
    additive = sink_specifications.additive
    prune_zero_sinks = sink_specifications.prune_zero_sinks
    sink_opt_ins = sink_specifications.sink_opt_ins
    metabolites_with_exchanges = find_metabolites_with_exchanges(model)
    prune_zero_sinks_2 = isnothing(prune_zero_sinks) ? String[] : string.(prune_zero_sinks)
    not_found_df = @chain metabolite_status_df begin
        @rsubset(:status == "not found", :additive == additive)
        @select(:metabolite)
    end
    unfound_metabolite_ids = String.(sort(unique(not_found_df.metabolite)))
    all_metabolite_ids = String.(sort(unique(metabolite_status_df.metabolite)))

    # for metabolite_id in unfound_metabolite_ids
    #     if metabolite_id in metabolites_with_exchanges
    #         # println("Skipping sinks for $metabolite_id which has an exchange.")
    #         continue
    #     end
    #     sink_name = "R_REVSK_$metabolite_id"
    #     if does_manual_prune_list_match_sink_name(sink_name, sink_opt_ins) ||
    #        sink_name ∉ prune_zero_sinks_2
    #         sink = Reaction(
    #             name = sink_name,
    #             stoichiometry = Dict("M_$(metabolite_id)" => -1.0),
    #             lower_bound = -1000.0,
    #             upper_bound = 1000.0,
    #         )
    #         model.reactions[sink_name] = sink
    #         push!(added_sink_ids, sink_name)
    #     else
    #         # println("Skipping zero flux sink $sink_name")
    #     end
    # end

    sink_ids_to_add = [
        (metabolite_id, "R_REVSK_$metabolite_id") for
        metabolite_id in all_metabolite_ids if should_sink_id_be_included(
            "R_REVSK_$metabolite_id";
            unfound_metabolite_ids = unfound_metabolite_ids,
            sink_opt_ins = sink_opt_ins,
            prune_zero_sinks = prune_zero_sinks_2,
        )
    ]

    final_added_sink_ids = []
    for (metabolite_id, sink_id) in sink_ids_to_add
        sink = Reaction(
            name = sink_id,
            stoichiometry = Dict("M_$(metabolite_id)" => -1.0),
            lower_bound = -1000.0,
            upper_bound = 1000.0,
        )
        model.reactions[sink_id] = sink
        push!(final_added_sink_ids, sink_id)
    end
    return final_added_sink_ids
end

"""
    print_metabolite_bounds_on_constraint_tree(ct::C.ConstraintTree)

This function is for debugging. Prints metabolite bounds on the given ConstraintTree `:flux_stoichiometry` branch to ensure they were added.

# Arguments
1. `ct::C.ConstraintTree`: ConstraintTree to print values from
"""
function print_metabolite_bounds_on_constraint_tree(ct::C.ConstraintTree)
    function walk(tree, path = "")
        for (key, node) in pairs(tree)
            current_path = isempty(path) ? string(key) : "$path.$key"
            if node isa C.Constraint
                println("Symbol: $current_path | Bounds: $(node.bound)")
            elseif node isa C.ConstraintTree
                walk(node, current_path)
            end
        end
    end
    walk(ct.flux_stoichiometry)
end

end

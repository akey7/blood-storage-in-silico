module MetaboliteBounds

using CSV
using DataFrames
using DataFramesMeta
import AbstractFBCModels as A
import ConstraintTrees as C
using COBREXA

export load_metabolite_bounds,
    query_metabolite_bounds,
    find_metabolite_matches,
    add_metabolite_bounds_to_constraint_tree!,
    print_metabolite_bounds_on_constraint_tree

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
    println(
        "Found $found_count, in exchange $in_exchange_count, not found $not_found_count",
    )
    return metabolite_status_df
end

function add_metabolite_bounds_to_constraint_tree!(
    ct::C.ConstraintTree,
    metabolite_bounds_df::DataFrame,
    additive::String,
    final_time::Int64;
    metabolites_to_ignore::Union{Vector{String},Nothing} = nothing,
)
    metabolites_to_ignore_2 = !isnothing(metabolites_to_ignore) ? metabolites_to_ignore : []
    for k in keys(ct.flux_stoichiometry)
        short_metabolite_id = string(k)[3:end]
        bounds = query_metabolite_bounds(
            metabolite_bounds_df,
            additive,
            short_metabolite_id,
            final_time,
        )
        if isnothing(bounds) || short_metabolite_id ∉ metabolites_to_ignore_2
            ct.flux_stoichiometry[k].bound = C.EqualTo(0.0)
        else
            lb, ub = bounds
            ct.flux_stoichiometry[k].bound = C.Between(lb, ub)
        end
    end

    # Just return something, even though this was modified in place.
    return ct
end

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

module PruningOptimizations

using DataFrames
using DataFramesMeta
using JuMP
using COBREXA
using HiGHS
import AbstractFBCModels as A
import ConstraintTrees as C

export case_3_constraint_tree!,
    optimize_case_3, case_1_constraint_tree!, optimize_case_1, analyze_pruning_optimization

@doc raw"""
    case_3_constraint_tree!(model::A.AbstractFBCModel, metabolite_status_df::DataFrame, additive::AbstractString)

Creates an objective and associated using the the model's `ConstraintTree` to prune fluxes according to Case 3 in the Bordbar (2016) paper.

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
2. `objective::C.LinearValue`: Objective to optimize the constraint tree for. This can be the objective for the `ConstraintTree` passed as the first argument, , and accessed as `ct.objective.value` at invocation time.

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

# A custom ConstratintTrees bound struct used for indicator variables.
# From: https://cobrexa.github.io/ConstraintTrees.jl/stable/3-mixed-integer-optimization/
# See also https://jump.dev/JuMP.jl/stable/tutorials/linear/sudoku/#Mixed-integer-linear-programming-formulation
mutable struct IntegerFromTo <: C.Bound
    from::Int
    to::Int
end

@doc raw"""
    case_1_constraint_tree(model::A.AbstractFBCModel)

Creates an objective and associated using the model's `ConstraintTree` to prune fluxes according to Case 1 in the Bordbar (2016) paper.

``\min \sum_{i=1}^{m} 1_{\Delta x_i \neq 0}``

# Arguments
1. `model::A.AbstractFBCModel`: Model with the `ConstraintTree`

# Returns
`ConstraintTree`

Returns the modified `ConstraintTree` with the proper objective for optimization.
"""
function case_1_constraint_tree!(model::A.AbstractFBCModel)
    ct = flux_balance_constraints(model)
    sink_ids =
        [flux_id for (flux_id, _) in ct.fluxes if occursin("R_UNKNOWN_SK", string(flux_id))]
    indicator_bounds = [IntegerFromTo(0, 1) for _ in eachindex(sink_ids)]
    indicator_variables =
        :indicators^C.variables(keys = sink_ids, bounds = indicator_bounds)
    new_ct = ct + indicator_variables
    new_ct.objective = C.Constraint(
        C.sum(new_ct.indicators[Symbol(sink_id)].value for sink_id in sink_ids; init = 0.0),
    )
    return new_ct
end

"""
    jump_constraint(m, x, v::C.Value, b::C.EqualTo)

Taken from [Example: Mixed integer optimization (MILP)](https://cobrexa.github.io/ConstraintTrees.jl/stable/3-mixed-integer-optimization/#Example:-Mixed-integer-optimization-(MILP))

Sets an `EqualTo` constraint in a JuMP model. Part of a multi-dispatch function with 3 methods.

# Arguments
1. `m`: JuMP model
2. `x`: Reference to variable on which constraint will be set.
3. `v::C.Value`: Value to set the constraint to
4. `b::C.EqualTo`: The `C.EqualTo` bound

# Returns
The specified JuMP constraint
"""
function jump_constraint(m, x, v::C.Value, b::C.EqualTo)
    JuMP.@constraint(m, C.substitute(v, x) == b.equal_to)
end

"""
    jump_constraint(m, x, v::C.Value, b::C.Between)

Taken from [Example: Mixed integer optimization (MILP)](https://cobrexa.github.io/ConstraintTrees.jl/stable/3-mixed-integer-optimization/#Example:-Mixed-integer-optimization-(MILP))

Sets an `Between` constraint in a JuMP model. Part of a multi-dispatch function with 3 methods.

# Arguments
1. `m`: JuMP model
2. `x`: Reference to variable on which constraint will be set.
3. `v::C.Value`: Value to set the constraint to
4. `b::C.Between`: The `C.Between` bound

# Returns
The specified JuMP constraint
"""
function jump_constraint(m, x, v::C.Value, b::C.Between)
    isinf(b.lower) || JuMP.@constraint(m, C.substitute(v, x) >= b.lower)
    isinf(b.upper) || JuMP.@constraint(m, C.substitute(v, x) <= b.upper)
end

"""
    jump_constraint(m, x, v::C.Value, b::IntegerFromTo)

Taken from [Example: Mixed integer optimization (MILP)](https://cobrexa.github.io/ConstraintTrees.jl/stable/3-mixed-integer-optimization/#Example:-Mixed-integer-optimization-(MILP)). See also [Mixed integer linear programming formulation](https://jump.dev/JuMP.jl/stable/tutorials/linear/sudoku/#Mixed-integer-linear-programming-formulation)

Sets an `IntegerFromTo` constraint in a JuMP model. Part of a multi-dispatch function with 3 methods.

# Arguments
1. `m`: JuMP model
2. `x`: Reference to variable on which constraint will be set.
3. `v::C.Value`: Value to set the constraint to
4. `b::C.IntegerFromTo`: The `C.Between` bound

# Returns
The specified JuMP constraint.
"""
function jump_constraint(m, x, v::C.Value, b::IntegerFromTo)
    # var = JuMP.@variable(m, binary = true)  # Appears to generate same results as integer = true setup
    var = JuMP.@variable(m, integer = true)
    JuMP.@constraint(m, var >= b.from)
    JuMP.@constraint(m, var <= b.to)
    JuMP.@constraint(m, C.substitute(v, x) == var)
end

"""
    optimize_case_1(ct::C.ConstraintTree, objective::C.Value)

Create a JuMP model with the given Case 1 `ConstraintTree` and optimize it to find zero flux reactions to prune.

# Arguments
1. `ct::C.ConstraintTree`: `ConstraintTree` with Case 3 objective.
2. `objective::C.LinearValue`: Objective to optimize the constraint tree for. This can be the objective for the `ConstraintTree` passed as the first argument, and accessed as `ct.objective.value` at invocation time.

# Returns
`C.Tree{Float64}`

`C.Tree{Float64}` with the optimization results substituted in. These results can be used to prune a model.
"""
function optimize_case_1(ct::C.ConstraintTree, objective::C.Value)
    jump_model = JuMP.Model(HiGHS.Optimizer)
    JuMP.@variable(jump_model, x[1:C.variable_count(ct)])
    JuMP.@objective(jump_model, JuMP.MIN_SENSE, C.substitute(objective, x))
    C.traverse(ct) do c
        isnothing(c.bound) || jump_constraint(jump_model, x, c.value, c.bound)
    end
    JuMP.set_silent(jump_model)
    JuMP.optimize!(jump_model)
    if is_solved_and_feasible(jump_model)
        println("Case 1 optimization success!")
        result_ct = deepcopy(ct)
        var_values = JuMP.value.(jump_model[:x])
        solution_tree = C.substitute_values(result_ct, var_values)
        return solution_tree
    else
        println("OH NO CASE 1 OPTIMIZATION FAILED!")
        return nothing
    end
end

"""
    analyze_pruning_optimization(pruning_optimization_result::C.Tree{Float64})

Analyze the results of the Case 3 or Case 1 optimization to make lists of of sinks added for unmeasured metabolites that have zero flux and non-zero flux. Also gathers these results into a DataFrame for easier manual inspection.

# Argument
1. `pruning_optimization_result::C.Tree{Float64}`: Case 3 optimization result.

# Returns
`Tuple{Vector{String},Vector{String},DataFrame}`

Tuple of reaction ids for zero flux Case 3 sinks, non-zero flux Case 3 sinks, and a status DataFrame for manual inspection.
"""
function analyze_pruning_optimization(pruning_optimization_result::C.Tree{Float64})
    zero_sinks = [
        k for (k, v) in pruning_optimization_result.fluxes if
        isapprox(v, 0.0) && contains(string(k), "R_UNKNOWN_SK")
    ]
    nonzero_sinks = [
        k for (k, v) in pruning_optimization_result.fluxes if
        !isapprox(v, 0.0) && contains(string(k), "R_UNKNOWN_SK")
    ]
    sink_status_rows = []
    for zero_sink in zero_sinks
        row = (sink = zero_sink, is_non_zero = false)
        push!(sink_status_rows, row)
    end
    for nonzero_sink in nonzero_sinks
        row = (sink = nonzero_sink, is_non_zero = true)
        push!(sink_status_rows, row)
    end
    unordered_df = DataFrame(sink_status_rows)
    sink_status_df = @orderby(unordered_df, :is_non_zero, :sink)
    return zero_sinks, nonzero_sinks, sink_status_df
end

end

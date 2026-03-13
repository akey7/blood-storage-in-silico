module PruningOptimizations

using DataFrames
using DataFramesMeta
using JuMP
using COBREXA
using HiGHS
import AbstractFBCModels as A
import ConstraintTrees as C

export case_3_constraint_tree,
    optimize_case_3,
    case_1_constraint_tree,
    optimize_case_1,
    analyze_pruning_optimization,
    check_case_1_optimization_results,
    list_non_zeros,
    inspect_results

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
function case_3_constraint_tree(
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
1. `model::A.AbstractFBCModel`: Model to extract the `ConstraintTree` from.

# Returns
`ConstraintTree`

Returns a new `ConstraintTree` derived from the given model with the proper objective for optimization.
"""
function case_1_constraint_tree(model::A.AbstractFBCModel)
    ct = flux_balance_constraints(model)
    sink_ids = [id for (id, _) in ct.fluxes if occursin("R_UNKNOWN_SK", string(id))]
    indicator_vars =
        :indicators^C.variables(
            keys = [Symbol("ind_", id) for id in sink_ids],
            bounds = [IntegerFromTo(0, 1) for _ in sink_ids],
        )
    full_ct = ct + indicator_vars
    BIG_M = 1000.0
    couplings = C.ConstraintTree()
    for id in sink_ids
        ind_id = Symbol("ind_", id)
        v = full_ct.fluxes[id].value
        z = full_ct.indicators[ind_id].value
        couplings[Symbol("up_", id)] = C.Constraint(v - BIG_M * z, (-BIG_M, 0.0))
        couplings[Symbol("lo_", id)] = C.Constraint(v + BIG_M * z, (0.0, BIG_M))
    end
    final_ct = full_ct + :coupling^couplings

    # # Force the first sink to always be on so that at least one indicator is 1
    # final_ct.fluxes[sink_ids[1]].bound = C.Between(0.1, 1000.0)

    final_ct.objective = C.Constraint(
        sum(full_ct.indicators[Symbol("ind_", id)].value for id in sink_ids),
        nothing, # No bound, this is an objective
    )

    for (id, flux) in final_ct.fluxes
        if flux.bound == (-Inf, Inf)
            println("Cleanup: $id was unbounded, setting to finite bounds")
            flux.bound = C.Between(-1000.0, 1000.0)
        end
    end

    return final_ct
end

"""
    optimize_case_1(ct::C.ConstraintTree, objective::C.Value)

Create a JuMP model with the given Case 1 `ConstraintTree` and optimize it to find zero flux reactions to prune. If the model is infeasible, attempts to list constraints that make the model infeasible. To help with potential debugging, writes the JuMP model diagnostics to `output/debug_model.lp`

This function contains multiple inner functions to help with translating `ConstraintTree` constraints to JuMP constraints:
1. `register_var!()`: Makes a new JuMP variable
2. `find_all_indices!()`: Safely find indices by iterating keys only.
3. `process_tree!()`: Traverse `ConstraintTree` to find variable definitions and apply bounds
4. `to_jump()`: Helper to convert `C.Value` to JuMP `AffExpr`
5. `add_constraints!()`: Add Constraints recursively from the `ConstraintTree`

# Arguments
1. `ct::C.ConstraintTree`: `ConstraintTree` with Case 1 objective.
2. `objective::C.LinearValue`: Objective to optimize the constraint tree for. This can be the objective for the `ConstraintTree` passed as the first argument, and accessed as `ct.objective.value` at invocation time.

# Returns
`C.Tree{Float64}`

`C.Tree{Float64}` with the optimization results substituted in. These results can then be used to prune a model.
"""
function optimize_case_1(ct::C.ConstraintTree, objective::C.Value)
    jump_model = JuMP.Model(HiGHS.Optimizer)
    JuMP.set_optimizer_attribute(jump_model, "mip_feasibility_tolerance", 1e-8)
    JuMP.set_optimizer_attribute(jump_model, "primal_feasibility_tolerance", 1e-8)
    jump_vars = Dict{Int,JuMP.VariableRef}()

    function register_var!(idx::Int)
        if !haskey(jump_vars, idx)
            jump_vars[idx] = JuMP.@variable(jump_model)
        end
    end

    function find_all_indices!(val)
        val isa C.Value || return
        for idx in keys(val.idxs)
            register_var!(idx)
        end
    end

    function process_tree!(subtree)
        subtree isa C.ConstraintTree || return
        MAX_BOUND = 1000.0
        for (name, entry) in subtree
            if entry isa C.Constraint
                find_all_indices!(entry.value)
            elseif entry isa C.ConstraintTree
                process_tree!(entry)
            elseif hasproperty(entry, :index) && hasproperty(entry, :bound)
                idx = entry.index
                register_var!(idx)
                v = jump_vars[idx]
                bound = entry.bound
                if bound isa IntegerFromTo
                    JuMP.set_lower_bound(v, Float64(bound.lower))
                    JuMP.set_upper_bound(v, Float64(bound.upper))
                    # JuMP.set_integer(v)
                    JuMP.set_binary(v)
                elseif bound isa C.Between
                    bound_lower = isinf(bound.lower) ? -MAX_BOUND : bound.lower
                    bound_upper = isinf(bound.upper) ? MAX_BOUND : bound.upper
                    JuMP.set_lower_bound(v, bound_lower)
                    JuMP.set_upper_bound(v, bound_upper)
                else
                    bound_type = typeof(bound)
                    @error "Case 1 optimization: unknown bound type $bound_type for $v, stopping"
                end
            end
        end
    end

    find_all_indices!(objective)
    process_tree!(ct)

    function to_jump(val::C.Value)
        expr = JuMP.AffExpr(0.0)
        for idx in keys(val.idxs)
            coeff = val.idxs[idx]
            JuMP.add_to_expression!(expr, coeff, jump_vars[idx])
        end
        return expr
    end

    jump_constraints = Dict{String,JuMP.ConstraintRef}()

    function add_constraints!(subtree, prefix = "")
        subtree isa C.ConstraintTree || return
        for (name, entry) in subtree
            full_name = isempty(prefix) ? string(name) : "$(prefix).$(name)"
            if entry isa C.Constraint && string(name) != "objective"
                expr = to_jump(entry.value)
                b = entry.bound

                if b isa C.Between
                    jump_constraints[full_name] = JuMP.@constraint(
                        jump_model,
                        b.lower <= expr <= b.upper,
                        base_name=full_name
                    )
                elseif b isa Float64
                    jump_constraints[full_name] =
                        JuMP.@constraint(jump_model, expr == b, base_name=full_name)
                elseif b isa C.EqualTo
                    jump_constraints[full_name] = JuMP.@constraint(
                        jump_model,
                        expr == b.equal_to,
                        base_name=full_name
                    )
                elseif b isa IntegerFromTo
                    jump_constraints[full_name] = JuMP.@constraint(
                        jump_model,
                        b.from <= expr <= b.to,
                        base_name=full_name
                    )
                else
                    constraint_type = typeof(b)
                    @error "Case 1 optimization: unknown constraint type $constraint_type for $full_name"
                end
            elseif entry isa C.ConstraintTree
                add_constraints!(entry, full_name)
            end
        end
    end
    add_constraints!(ct)
    JuMP.@objective(jump_model, JuMP.MIN_SENSE, to_jump(objective))
    jump_model_filename = joinpath("output", "debug_model.lp")
    write_to_file(jump_model, jump_model_filename)

    # JuMP.optimize!(jump_model)
    # status = JuMP.termination_status(jump_model)
    # if status in [JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL]
    #     values_dict = Dict(idx => JuMP.value(v) for (idx, v) in jump_vars)
    #     return C.substitute_values(ct, values_dict)
    # elseif status == JuMP.MOI.DUAL_INFEASIBLE
    #     println("--- Model is $status ---")
    #     println("There are some values in the model, so the model likely has something unbounded. Here is what we know")
    #     for (tree_idx, jump_var_ref) in jump_vars
    #         val = JuMP.value(jump_var_ref)
    #         if abs(val) > 1e-6
    #             println("  Tree Index [$tree_idx]: $val")
    #         end
    #     end
    #     error("Optimization failed: Model is infeasible.")
    # elseif status == JuMP.MOI.INFEASIBLE
    #     println("--- Model is Infeasible. Starting Conflict Analysis ---")
    #     JuMP.compute_conflict!(jump_model)
    #     println("The following constraints contribute to the conflict:")
    #     for (name, con) in jump_constraints
    #         if JuMP.get_attribute(con, JuMP.MOI.ConstraintConflictStatus()) ==
    #            JuMP.MOI.IN_CONFLICT
    #             println(" - $con")
    #         end
    #     end
    #     error("Optimization failed: Model is infeasible.")
    # else
    #     error(
    #         "Optimization failed with the following status and no further information is available: $status",
    #     )
    # end

    try
        JuMP.optimize!(jump_model)
        status = JuMP.termination_status(jump_model)
        if JuMP.has_values(jump_model)
            max_idx = maximum(keys(jump_vars))
            values_vector = zeros(Float64, max_idx)
            for (idx, v) in jump_vars
                values_vector[idx] = JuMP.value(v)
            end
            results = C.substitute_values(ct, values_vector)
            if status != JuMP.MOI.OPTIMAL
                @warn "Solver finished with non-optimal status: $status. Returning partial results."
            end
            return results
        else
            @error "Solver finished with status $status but no values were returned."
            if status == JuMP.MOI.INFEASIBLE
                JuMP.compute_conflict!(jump_model)
                println("The following constraints contribute to the conflict:")
                for (name, con) in jump_constraints
                    if JuMP.get_attribute(con, JuMP.MOI.ConstraintConflictStatus()) ==
                    JuMP.MOI.IN_CONFLICT
                        println(" - $con")
                    end
                end
            end
            return nothing
        end
    catch e
        println("\n!!! Optimization or Substitution Crashed !!!")
        println("Error type: ", typeof(e))
        if JuMP.has_values(jump_model)
            println("Emergency Value Dump")
            vars = JuMP.all_variables(jump_model)
            vals = value.(vars)
            for v in vals
                println(v)
            end
        end
        rethrow(e)
    end
end

function inspect_results(tree, prefix="", threshold=1e-6)
    # Check if the current node is a leaf (Float64)
    if tree isa Float64
        if abs(tree) > threshold
            println("$prefix: $tree")
        end
        return
    end
    for (name, subtree) in tree
        new_prefix = isempty(prefix) ? string(name) : "$prefix.$name"
        inspect_results(subtree, new_prefix, threshold)
    end
end

"""
    analyze_pruning_optimization(pruning_optimization_result::C.Tree{Float64})

Analyze the results of the Case 3 or Case 1 optimization to make lists of of sinks added for unmeasured metabolites that have zero flux and non-zero flux. Also gathers these results into a DataFrame for easier manual inspection.

# Argument
1. `pruning_optimization_result::C.Tree{Float64}`: Case 3 optimization result.

# Returns
`Tuple{Vector{String},Vector{String},DataFrame}`

Tuple of reaction ids for zero flux sinks, non-zero flux sinks, and a status DataFrame for manual inspection.
"""
function analyze_pruning_optimization(pruning_optimization_result::C.Tree{Float64})
    pruning_data = pruning_optimization_result.fluxes
    zero_sinks = [
        k for (k, v) in pruning_data if
        isapprox(v, 0.0) && contains(string(k), "R_UNKNOWN_SK")
    ]
    nonzero_sinks = [
        k for (k, v) in pruning_data if
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
    n_nonzero_sinks = length(nonzero_sinks)
    n_zero_sinks = length(zero_sinks)
    if n_nonzero_sinks < 1 && n_zero_sinks > 0
        @warn "Empty nonzero sinks, $n_zero_sinks zero sinks"
    end
    unordered_df = DataFrame(sink_status_rows)
    sink_status_df = @orderby(unordered_df, :is_non_zero, :sink)
    return zero_sinks, nonzero_sinks, sink_status_df
end

function check_case_1_optimization_results(ct::C.Tree{Float64})
    nonzero_indicator_ids = []
    zero_indicator_ids = []
    function walk(tree, path = "")
        for (key, node) in pairs(tree)
            current_path = isempty(path) ? string(key) : "$path.$key"
            if node isa Float64
                if !isapprox(node, 0.0)
                    push!(nonzero_indicator_ids, current_path)
                else
                    push!(zero_indicator_ids, current_path)
                end
            elseif node isa C.ConstraintTree
                walk(node, current_path)
            end
        end
    end
    walk(ct.indicators)
    n_zero_indicators = length(zero_indicator_ids)
    if length(nonzero_indicator_ids) < 1
        @warn "Did not find any non-zero indicators, but found $n_zero_indicators zero indicators."
    end
    return nonzero_indicator_ids
end

end

module PruningOptimizations

using DataFrames
using DataFramesMeta
using JuMP
using COBREXA
using HiGHS
import AbstractFBCModels as A
import ConstraintTrees as C
using Printf

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

"""
    jump_constraint(m, x, v::C.Value, b::C.EqualTo)

Attach a ConstraintTrees equality bound to a JuMP model.
"""
function jump_constraint(m, x, v::C.Value, b::C.EqualTo)
    @constraint(m, C.substitute(v, x) == b.equal_to)
end

"""
    jump_constraint(m, x, v::C.Value, b::C.Between)

Attach a ConstraintTrees interval bound to a JuMP model.
"""
function jump_constraint(m, x, v::C.Value, b::C.Between)
    isinf(b.lower) || @constraint(m, C.substitute(v, x) >= b.lower)
    isinf(b.upper) || @constraint(m, C.substitute(v, x) <= b.upper)
end

"""
    bound_big_m(bound; fallback = 1000.0)

Choose a big-M from a sink bound when possible, otherwise use `fallback`.
"""
function bound_big_m(bound; fallback::Float64 = 1000.0)
    if bound isa C.Between
        vals = Float64[]
        isinf(bound.lower) || push!(vals, abs(bound.lower))
        isinf(bound.upper) || push!(vals, abs(bound.upper))
        return isempty(vals) ? fallback : max(maximum(vals), 1e-9)
    elseif bound isa C.EqualTo
        return max(abs(bound.equal_to), 1e-9)
    else
        return fallback
    end
end

"""
    optimize_case_1(
        ct::A.ConstraintTree;
        optimizer = HiGHS.Optimizer,
        fallback_M::Float64 = 1000.0,
        force_first_sink_on::Bool = false,
        force_first_sink_lb::Float64 = 0.1,
        silent::Bool = true,
        write_lp_path::Union{Nothing,String} = "output/debug_case1.lp",
    )

JuMP MILP for Bordbar (2016) Case 1:

``\min \sum_{i=1}^{m} 1_{\Delta x_i \neq 0}``

a sum of binary indicators, with each indicator `i` determines whether sink reaction `i` is allowed to carry flux.

# Arguments
1. `ct::C.ConstraintTree`: ConstraintTree with sinks and dx/dt metabolites bounds added.
2. `optimizer = HiGHS.Optimizer`: Reference to an optimizer.
3. `fallback_M::Float64 = 1000.0`: If a sink bound is not found when creating indocator/coupling constriants, this is the fallback value.
4. `force_first_sink_on::Bool = false`: A debugging option. If `true`, forcibly sets the first sink to have non-zero flux, which will force the corresponding indicator to 1. Defaults to `false`, which does not force any sinks, and which should be used for general sinnk pruning.
5. `force_first_sink_lb::Float64 = 0.1`: A non-zero lower bound to force the first sink on with if `force_first_sink_on` is `true`.
6. `silent::Bool = true`: If `true`, the optimizer output is silenced.
7. `write_lp_path::Union{Nothing,String} = "output/debug_case1.lp"`: A filename to write the JuMP model to for debugging. If `nothing`, does not write the debugging file.

# Returns

Returns a named tuple with:
1. `solution_tree`: base ConstraintTree with continuous variables substituted
2. `indicator_values`: `Dict{Symbol,Float64}` mapping sink id => binary value
3. `sink_ids`: Sink ids
4. `jump_model`: JuMP model
"""
function optimize_case_1(
    ct::C.ConstraintTree;
    optimizer = HiGHS.Optimizer,
    fallback_M::Float64 = 1000.0,
    force_first_sink_on::Bool = false,
    force_first_sink_lb::Float64 = 0.1,
    silent::Bool = true,
    write_lp_path::Union{Nothing,String} = "output/debug_case1.lp",
)
    sink_ids = [id for (id, _) in ct.fluxes if occursin("R_UNKNOWN_SK", string(id))]
    isempty(sink_ids) && error("No sink reactions matching `R_UNKNOWN_SK` were found.")
    jump_model = JuMP.Model(optimizer)
    silent && JuMP.set_silent(jump_model)
    x = Vector{JuMP.VariableRef}(undef, C.variable_count(ct))
    for i in eachindex(x)
        x[i] = @variable(jump_model, base_name = "x_$i")
    end
    C.traverse(ct) do c
        isnothing(c.bound) || jump_constraint(jump_model, x, c.value, c.bound)
    end
    @variable(jump_model, z[sink_ids], Bin)

    # Sink (vi) indicator (zi) coupling
    #
    #       -M_i * z_i <= v_i <= M_i * z_i
    #
    # If z_i = 0, then v_i = 0.
    # If z_i = 1, then v_i is allowed within ±M_i.
    #
    for id in sink_ids
        v_expr = C.substitute(ct.fluxes[id].value, x)
        M_i = bound_big_m(ct.fluxes[id].bound; fallback = fallback_M)

        @constraint(jump_model, v_expr <=  M_i * z[id])
        @constraint(jump_model, v_expr >= -M_i * z[id])
    end
    if force_first_sink_on
        forced_id = sink_ids[1]
        forced_v = C.substitute(ct.fluxes[forced_id].value, x)
        @constraint(jump_model, forced_v >= force_first_sink_lb)
    end
    @objective(jump_model, Min, sum(z[id] for id in sink_ids))
    if !isnothing(write_lp_path)
        mkpath(dirname(write_lp_path))
        write_to_file(jump_model, write_lp_path)
    end
    JuMP.optimize!(jump_model)
    status = JuMP.termination_status(jump_model)
    if !(status in (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL))
        error("Optimization failed with termination status: $status")
    end
    solved_values = JuMP.value.(x)
    solution_tree = C.substitute_values(ct, solved_values)
    indicator_values = Dict(id => JuMP.value(z[id]) for id in sink_ids)

    # Begin diagnostics
    @printf("objective = %.12f\n", JuMP.objective_value(jump_model))
    for id in sink_ids
        v = solution_tree.fluxes[id]
        zi = indicator_values[id]
        if !(isapprox(v, 0.0) && isapprox(zi, 0.0))
            @printf("%s   flux = %.12f   indicator = %.12f\n", string(id), v, zi)
        end
    end

    if force_first_sink_on
        forced_id = sink_ids[1]
        @printf(
            "FORCED %s   flux = %.12f   indicator = %.12f\n",
            string(forced_id),
            solution_tree.fluxes[forced_id],
            indicator_values[forced_id],
        )
    end
    # End diagnostics

    return (
        solution_tree = solution_tree,
        indicator_values = indicator_values,
        sink_ids = sink_ids,
        jump_model = jump_model,
    )
end

# If needed, optimizer debugging code.
# status = JuMP.termination_status(jump_model)
#     if status in [JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL]
#         solved_values = JuMP.value.(jump_model[:x])
#         solution_tree = C.substitute_values(cs, solved_values)

#         # Diagnostics
#         @printf("objective = %.12f\n", JuMP.objective_value(jump_model))
#         for id in sink_ids
#             z = solution_tree.indicators[Symbol("ind_", id)]
#             v = solution_tree.fluxes[id]
#             if !(isapprox(v, 0.0) && isapprox(z, 0.0))
#                 @printf("%s   flux = %.12f   indicator = %.12f\n", string(id), v, z)
#             end
#         end
#         forced_id = sink_ids[1]
#         forced_ind = Symbol("ind_", forced_id)
#         @printf(
#             "FORCED %s   flux = %.12f   indicator = %.12f\n",
#             string(forced_id),
#             solution_tree.fluxes[forced_id],
#             solution_tree.indicators[forced_ind],
#         )
#         # End diagnsotics

#         return solution_tree
#     elseif status == JuMP.MOI.INFEASIBLE
#         println("--- Model is Infeasible. Starting Conflict Analysis ---")
#         JuMP.compute_conflict!(jump_model)
#         println("The following constraints contribute to the conflict:")
#         for (name, con) in jump_constraints
#             if JuMP.get_attribute(con, JuMP.MOI.ConstraintConflictStatus()) ==
#                JuMP.MOI.IN_CONFLICT
#                 println(" - $con")
#             end
#         end
#         error("Optimization failed: Model is infeasible.")
#     else
#         error(
#             "Optimization failed with the following status and no further information is available: $status",
#         )
#     end

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

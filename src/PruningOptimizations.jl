module PruningOptimizations

using DataFrames
using DataFramesMeta
using JuMP
using COBREXA
using HiGHS
import AbstractFBCModels as A
import ConstraintTrees as C
using Printf
import MathOptInterface as MOI

export optimize_case_1,
    analyze_pruning_optimization,
    print_sinks_in_model,
    optimization_failure_analysis,
    optimize_case_4,
    optimize_case_3,
    find_sinks_on_ct

"""
    jump_constraint(m, x, v::C.Value, b::C.EqualTo; base_name::String)

Returns a constraint based off a ConstraintTrees equality bound to a JuMP model.

# Arguments
1. `m`: JuMP model onto which the
2. `x`: JuMP variable(s)
3. `v`: `ConstraintTree` value
4. `b::C.EqualTo`: Equality bound
5. `base_name::String`: The base name to use for the variable(s)

# Returns
`JuMP.ConstraintRef`

Returns the new constraint to attach to the JuMP model.
"""
function jump_constraint(m, x, v::C.Value, b::C.EqualTo; base_name::String)
    JuMP.@constraint(m, C.substitute(v, x) == b.equal_to, base_name = "$(base_name)_eq")
end

"""
    jump_constraint(m, x, v::C.Value, b::C.Between; base_name::String)

Returns a constraint based off a ConstraintTrees interval bound to a JuMP model.

# Arguments
1. `m`: JuMP model onto which the
2. `x`: JuMP variable(s)
3. `v`: `ConstraintTree` value
4. `b::C.Between`: Between bounds
5. `base_name::String`: The base name to use for the variable(s)

# Returns
`JuMP.ConstraintRef`

Returns the new constraint to attach to the JuMP model.
"""
function jump_constraint(m, x, v::C.Value, b::C.Between; base_name::String)
    isinf(b.lower) ||
        JuMP.@constraint(m, C.substitute(v, x) >= b.lower, base_name = "$(base_name)_lb")
    isinf(b.upper) ||
        JuMP.@constraint(m, C.substitute(v, x) <= b.upper, base_name = "$(base_name)_ub")
end

"""
    bound_big_m(bound; fallback = 1000.0)

Choose a big-M from a sink bound when possible, otherwise use `fallback`. Only support `C.Between` and `C.EqualTo` bounds.

# Arguments
1. `bound`: The bound of the sink.
2. `fallback`: In case the bound is not `C.Between` or `C.EqualTo`, this is what is returned.

# Returns
`Float64`

The bound to use as big-M in the coupling constraint.
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
        print_objective_value::Bool = false,
    )

JuMP MILP for Bordbar (2016) Case 1:

``\\min \\sum_{i=1}^{m} 1_{\\Delta x_i \\neq 0}``

a sum of binary indicators, with each indicator `i` determines whether sink reaction `i` is allowed to carry flux.

# Arguments
1. `ct::C.ConstraintTree`: ConstraintTree with sinks and dx/dt metabolites bounds added.
2. `optimizer = HiGHS.Optimizer`: Reference to an optimizer.
3. `fallback_M::Float64 = 1000.0`: If a sink bound is not found when creating indocator/coupling constriants, this is the fallback value.
4. `force_first_sink_on::Bool = false`: A debugging option. If `true`, forcibly sets the first sink to have non-zero flux, which will force the corresponding indicator to 1. Defaults to `false`, which does not force any sinks, and which should be used for general sinnk pruning.
5. `force_first_sink_lb::Float64 = 0.1`: A non-zero lower bound to force the first sink on with if `force_first_sink_on` is `true`.
6. `silent::Bool = true`: If `true`, the optimizer output is silenced.
7. `write_lp_path::Union{Nothing,String} = "output/debug_case1.lp"`: A filename to write the JuMP model to for debugging. If `nothing`, does not write the debugging file.
8. `print_objective_value::Bool = false`: If `true` prints the objective value.

# Returns
`Tuple{Symbol,Union{ConstraintTree,Vector{String}}}`

Returns a tuple with two elements
1. A symbol, `:ok` or `:fail`
2. If the symbol is `:ok`, the second element is a `NamedTuple` with the fields listed below. If the symbol is `:fail`, the second element is a `Vector{String}` of conflicting constraint names or a message that no further information is available.

Elements of the successful named tuple:
1. `solution_tree`: base ConstraintTree with continuous variables substituted
2. `indicator_values`: `Dict{Symbol,Float64}` mapping sink id to binary value
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
    print_objective_value::Bool = false,
)
    sink_ids = [id for (id, _) in ct.fluxes if occursin("R_REVSK", string(id))]
    isempty(sink_ids) && error("No sink reactions matching `R_REVSK`` were found.")
    jump_model = JuMP.Model(optimizer)
    silent && JuMP.set_silent(jump_model)
    x = Vector{JuMP.VariableRef}(undef, C.variable_count(ct))
    for i in eachindex(x)
        x[i] = JuMP.@variable(jump_model, base_name = "x_$i")
    end
    C.itraverse(ct) do path, con
        ct_path = join(path, ".")
        isnothing(con.bound) ||
            jump_constraint(jump_model, x, con.value, con.bound, base_name = ct_path)
    end

    # Sink (vi) indicator (zi) coupling
    #
    #       -M_i * z_i <= v_i <= M_i * z_i
    #
    # If z_i = 0, then v_i = 0.
    # If z_i = 1, then v_i is allowed within ±M_i.
    #
    @variable(jump_model, z[sink_ids], Bin)
    for id in sink_ids
        v_expr = C.substitute(ct.fluxes[id].value, x)
        M_i = bound_big_m(ct.fluxes[id].bound; fallback = fallback_M)
        JuMP.@constraint(jump_model, v_expr <= M_i * z[id], base_name = "big_m_$(id)_ub")
        JuMP.@constraint(jump_model, v_expr >= -M_i * z[id], base_name = "big_m_$(id)_lb")
    end
    if force_first_sink_on
        forced_id = sink_ids[1]
        forced_v = C.substitute(ct.fluxes[forced_id].value, x)
        JuMP.@constraint(jump_model, forced_v >= force_first_sink_lb)
    end
    JuMP.@objective(jump_model, Min, sum(z[id] for id in sink_ids))
    if !isnothing(write_lp_path)
        mkpath(dirname(write_lp_path))
        write_to_file(jump_model, write_lp_path)
    end
    JuMP.optimize!(jump_model)
    status = JuMP.termination_status(jump_model)
    if status in [JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL]
        solved_values = JuMP.value.(x)
        solution_tree = C.substitute_values(ct, solved_values)
        indicator_values = Dict(id => JuMP.value(z[id]) for id in sink_ids)
        if print_objective_value
            @printf("objective = %.12f\n", JuMP.objective_value(jump_model))
            for id in sink_ids
                v = solution_tree.fluxes[id]
                zi = indicator_values[id]
                if !(isapprox(v, 0.0) && isapprox(zi, 0.0))
                    @printf("%s   flux = %.12f   indicator = %.12f\n", string(id), v, zi)
                end
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
        result = (
            solution_tree = solution_tree,
            indicator_values = indicator_values,
            sink_ids = sink_ids,
            jump_model = jump_model,
        )
        return :ok, result
    elseif status == JuMP.MOI.INFEASIBLE
        conflicted_constraints = optimization_failure_analysis(jump_model)
        return :fail, conflicted_constraints
    else
        return :fail, ["No further information is available."]
    end
end

"""
    optimization_failure_analysis(jump_model::JuMP.Model)

Gathers names of conflicted constraints in the provided JuMP model.

# Arguments
1. `jump_model::JuMP.Model`: Broken JuMP model

# Returns
`Vector{String}`

Returns names of detected conflicted constraints.
"""
function optimization_failure_analysis(jump_model::JuMP.Model)
    constraints = JuMP.ConstraintRef[]
    for (F, S) in JuMP.list_of_constraint_types(jump_model)
        for con in JuMP.all_constraints(jump_model, F, S)
            push!(constraints, con)
        end
    end
    JuMP.compute_conflict!(jump_model)
    model_conflict_status = JuMP.get_attribute(jump_model, MOI.ConflictStatus())
    conflicted_constraints = [
        JuMP.name(con) for con in constraints if
        MOI.get(jump_model, MOI.ConstraintConflictStatus(), con) == MOI.IN_CONFLICT
    ]
    return conflicted_constraints
end

"""
    analyze_pruning_optimization(optimization_result; atol::Float64 = 1.0e-6)

Classify sinks from an optimization result. This is the second element of the tuple returned from one of the optimization functions.

# Arguments
1. `optimize_case_1_result`: Result from one of the optimization functions.
2. `atol::Float64 = 1e-9`: Tolerance for approximate zero comparisons.

# Returns
Named tuple with fields:

1. `prune::Vector{Symbol}`: Sinks to be pruned because the carry no flux.
2. `keep::Vector{Symbol}`: Sinks to keep because they carry flux.
"""
function analyze_pruning_optimization(optimization_result; atol::Float64 = 1.0e-6)
    solution_tree = optimization_result.solution_tree
    sink_ids = optimization_result.sink_ids
    prune = []
    keep = []
    for sink_id in sink_ids
        flux = solution_tree.fluxes[sink_id]
        if isapprox(flux, 0.0; atol = atol)
            push!(prune, sink_id)
        else
            push!(keep, sink_id)
        end
    end

    # n_keep = length(keep)
    # if n_keep == 0
    #     @warn "Case 1 optimization found no sinks to keep."
    # end

    return (prune = prune, keep = keep)
end

"""
    find_sinks_on_ct(ct::C.ConstraintTree)

Find the sink reactions on the given constraint tree.

# Arguments
1. `ct::C.ConstraintTree`: Constraint tree to search.

# Returns
`Vector{Symbol}`

A vector of sink reaction ids.
"""
function find_sinks_on_ct(ct::C.ConstraintTree)
    sink_ids = [id for (id, _) in ct.fluxes if occursin("R_REVSK_", string(id))]
    return sink_ids
end

"""
    optimize_case_3(original_ct::C.ConstraintTree, unmeasured_metabolite_ids::Vector{Symbol}; optimizer = HiGHS.Optimizer, silent::Bool = true, write_lp_path::Union{Nothing,String} = "output/debug_case3.lp", print_objective_value::Bool = false)

Optimize a constraint tree for Case 3 pruning as described in the Bordbar (2016) paper. That means using JuMP to optimize for the following objective:

``\\min \\sum_{i=1}^{m} \\lvert \\Delta x_i \\rvert + \\sum_{j=1}^{n} \\lvert v_j \\rvert``

Where ``|\\Delta x_i|`` denotes magnitude of the rate of change of the unmeasured metabolites and ``|v_j|`` is the magnitude of the reaction fluxes in the network.

# Arguments
1. `original_ct::C.ConstraintTree`: The constraint tree to be optimized.
2. `unmeasured_metabolite_ids::Vector{Symbol}`: Metabolite ids that are unmeasured.
3. `optimizer = HiGHS.Optimizer`: Reference to an optimizer which defaults to HiGHS
4. `silent::Bool = true`: If left at the default of `true`, suppresses optimizer logging. Note: this logging is generally unnecessary since comprehensive diagnostics are performed in the event of an error.
5. `write_lp_path::Union{Nothing,String} = "output/debug_case3.lp"`
6. `print_objective_value::Bool = false`: If `true`, prints a logging message of the objective value after optimization.

# Returns
`Tuple{Symbol,Union{Vector{String},NamedTuple}}`

Returns a tuple with two fields:
1. `:ok` or `:fail`: The status of the optimization.
2. If first element is `:ok`, this is another named tuple as specified below. If `:fail`, this is a vector of strings with either (1) the names of broken constraints or (2) an error message with as much detail as possible.

Upon success, the named tuple returned as the second element has the following fields:
1. `solution_tree`: A constraint tree with the optimized values in it. These optimized values can be used for sink pruning.
2. `sink_ids`: The ids of the reactions that were optimized as sinks.
3. `jump_model`: The JuMP model that was optimized or failed to optimize.
"""
function optimize_case_3(
    original_ct::C.ConstraintTree,
    unmeasured_metabolite_ids::Vector{Symbol};
    optimizer = HiGHS.Optimizer,
    silent::Bool = true,
    write_lp_path::Union{Nothing,String} = "output/debug_case3.lp",
    print_objective_value::Bool = false,
)
    ct = deepcopy(original_ct)
    isempty(unmeasured_metabolite_ids) && error("No unmeasured metabolites were found")
    sink_flux_ids = find_sinks_on_ct(ct)
    isempty(sink_flux_ids) && error("No sink fluxes were found")
    jump_model = JuMP.Model(optimizer)
    silent && JuMP.set_silent(jump_model)
    JuMP.@variable(jump_model, x[1:C.variable_count(ct)])
    JuMP.@variable(jump_model, t_dx[1:length(unmeasured_metabolite_ids)] >= 0)
    JuMP.@variable(jump_model, t_sink[1:length(sink_flux_ids)] >= 0)
    delta_x_exprs = [
        C.substitute(ct.flux_stoichiometry[met_id].value, x) for
        met_id in unmeasured_metabolite_ids
    ]
    sink_flux_exprs =
        [C.substitute(ct.fluxes[sink_id].value, x) for sink_id in sink_flux_ids]
    for i in eachindex(delta_x_exprs)
        met_id = unmeasured_metabolite_ids[i]
        base = "abs_dx_" * string(met_id)
        JuMP.@constraint(
            jump_model,
            t_dx[i] >= delta_x_exprs[i],
            base_name = base * "_pos",
        )
        JuMP.@constraint(
            jump_model,
            t_dx[i] >= -delta_x_exprs[i],
            base_name = base * "_neg",
        )
    end
    for j in eachindex(sink_flux_exprs)
        sink_id = sink_flux_ids[j]
        base = "abs_sink_" * string(sink_id)
        JuMP.@constraint(
            jump_model,
            t_sink[j] >= sink_flux_exprs[j],
            base_name = base * "_pos",
        )
        JuMP.@constraint(
            jump_model,
            t_sink[j] >= -sink_flux_exprs[j],
            base_name = base * "_neg",
        )
    end
    C.itraverse(ct) do path, con
        ct_path = join(string.(path), ".")
        b = con.bound
        isnothing(b) && return
        val = C.substitute(con.value, x)
        if b isa C.EqualTo
            JuMP.@constraint(jump_model, val == b.equal_to, base_name = ct_path)
        elseif b isa C.Between
            if !isinf(b.lower)
                JuMP.@constraint(jump_model, val >= b.lower, base_name = ct_path * "_lb")
            end
            if !isinf(b.upper)
                JuMP.@constraint(jump_model, val <= b.upper, base_name = ct_path * "_ub")
            end
        else
            throw(ArgumentError("Unsupported bound type: $(typeof(b))"))
        end
    end
    JuMP.@objective(jump_model, Min, sum(t_dx) + sum(t_sink),)
    if print_objective_value
        println("Case 3 objective: min sum(abs(Δx_i)) + sum(abs(v_j))")
        println("n_unmeasured_metabolites = $(length(unmeasured_metabolite_ids))")
        println("n_sink_fluxes = $(length(sink_flux_ids))")
    end
    if !isnothing(write_lp_path)
        try
            JuMP.write_to_file(jump_model, write_lp_path)
        catch err
            @warn "Could not write model to file $write_lp_path" exception = err
        end
    end
    JuMP.optimize!(jump_model)
    status = JuMP.termination_status(jump_model)
    if status in [JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL]
        solution_tree = C.substitute_values(ct, JuMP.value.(jump_model[:x]))
        result = (
            solution_tree = solution_tree,
            sink_ids = sink_flux_ids,
            jump_model = jump_model,
        )
        return :ok, result
    elseif status == JuMP.MOI.INFEASIBLE
        conflicted_constraints = optimization_failure_analysis(jump_model)
        return :fail, conflicted_constraints
    else
        return :fail, ["status=$status No further information is available."]
    end
end

"""
    print_sinks_in_model(fba_model::A.AbstractFBCModel)

This is a diagnostic helper function to print the reaction ids of sinks in the given `A.AbstractFBCModel` to ensure sinks were indeed added. Also prints the bounds of the reactions.

# Arguments
1. `fba_model::A.AbstractFBCModel`: Model to analyze.
"""
function print_sinks_in_model(fba_model::A.AbstractFBCModel)
    for (rxn_id, rxn) in fba_model.reactions
        if occursin("R_REVSK", rxn_id)
            println(rxn_id, ": ", rxn.lower_bound, ", ", rxn.upper_bound)
        end
    end
end

end

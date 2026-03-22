using CSV
using DataFrames
using DataFramesMeta
using COBREXA
using HiGHS

include("src/UfbaSampler.jl")
using .UfbaSampler
include("src/FbaModelBuilder.jl")
using .FbaModelBuilder
include("src/PruningOptimizations.jl")
using .PruningOptimizations
include("src/MetaboliteBounds.jl")
using .MetaboliteBounds

@info "Loading base RBC-GEM"
base_rbc_gem = load_base_rbc_gem()

@info "Loading metabolite bounds"
metabolite_bounds_df = load_metabolite_bounds()

@info "Loading flux bounds overrides"
flux_bounds_overrides_df = load_flux_bounds_overrides()

@info "Loading metabolite measurement opt outs"
metabolites_to_ignore = load_metabolite_measurement_opt_outs()

@info "Create FBA model and map metabolites onto that model"
fba_model, _ = create_fba_model(
    base_rbc_gem;
    exchanges = default_exchanges(),
    flux_bounds_overrides_df = flux_bounds_overrides_df,
)

@info "Reference additive and time point"
additive = "01-Ctrl AS3"
final_time = 2
println("additive: $additive, final_time: $final_time")

# @info "Adding sinks to model"
# metabolite_status_df =
#     find_metabolite_matches(fba_model, metabolite_bounds_df, additive, final_time)
# first_sink_specifications = (
#     metabolite_status_df = metabolite_status_df,
#     additive = additive,
#     prune_zero_sinks = nothing,
#     sink_opt_outs = nothing,
# )
# first_added_sink_ids =
#     add_sinks_for_unmatched_metabolites!(fba_model, first_sink_specifications)

# print_sinks_in_model(fba_model)

# @info "Add metabolite bounds to ConstraintTree"
# case1_ct = flux_balance_constraints(fba_model)
# metabolites_to_ignore = ["g6p_c", "glc__D_c", "pyr_e", "lac__L_e"]
# metabolites_to_ignore = nothing
# add_metabolite_bounds_to_constraint_tree!(
#     case1_ct,
#     metabolite_bounds_df,
#     additive,
#     final_time;
#     metabolites_to_ignore = metabolites_to_ignore,
# )
# print_metabolite_bounds_on_constraint_tree(case1_ct)

# @info "Optimize for Case 1"
# optimize_case_1_result = optimize_case_1(
#     case1_ct;
#     force_first_sink_on = true,
#     force_first_sink_lb = 0.1,
#     print_objective_value = true,
# )
# optimize_case_1_ok_fail, optimize_case_1_result = optimize_case_1(case1_ct)
# if optimize_case_1_ok_fail == :fail
#     display(optimize_case_1_result)
#     error("Case 1 optimization failed. Conflicting constraints are listed above. Stopping.")
# end
# case_1_analysis = analyze_case_1_pruning_optimization(optimize_case_1_result)

# @info "Prune zero sinks according to Case 1"
# second_fba_model, _ = create_fba_model(
#     base_rbc_gem;
#     exchanges = default_exchanges(),
#     flux_bounds_overrides_df = flux_bounds_overrides_df,
# )
# prune_zero_sinks = case_1_analysis.prune
# second_sink_specifications = (
#     metabolite_status_df = metabolite_status_df,
#     additive = additive,
#     prune_zero_sinks = prune_zero_sinks,
#     sink_opt_outs = nothing,
# )
# second_added_sink_ids =
#     add_sinks_for_unmatched_metabolites!(second_fba_model, second_sink_specifications)
# println("Added the following sinks")
# display(second_added_sink_ids)
# second_ct = flux_balance_constraints(second_fba_model)
# add_metabolite_bounds_to_constraint_tree!(
#     second_ct,
#     metabolite_bounds_df,
#     additive,
#     final_time;
#     metabolites_to_ignore = metabolites_to_ignore,
# )
# print_metabolite_bounds_on_constraint_tree(second_ct)

function find_zero_fluxes(solution_tree)
    zero_fluxes = []
    for (reaction_id, flux) in first_solution_tree.fluxes
        if isapprox(flux, 0.0)
            push!(zero_fluxes, reaction_id)
        end
    end
    n_zero_fluxes = length(zero_fluxes)
    return zero_fluxes, n_zero_fluxes
end

@info "First test case: no sinks, no metabolite bounds"
first_ct = flux_balance_constraints(fba_model)
first_status, first_solution_tree =
    optimize_constraint_tree(first_ct, first_ct.objective.value)
println("First result: $first_status")
if !isnothing(first_solution_tree)
    first_zero_fluxes, first_n_zero_fluxes = find_zero_fluxes(first_solution_tree)
    println("n_zero_fluxes: $first_n_zero_fluxes")
else
    println("Failed so no values to display")
end

# @info "Test basic FBA with no sinks and no metabolite bounds"
# result = optimize_constraint_tree(second_ct, second_ct.objective.value)
# if !isnothing(result)
#     display(result)
# end

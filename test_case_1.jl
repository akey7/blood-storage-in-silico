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

@info "Create FBA model and map metabolites onto that model"
fba_model, _ = create_fba_model(
    base_rbc_gem;
    exchanges = default_exchanges(),
    flux_bounds_overrides_df = flux_bounds_overrides_df,
)

@info "Setting up the following model"
additive = "01-Ctrl AS3"
final_time = 3
println("additive: $additive, final_time: $final_time")

metabolite_status_df =
    find_metabolite_matches(fba_model, metabolite_bounds_df, additive, final_time)
first_sink_specifications = (
    metabolite_status_df = metabolite_status_df,
    additive = additive,
    prune_zero_sinks = nothing,
    sink_opt_outs = nothing,
)
first_added_sink_ids =
    add_sinks_for_unmatched_metabolites!(fba_model, first_sink_specifications)

# Ensure sinks were added by printing them
print_sinks_in_model(fba_model)

@info "Add metabolite bounds to ConstraintTree"
case1_ct = flux_balance_constraints(fba_model)
# case1_metabolites_to_ignore = ["g6p_c", "glc__D_c", "pyr_e", "lac__L_e"]
add_metabolite_bounds_to_constraint_tree!(
    case1_ct,
    metabolite_bounds_df,
    additive,
    final_time,
)
# print_metabolite_bounds_on_constraint_tree(case1_ct)

@info "Optimize constraint tree"
optimize_case_1_result = optimize_case_1(
    case1_ct;
    force_first_sink_on = true,
    force_first_sink_lb = 0.1,
    print_objective_value = true,
)
# optimize_case_1_result = optimize_case_1(case1_ct)
case_1_analysis = analyze_case_1_pruning_optimization(optimize_case_1_result)

@info "Prune zero sinks"
second_fba_model, _ = create_fba_model(
    base_rbc_gem;
    exchanges = default_exchanges(),
    flux_bounds_overrides_df = flux_bounds_overrides_df,
)
prune_zero_sinks = case_1_analysis.prune
second_sink_specifications = (
    metabolite_status_df = metabolite_status_df,
    additive = additive,
    prune_zero_sinks = prune_zero_sinks,
    sink_opt_outs = nothing,
)
second_added_sink_ids =
    add_sinks_for_unmatched_metabolites!(second_fba_model, second_sink_specifications)
println("Added the following sinks")
display(second_added_sink_ids)
second_ct = flux_balance_constraints(second_fba_model)
add_metabolite_bounds_to_constraint_tree!(
    second_ct,
    metabolite_bounds_df,
    additive,
    final_time,
)
print_metabolite_bounds_on_constraint_tree(second_ct)

# @info "Case 1: FBA of pruned model"
# second_ct_solution_tree = optimized_values(second_ct; optimizer = HiGHS.Optimizer)
# if isnothing(second_ct_solution_tree)
#     println("Simple optimization failed")
# else
#     println("Simple optimization succeeded!")
#     display(second_ct_solution_tree.fluxes)
# end

# @info "Custom FBA of pruned and bounded ConstraintTree"
# result = optimize_constriant_tree(second_ct, second_ct.objective.value)
# if !isnothing(result)
#     display(result)
# end

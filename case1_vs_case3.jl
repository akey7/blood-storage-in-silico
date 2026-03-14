using CSV
using DataFrames
using DataFramesMeta
using HiGHS

include("src/UfbaSampler.jl")
using .UfbaSampler
include("src/FbaModelBuilder.jl")
using .FbaModelBuilder
include("src/PruningOptimizations.jl")
using .PruningOptimizations
include("src/MetaboliteBounds.jl")
using .MetaboliteBounds

@info "Loading metabolite bounds"
metabolite_bounds_df = load_metabolite_bounds()

@info "Loading flux bounds overrides"
flux_bounds_overrides_df = load_flux_bounds_overrides()

@info "Create FBA model and map metabolites onto that model"
fba_model, fba_model_metabolites_df = create_fba_model(
    load_base_rbc_gem();
    exchanges = default_exchanges(),
    flux_bounds_overrides_df = flux_bounds_overrides_df,
)
additive = "01-Ctrl AS3"
final_time = 2
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
# for (rxn_id, rxn) in fba_model.reactions
#     if occursin("R_UNKNOWN_SK", rxn_id)
#         println(rxn_id, ": ", rxn.lower_bound, ", ", rxn.upper_bound)
#     end
# end

@info "Case 1: construct constraint tree"
case1_ct = case_1_constraint_tree(fba_model)
# case1_metabolites_to_ignore = ["g6p_c", "glc__D_c", "pyr_e", "lac__L_e"]
# add_metabolite_bounds_to_constraint_tree!(
#     case1_ct,
#     metabolite_bounds_df,
#     additive,
#     final_time,
# )
@info "Case 1: Optimize constraint tree"

# print_metabolite_bounds_on_constraint_tree(case1_ct)  # Disabled, only for debugging
# case1_optimization_tree = optimize_case_1_v2(case1_ct, case1_ct.objective.value)
case1_optimization_tree =
    milp_optimized_vars(case1_ct, case1_ct.objective.value, HiGHS.Optimizer)
display(case1_optimization_tree)

# inspect_results(case1_optimization_tree)
# nonzero_indicator_ids = check_case_1_optimization_results(case1_optimization_tree)
# display(first(nonzero_indicator_ids, 10))
# case1_zero_sinks, case1_nonzero_sinks, case1_sink_status_df =
#     analyze_pruning_optimization(case1_optimization_tree)
# case1_sink_status_filename = joinpath("output", "case1_vs_case3", "case1_sinks.csv")
# CSV.write(case1_sink_status_filename, case1_sink_status_df)
# println("Wrote $case1_sink_status_filename")

# @info "Case 3 optimization"
# case3_additive = "01-Ctrl AS3"
# case3_ct = case_3_constraint_tree(fba_model, metabolite_status_df, case3_additive)
# add_metabolite_bounds_to_constraint_tree!(
#     case3_ct,
#     metabolite_bounds_df,
#     additive,
#     final_time,
# )
# case3_pruning_optimization_result = optimize_case_3(case3_ct, case3_ct.objective.value)
# case3_zero_sinks, case3_nonzero_sinks, case3_sink_status_df =
#     analyze_pruning_optimization(case3_pruning_optimization_result)
# case3_sink_status_filename = joinpath("output", "case1_vs_case3", "case3_sinks.csv")
# CSV.write(case3_sink_status_filename, case3_sink_status_df)
# println("Wrote $case3_sink_status_filename")

# @info "Comparing Case 1 vs Case 3 sinks"
# left_df = @chain case1_sink_status_df begin
#     @rename(:case1_non_zero = :is_non_zero)
#     @rtransform(:sink_name = replace(string(:sink), "R_UNKNOWN_" => ""))
#     @select(:sink_name, :case1_non_zero)
# end
# right_df = @chain case3_sink_status_df begin
#     @rename(:case3_non_zero = :is_non_zero)
#     @rtransform(:sink_name = replace(string(:sink), "R_UNKNOWN_" => ""))
#     @select(:sink_name, :case3_non_zero)
# end
# comparison_df = @chain left_df begin
#     outerjoin(right_df; on = :sink_name)
#     @rtransform(:case1_case3_different = :case1_non_zero != :case3_non_zero)
#     @orderby(:sink_name)
# end
# comparison_filename = joinpath("output", "case1_vs_case3", "case1_vs_case3_comparison.csv")
# CSV.write(comparison_filename, comparison_df)
# println("Wrote $comparison_filename")

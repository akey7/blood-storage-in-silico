using CSV
using DataFrames
using DataFramesMeta

include("src/UfbaSampler.jl")
using .UfbaSampler
include("src/FbaModelBuilder.jl")
using .FbaModelBuilder
include("src/PruningOptimizations.jl")
using .PruningOptimizations

@info "Metabolite bounds"
metabolites_bounds_df = load_metabolite_bounds()
display(first(metabolites_bounds_df, 10))

@info "Flux bounds overrides"
flux_bounds_overrides_df = load_flux_bounds_overrides()

@info "Create FBA model and map metabolites onto that model"
fba_model, fba_model_metabolites_df = create_fba_model(
    load_base_rbc_gem();
    exchanges = default_exchanges(),
    flux_bounds_overrides_df = flux_bounds_overrides_df,
)
mapping_additive = "01-Ctrl AS3"
metabolite_status_df =
    find_metabolite_matches(fba_model, metabolites_bounds_df, mapping_additive, 2)
display(first(metabolite_status_df, 10))
first_sink_specifications = (
    metabolite_status_df = metabolite_status_df,
    additive = mapping_additive,
    prune_zero_sinks = nothing,
    sink_opt_outs = nothing,
)
first_added_sink_ids =
    add_sinks_for_unmatched_metabolites!(fba_model, first_sink_specifications)
display(first(first_added_sink_ids, 10))

@info "Case 1 optimization"
case1_ct = case_1_constraint_tree!(fba_model)
case1_pruning_optimization_result = optimize_case_1(case1_ct, case1_ct.objective.value)
case1_zero_sinks, case1_nonzero_sinks, case1_sink_status_df =
    analyze_pruning_optimization(case1_pruning_optimization_result, :case1)
case1_sink_status_filename = joinpath("output", "case1_vs_case3", "case1_sinks.csv")
CSV.write(case1_sink_status_filename, case1_sink_status_df)
println("Wrote $case1_sink_status_filename")

@info "Case 3 optimization"
case3_additive = "01-Ctrl AS3"
case3_ct = case_3_constraint_tree!(fba_model, metabolite_status_df, case3_additive)
case3_pruning_optimization_result = optimize_case_3(case3_ct, case3_ct.objective.value)
case3_zero_sinks, case3_nonzero_sinks, case3_sink_status_df =
    analyze_pruning_optimization(case3_pruning_optimization_result, :case3)
case3_sink_status_filename = joinpath("output", "case1_vs_case3", "case3_sinks.csv")
CSV.write(case3_sink_status_filename, case3_sink_status_df)
println("Wrote $case3_sink_status_filename")

@info "Comparing Case 1 vs Case 3 sinks"
left_df = @chain case1_sink_status_df begin
    @rename(:case1_non_zero = :is_non_zero)
    @rtransform(:sink_name = replace(string(:sink), "R_UNKNOWN_" => ""))
    @select(:sink_name, :case1_non_zero)
end
right_df = @chain case3_sink_status_df begin
    @rename(:case3_non_zero = :is_non_zero)
    @rtransform(:sink_name = replace(string(:sink), "R_UNKNOWN_" => ""))
    @select(:sink_name, :case3_non_zero)
end
comparison_df = @chain left_df begin
    outerjoin(right_df; on = :sink_name)
    @rtransform(:case1_case3_different = :case1_non_zero != :case3_non_zero)
    @orderby(:sink_name)
end
comparison_filename = joinpath("output", "case1_vs_case3", "case1_vs_case3_comparison.csv")
CSV.write(comparison_filename, comparison_df)
println("Wrote $comparison_filename")

using CSV
using DataFrames
using OrderedCollections
using YAML
using CairoMakie

include("src/UfbaSamplerAnalysisAndViz.jl")
using .UfbaSamplerAnalysisAndViz

num_threads = Threads.nthreads()
println("Num threads $num_threads")

@info "Loading reaction ids to strings..."
rxn_ids_to_strings_filename = joinpath("output", "rxn_ids_to_strings.yml")
rxn_ids_to_strings =
    YAML.load_file(rxn_ids_to_strings_filename; dicttype = OrderedDict{String,Any})

@info "Reading sampling file"
sampling_filename = joinpath("output", "ufba_sampling.csv")
sampling_df = CSV.read(sampling_filename, DataFrame)

# @info "Global mixed model analysis"
# global_mixed_model_test(sampling_df)

# @info "Per reaction additive, time tests"
# per_reaction_df = per_reaction_additive_time_test(sampling_df)
# display(first(per_reaction_df, 20))

@info "Heatmaps!"
effects_adj_df = reaction_additive_across_time_df(sampling_df)
display(first(effects_adj_df, 20))
effects_adj_filename = joinpath("output", "uFBA_heatmaps", "effects_adj.csv")
CSV.write(effects_adj_filename, effects_adj_df)
println("Wrote $effects_adj_filename")

# @info "Diagnosing uFBA run"
# diagnostic_df = diagnose_flux_stats(sampling_df)
# diagnostic_filename = joinpath("output", "ufba_diagnostics.csv")
# CSV.write(diagnostic_filename, diagnostic_df)
# println("Wrote $diagnostic_filename")

# @info "Writing median flux DataFrame"
# median_flux_filename = joinpath("output", "ufba_median_fluxes.csv")
# median_flux_df = calc_median_flux_df(sampling_df)
# CSV.write(median_flux_filename, median_flux_df)
# println("Wrote $median_flux_filename")

# @info "Writing flux vector DataMatrix"
# data_matrix_filename = joinpath("output", "flux_vector_data_matrix.csv")
# data_matrix_df = prepare_median_flux_vector_matrix(sampling_df)
# CSV.write(data_matrix_filename, data_matrix_df)
# println("Wrote $data_matrix_filename")

# @info "Reporting measured and unmeasured metabolites, with and without sinks"
# absolute_quant_long_filename = joinpath("output", "absolute_quant_long.csv")
# absolute_quant_long_df = CSV.read(absolute_quant_long_filename, DataFrame)
# fba_model_metabolites_filename = joinpath("output", "fba_model_metabolites.csv")
# fba_model_metabolites_df = CSV.read(fba_model_metabolites_filename, DataFrame)
# ufba_optimized_sinks_filename = joinpath("output", "ufba_sinks_optimized.csv")
# ufba_optimized_sinks_df = CSV.read(ufba_optimized_sinks_filename, DataFrame)

# measurements_and_sinks_report_df, measurements_and_sinks_report_by_model_df =
#     prepare_measurements_and_sinks_report_df(
#         absolute_quant_long_df,
#         fba_model_metabolites_df,
#         ufba_optimized_sinks_df,
#         sampling_df,
#     )
# measurements_and_sinks_report_filename =
#     joinpath("output", "measurements_and_sinks_report.csv")
# CSV.write(measurements_and_sinks_report_filename, measurements_and_sinks_report_df)
# println("Wrote $measurements_and_sinks_report_filename")
# measurements_and_sinks_report_by_model_filename =
#     joinpath("output", "measurements_and_sinks_report_by_model.csv")
# CSV.write(
#     measurements_and_sinks_report_by_model_filename,
#     measurements_and_sinks_report_by_model_df,
# )
# println("Wrote $measurements_and_sinks_report_by_model_filename")

# @info "Comparing control vs. treatment fluxes"
# comparison_result =
#     compare_flux_distributions(sampling_df; alpha = 0.01, interesting_cohen_effect_z = 2.0)
# interesting_vs_uninteresting_df = comparison_result.interesting_vs_uninteresting_df
# control_vs_treatments_df = comparison_result.interesting_df
# ranked_df = comparison_result.ranked_df
# display(interesting_vs_uninteresting_df)
# control_vs_treatments_filename = joinpath("output", "control_vs_treatment.csv")
# CSV.write(control_vs_treatments_filename, control_vs_treatments_df)
# println("Wrote $control_vs_treatments_filename")
# ranked_filename = joinpath("output", "control_vs_treatment_ranked.csv")
# CSV.write(ranked_filename, ranked_df)
# println("Wrote $ranked_filename")

# @info "Plotting uFBA histograms"
# plot_all_histograms_for_reactions(sampling_df, rxn_ids_to_strings; bins = 80)

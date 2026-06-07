using CSV
using XLSX
using DataFrames
using DataFramesMeta

include("src/UfbaSamplerAnalysis.jl")
using .UfbaSamplerAnalysis

include("src/RInterface.jl")
using .RInterface

num_threads = Threads.nthreads()
println("Num threads $num_threads")

@info "Reading sampling file and valid additive / time combinations"
sampling_results = load_sampling_results()
sampling_df = sampling_results.sampling_df
working_models_df = sampling_results.working_models_df

@info "Loading reaction strings, subsystems, categories"
reaction_ids_to_strings_filename =
    joinpath("output", "rxn_strings_subsystems_categories.csv")
reaction_ids_to_strings_df = CSV.read(reaction_ids_to_strings_filename, DataFrame)

# Disabling because very long for second dataset
# @info "Global mixed model analysis"
# global_mixed_model_test(sampling_df)

# Disabling because very long for second dataset
# @info "Per reaction additive, time tests"
# per_reaction_df = per_reaction_additive_time_test(sampling_df)
# display(first(per_reaction_df, 20))

# @info "Diagnosing uFBA run"
# diagnostic_df = diagnose_flux_stats(sampling_df)
# diagnostic_filename = joinpath("output", "ufba_diagnostics.csv")
# CSV.write(diagnostic_filename, diagnostic_df)
# println("Wrote $diagnostic_filename")

@info "Writing median flux DataFrame"
median_flux_filename = joinpath("output", "ufba_median_fluxes.csv")
median_flux_df = calc_median_flux_df(sampling_df)
CSV.write(median_flux_filename, median_flux_df)
println("Wrote $median_flux_filename")

# @info "Writing flux vector data matrices"
# data_matrix_path = joinpath("output", "flux_vector_data_matrices")
# write_all_flux_vector_matrices(sampling_df, data_matrix_path)

@info "Reporting measured and unmeasured metabolites, with and without sinks"
absolute_quant_long_filename = joinpath("output", "absolute_quant_long.csv")
absolute_quant_long_df = CSV.read(absolute_quant_long_filename, DataFrame)
fba_model_metabolites_filename = joinpath("output", "fba_model_metabolites.csv")
fba_model_metabolites_df = CSV.read(fba_model_metabolites_filename, DataFrame)
ufba_optimized_sinks_filename = joinpath("output", "ufba_sinks_optimized.csv")
ufba_optimized_sinks_df = CSV.read(ufba_optimized_sinks_filename, DataFrame)
measurements_and_sinks_report_df, measurements_and_sinks_report_by_model_df =
    prepare_measurements_and_sinks_report_df(
        absolute_quant_long_df,
        fba_model_metabolites_df,
        ufba_optimized_sinks_df,
        sampling_df,
    )
measurements_and_sinks_report_filename =
    joinpath("output", "measurements_and_sinks_report.csv")
CSV.write(measurements_and_sinks_report_filename, measurements_and_sinks_report_df)
println("Wrote $measurements_and_sinks_report_filename")
measurements_and_sinks_report_by_model_filename =
    joinpath("output", "measurements_and_sinks_report_by_model.csv")
CSV.write(
    measurements_and_sinks_report_by_model_filename,
    measurements_and_sinks_report_by_model_df,
)
println("Wrote $measurements_and_sinks_report_by_model_filename")

@info "Comparing control vs. treatment fluxes"
comparison_result_0 =
    compare_flux_distributions(sampling_df, working_models_df; alpha = 0.01)
# comparison_result = remove_reaction_string_prefix(comparison_result_0)
# comparison_result_filename = joinpath("output", "reaction_treatment_comparison.xlsx")
# XLSX.writetable(
#     comparison_result_filename,
#     "score_ranking" => comparison_result.score_ranking_df,
#     "effects_wide" => comparison_result.effects_wide_df,
#     "significance_wide" => comparison_result.significance_wide_df,
#     "heatmap_rank" => comparison_result.heatmap_rank_df;
#     overwrite = true,
# )
# println("Wrote $comparison_result_filename")

@info "Combining reactions, metabolites, and measurements report"
interesting_metabolites_filename = joinpath("input", "interesting_metabolites.csv")
interesting_metabolites_df = CSV.read(interesting_metabolites_filename, DataFrame)
fba_reactions_metabolites_filename =
    joinpath("output", "fba_model_reactions_metabolites.csv")
metabolite_ids_names_filename = joinpath("input", "Metabolite Id to Name Map.csv")
fba_reactions_metabolites_df = CSV.read(fba_reactions_metabolites_filename, DataFrame)
metabolite_ids_names_df = CSV.read(metabolite_ids_names_filename, DataFrame)
reactions_metabolites_result = reactions_metabolites_report_dfs(
    fba_reactions_metabolites_df,
    metabolite_ids_names_df,
    reaction_ids_to_strings_df,
    measurements_and_sinks_report_df,
    interesting_metabolites_df,
)
reactions_metabolites_filename =
    joinpath("output", "reactions_metabolites_measurements.xlsx")
reactions_metabolites_df = reactions_metabolites_result.reactions_metabolites_df
reactions_measured_df = reactions_metabolites_result.reactions_measured_df
subsystems_measured_df = reactions_metabolites_result.subsystems_measured_df
categories_measured_df = reactions_metabolites_result.categories_measured_df
XLSX.writetable(
    reactions_metabolites_filename,
    "reactions_metabolites" => reactions_metabolites_df,
    "reactions_measured" => reactions_measured_df,
    "subsystems_measured" => subsystems_measured_df,
    "categories_measured" => categories_measured_df,
    overwrite = true,
)
println("Wrote $reactions_metabolites_filename")

@info "Determining (interesting reaction)/treatment pairs"
interesting_reactions_df = reactions_metabolites_result.interesting_reactions_df
control_vs_treatments_df = comparison_result_0.control_vs_treatments_df  # Retain R_ prefix for reaction ids
filter_reactions_treatments_result = filter_reactions_treatments_df(
    interesting_reactions_df,
    control_vs_treatments_df,
    reaction_ids_to_strings_df;
    reaction_classifications = nothing,
)
reactions_treatments_filename = joinpath("output", "analysis_reactions_treatments.xlsx")
XLSX.writetable(
    reactions_treatments_filename,
    "reactions_treatments" =>
        filter_reactions_treatments_result.interesting_reactions_treatments_df,
    "complete_conditions" => filter_reactions_treatments_result.complete_conditions_df,
    "interesting_reactions" => interesting_reactions_df,
    overwrite = true,
)
println("Wrote $reactions_treatments_filename")

@info "Reaction correlations"
corr_1_dict = reaction_correlations_one_additive_one_time(
    sampling_df,
    [:inner_reaction],
    working_models_df,
)
export_correlation_dict_for_r(corr_1_dict, "output")
println("Wrote matrices for R")

using CSV
using XLSX
using YAML
using DataFrames
using CairoMakie
using OrderedCollections

include("src/UfbaSamplerAnalysis2.jl")
using .UfbaSamplerAnalysis2

num_threads = Threads.nthreads()
println("Num threads $num_threads")

###########################################################
# FILTERING SAMPLING DATAFRAME FOR ANALYSIS PIPELINE      #
###########################################################

@info "Reading sampling file and valid additive / time combinations"
sampling_results = load_and_select_sampling_results()
sampling_df = sampling_results.sampling_df
working_models_df = sampling_results.working_models_df
working_model_filename = joinpath("output", "ufba_sampling_complete_additives.csv")
CSV.write(working_model_filename, working_models_df)
println("Wrote $working_model_filename")

@info "Loading reaction ids to strings..."
rxn_ids_to_strings_filename = joinpath("output", "rxn_ids_to_strings.yml")
rxn_ids_to_strings =
    YAML.load_file(rxn_ids_to_strings_filename; dicttype = OrderedDict{String,Any})

# Uncomment to plot from first dataset
@info "Plotting uFBA histogram and density plots"
plot_all_distributions_for_reactions(sampling_df, rxn_ids_to_strings; bins = 80)

# Select comment for the an interesting treatement and compare with control
# Especially ensure the control matches the first or second dataset.
# control_additive = "AS3"
# control_additive = "01-Ctrl AS3"
# treatment_additive = "adenine"
# treatment_additive = "arginine"
# @info "Plotting uFBA histogram and density plots for $control_additive vs $treatment_additive"
# plot_all_distributions_for_reactions(
#     sampling_df,
#     rxn_ids_to_strings;
#     control_additive = control_additive,
#     treatment_additive = treatment_additive,
#     bins = 80,
# )

###########################################################
# ANALYSIS FOR DATA VIZ PIPELINE                          #
###########################################################

@info "K-Means/PCA analysis and plots of median fluxes and Cohen's effects"
control_vs_treatments_signif_filename =
    joinpath("output", "analysis_control_vs_treatments_signif.csv")
control_vs_treatments_signif_df = CSV.read(control_vs_treatments_signif_filename, DataFrame)
prepared_treatments_result =
    prepare_treatment_effects_dfs(control_vs_treatments_signif_df; complete_only = true)
treatement_pca_result = pca_treatment_effects(prepared_treatments_result; n_pcs = 5)
treatment_k_means_df = k_means_treatment_effects(prepared_treatments_result)
median_fluxes_filename = joinpath("output", "ufba_median_fluxes.csv")
median_fluxes_df = CSV.read(median_fluxes_filename, DataFrame)
prepare_median_fluxes_result =
    prepare_median_fluxes_dfs(median_fluxes_df; complete_only = true)
fluxes_k_means_df = k_means_median_fluxes(prepare_median_fluxes_result)
fluxes_pca_result = pca_median_fluxes(prepare_median_fluxes_result)

# For first dataset
treatement_distances_df = treatment_distances_from_control(
    prepare_median_fluxes_result;
    control_additive = "01-Ctrl AS3",
)

# For second dataset
# treatement_distances_df =
#     treatment_distances_from_control(prepare_median_fluxes_result; control_additive = "AS3")

k_means_pca_filename = joinpath("output", "viz_k_means_pca_distance.xlsx")
XLSX.writetable(
    k_means_pca_filename,
    "flux_pca" => fluxes_pca_result.pca_df,
    "flux_pca_loadings" => fluxes_pca_result.loadings_df,
    "flux_k_means" => fluxes_k_means_df,
    "effects_pca" => treatement_pca_result.pca_df,
    "effects_pca_loadings" => treatement_pca_result.loadings_df,
    "effects_k_means" => treatment_k_means_df,
    "treatement_distances" => treatement_distances_df,
    overwrite = true,
)
println("Wrote $k_means_pca_filename")

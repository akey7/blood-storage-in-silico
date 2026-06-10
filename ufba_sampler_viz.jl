using CSV
using XLSX
using YAML
using DataFrames
using CairoMakie
using PlotlyJS
using OrderedCollections

include("src/UfbaSamplerViz.jl")
using .UfbaSamplerViz

num_threads = Threads.nthreads()
println("Num threads $num_threads")

# @info "Reading sampling file and valid additive / time combinations"
# sampling_results = load_sampling_results()
# sampling_df = sampling_results.sampling_df
# working_models_df = sampling_results.working_models_df

# @info "Loading reaction ids to strings..."
# rxn_ids_to_strings_filename = joinpath("output", "rxn_ids_to_strings.yml")
# rxn_ids_to_strings =
#     YAML.load_file(rxn_ids_to_strings_filename; dicttype = OrderedDict{String,Any})

# Uncomment to plot from first dataset
# @info "Plotting uFBA histogram and density plots"
# plot_all_distributions_for_reactions(sampling_df, rxn_ids_to_strings; bins = 80)

# Uncomment to select an interesting treatement and compare with control
# control_additive = "AS3"
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

# @info "3D histogram/KDE plot things"
# long_sampling_df = pivot_sampling_df_long(sampling_df)
# flux_df = @rsubset(long_sampling_df, :additive == "01-Ctrl AS3", :reaction_id == "R_ORNDC")
# display(first(flux_df, 10))
# p_kde = stacked_flux_kde_3d(flux_df)
# p_kde_filename = joinpath("output", "uFBA_3d_histograms", "line_kde_3d.html")
# savefig(p_kde, p_kde_filename)
# println("Wrote $p_kde_filename")

@info "Running k-means and plotting PCAs of Cohen's effects"
control_vs_treatments_signif_filename =
    joinpath("output", "analysis_control_vs_treatments_signif.csv")
control_vs_treatments_signif_df = CSV.read(control_vs_treatments_signif_filename, DataFrame)
prepared_treatments_result = prepare_treatment_effects_dfs(control_vs_treatments_signif_df)
treatement_pca_result = pca_treatment_effects(prepared_treatments_result)
treatement_pca_filename = joinpath("output", "viz_treatment_cluster_pca.xlsx")
treatment_k_means_df = k_means_treatment_effects(prepared_treatments_result)
XLSX.writetable(
    treatement_pca_filename,
    "pca" => treatement_pca_result.pca_df,
    "loadings" => treatement_pca_result.loadings_df,
    "k_means" => treatment_k_means_df,
    overwrite = true,
)
println("Wrote $treatement_pca_filename")

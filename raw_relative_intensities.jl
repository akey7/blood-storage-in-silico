using CairoMakie
using CSV

include("src/RawRelativeIntensities.jl")
using .RawRelativeIntensities

@info "Reading relative intensities"
long_df = load_relative_intensities()
additives = sort(unique(long_df.Additive))
# additives = ["01-Ctrl AS3"]
# for additive in additives
#     pca_result = pca_relative_intensities(relative_intensities_df, additive)
#     fig = plot_pca_panels(pca_result, "Raw Intensity PCA $additive")
#     fig_filename = joinpath("output", "pca_relative_intensity_plots", "$additive.png")
#     save(fig_filename, fig)
#     @info "Wrote $fig_filename"
# end

@info "Displaying 3D"
for additive in additives
    pca_result_3d = pca_relative_intensities(long_df, additive)
    df_filename = joinpath("output", "pca_plot_dfs", "Wide df for $additive.csv")
    CSV.write(df_filename, pca_result_3d.wide_df)
    println("Wrote $df_filename")
    display_pca_scores_3d(long_df, pca_result_3d, additive)
    println("3D for $additive. Press enter to continue")
    readline()
end

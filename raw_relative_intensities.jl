using CairoMakie

include("src/RawRelativeIntensities.jl")
using .RawRelativeIntensities

@info "Reading relative intensities"
relative_intensities_df = load_relative_intensities()
additives = sort(unique(relative_intensities_df.Additive))
# for additive in additives
#     pca_result = pca_relative_intensities(relative_intensities_df, additive)
#     fig = plot_pca_panels(pca_result, "Raw Intensity PCA $additive")
#     fig_filename = joinpath("output", "pca_relative_intensity_plots", "$additive.png")
#     save(fig_filename, fig)
#     @info "Wrote $fig_filename"
# end

@info "Displaying 3D"
for additive in additives
    pca_result_3d = pca_relative_intensities(relative_intensities_df, additive)
    display_pca_scores_3d(pca_result_3d, additive)
    println("3D for $additive. Press enter to continue")
    readline()
end

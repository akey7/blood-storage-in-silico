using CairoMakie

include("src/RawRelativeIntensities.jl")
using .RawRelativeIntensities

relative_intensities_df = load_relative_intensities()
pca_result = pca_relative_intensities(relative_intensities_df, "01-Ctrl AS3")
fig = plot_pca_panels(pca_result, "Raw PCA 01-Ctrl AS3")
fig_filename = joinpath("output", "pca_relative_intensity_plots", "control.png")
save(fig_filename, fig)
@info "Wrote $fig_filename"

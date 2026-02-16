using CairoMakie
using CSV

include("src/RawRelativeIntensities.jl")
using .RawRelativeIntensities

@info "Reading relative intensities"
long_df = load_relative_intensities()
additives = sort(unique(long_df.Additive))

@info "Aggregating loadings"
loadings_report(long_df)

# @info "Displaying 3D"
# limits = calc_pca_scores_3d_limits(long_df)
# for additive in additives
#     pca_result_3d = pca_relative_intensities(long_df, additive)
#     df_filename = joinpath("output", "pca_plot_dfs", "Wide df for $additive.csv")
#     CSV.write(df_filename, pca_result_3d.wide_df)
#     println("Wrote $df_filename")
#     display_pca_scores_3d(limits, pca_result_3d, additive)
#     println("3D for $additive. Press enter to continue")
#     readline()
# end

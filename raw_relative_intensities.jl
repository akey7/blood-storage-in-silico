using CairoMakie
using CSV

include("src/RawRelativeIntensities.jl")
using .RawRelativeIntensities

num_threads = Threads.nthreads()
println("Num threads $num_threads")

# @info "Reading ORIGINAL relative intensities"
# long_df = load_relative_intensities()
# additives = sort(unique(long_df.Additive))

@info "Reading SECOND SET OF relative intensities"
long_df = load_relative_intensities_2()
additives = sort(unique(long_df.Additive))

@info "Aggregating loadings"
loadings_df = pca_loadings_report(long_df)
loadings_filename = joinpath("output", "relative_pca_loadings.csv")
CSV.write(loadings_filename, loadings_df)
println("Wrote $loadings_filename")

@info "Plotting and Saving Single-Additive 2D PCAs"
plot_single_additive_2d_pcas(long_df)

@info "Plotting and Saving Additive Pair 2D PCAs"
plot_additive_pair_2d_pcas(long_df)

# @info "Displaying 3D PCAs"
# limits = calc_pca_scores_limits(long_df)
# for additive in additives
#     pca_result_3d = pca_relative_intensities(long_df, additive)
#     df_filename = joinpath("output", "pca_plot_dfs", "Wide df for $additive.csv")
#     CSV.write(df_filename, pca_result_3d.wide_df)
#     println("Wrote $df_filename")
#     display_pca_scores_3d(limits, pca_result_3d, additive)
#     println("3D for $additive. Press enter to continue")
#     readline()
# end

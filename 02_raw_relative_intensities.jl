using CairoMakie
using CSV

include("src/RawRelativeIntensities.jl")
using .RawRelativeIntensities

num_threads = Threads.nthreads()
println("Num threads $num_threads")

@info "Reading ORIGINAL relative intensities"
long_df = load_relative_intensities()
additives = sort(unique(long_df.Additive))

# @info "Reading SECOND SET OF relative intensities"
# long_df = load_relative_intensities_2()
# additives = sort(unique(long_df.Additive))

@info "Aggregating loadings"
loadings_df = pca_loadings_report(long_df)
loadings_filename = joinpath("output", "relative_pca_loadings.csv")
CSV.write(loadings_filename, loadings_df)
println("Wrote $loadings_filename")

@info "Plotting and Saving Single-Additive 2D PCAs"
plot_single_additive_2d_pcas(long_df)

# This will take too long for the 77 additive dataset, specify pairs explicitly?
# @info "Plotting and Saving Additive Pair 2D PCAs"
# plot_additive_pair_2d_pcas(long_df)

using CSV
using DataFrames

include("src/UfbaSamplerViz.jl")
using .UfbaSamplerViz

sampling_filename = joinpath("output", "ufba_sampling.csv")
sampling_df = CSV.read(sampling_filename, DataFrame)
plot_all_histograms(sampling_df)

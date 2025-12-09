using CSV
using DataFrames
using OrderedCollections
using YAML

include("src/UfbaSamplerViz.jl")
using .UfbaSamplerViz

@info "Reaction ids to strings..."
rxn_ids_to_strings_filename = joinpath("output", "rxn_ids_to_strings.yml")
rxn_ids_to_strings = YAML.load_file(rxn_ids_to_strings_filename; dicttype = OrderedDict{String,Any})
for (k, v) in rxn_ids_to_strings
    println(k, ": ", v)
end
# @info "Plotting histograms..."
# sampling_filename = joinpath("output", "ufba_sampling.csv")
# sampling_df = CSV.read(sampling_filename, DataFrame)
# plot_all_histograms(sampling_df)

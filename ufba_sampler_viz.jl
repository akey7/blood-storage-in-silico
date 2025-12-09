using CSV
using DataFrames
using OrderedCollections
using YAML

include("src/UfbaSamplerViz.jl")
using .UfbaSamplerViz

@info "Loading reaction ids to strings..."
rxn_ids_to_strings_filename = joinpath("output", "rxn_ids_to_strings.yml")
rxn_ids_to_strings =
    YAML.load_file(rxn_ids_to_strings_filename; dicttype = OrderedDict{String,Any})
# for (k, v) in rxn_ids_to_strings
#     println(k, ": ", v)
# end

@info "Reading sampling file"
sampling_filename = joinpath("output", "ufba_sampling.csv")
sampling_df = CSV.read(sampling_filename, DataFrame)

# @info "Plotting histograms, version 1..."
# plot_all_histograms_v1(sampling_df, rxn_ids_to_strings)

plot_all_histograms_for_reactions(sampling_df, rxn_ids_to_strings)

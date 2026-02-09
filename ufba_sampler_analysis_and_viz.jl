using CSV
using DataFrames
using OrderedCollections
using YAML

include("src/UfbaSamplerAnalysisAndViz.jl")
using .UfbaSamplerAnalysisAndViz

@info "Loading reaction ids to strings..."
rxn_ids_to_strings_filename = joinpath("output", "rxn_ids_to_strings.yml")
rxn_ids_to_strings =
    YAML.load_file(rxn_ids_to_strings_filename; dicttype = OrderedDict{String,Any})

@info "Reading sampling file"
sampling_filename = joinpath("output", "ufba_sampling.csv")
sampling_df = CSV.read(sampling_filename, DataFrame)

diagnostic_df = diagnose_flux_stats(sampling_df)
diagnostic_filename = joinpath("output", "ufba_diagnostics.csv")
CSV.write(diagnostic_filename, diagnostic_df)
println("Wrote $diagnostic_filename")

plot_all_histograms_for_reactions(sampling_df, rxn_ids_to_strings)

interesting_df = interesting_reactions_and_times(sampling_df)
interesting_filename = joinpath("output", "interesting_ufba_reactions_times.csv")
CSV.write(interesting_filename, interesting_df)
println("Wrote $interesting_filename")

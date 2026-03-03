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

@info "Diagnosing uFBA run"
diagnostic_df = diagnose_flux_stats(sampling_df)
diagnostic_filename = joinpath("output", "ufba_diagnostics.csv")
CSV.write(diagnostic_filename, diagnostic_df)
println("Wrote $diagnostic_filename")

@info "Mapping metabolites to sinks"
net_sink_fluxes_df = net_sink_fluxes(sampling_df)
net_sink_flux_filename = joinpath("output", "net_sink_fluxes.csv")
CSV.write(net_sink_flux_filename, net_sink_fluxes_df)
println("Wrote $net_sink_flux_filename")

# @info "Writing median flux DataFrame"
# median_flux_filename = joinpath("output", "ufba_median_fluxes.csv")
# median_flux_df = calc_median_flux_df(sampling_df)
# CSV.write(median_flux_filename, median_flux_df)
# println("Wrote $median_flux_filename")

# @info "Writing flux vector DataMatrix"
# data_matrix_filename = joinpath("output", "flux_vector_data_matrix.csv")
# data_matrix_df = prepare_median_flux_vector_matrix(sampling_df)
# CSV.write(data_matrix_filename, data_matrix_df)
# println("Wrote $data_matrix_filename")

# @info "Reporting measured and unmeasured metabolites, with and without sinks"
# absolute_quant_long_filename = joinpath("output", "absolute_quant_long.csv")
# absolute_quant_long_df = CSV.read(absolute_quant_long_filename, DataFrame)
# measurements_and_sinks_report_df =
#     prepare_measurements_and_sinks_report_df(sink_map, absolute_quant_long_df)
# display(first(measurements_and_sinks_report_df, 10))

# @info "Plotting uFBA histograms"
# plot_all_histograms_for_reactions(sampling_df, rxn_ids_to_strings)

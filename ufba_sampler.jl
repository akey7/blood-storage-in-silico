using CSV
using ArgParse

include("src/UfbaSampler.jl")
using .UfbaSampler

model = create_3p_model()
metabolites_bounds_df = load_metabolite_bounds()

# println(query_metabolite_bounds(metabolites_bounds_df, "01-Ctrl AS3", "cys__L_c", 2))

s = ArgParseSettings()
@add_arg_table! s begin
    "--nchains"
    help = "Number of chains during sampling"
    arg_type = Int64
    default = 10
end
n_chains = parse_args(s)["nchains"]

sampling_df, status_df =
    ufba_result =
        ufba_all_additives_all_times(model, metabolites_bounds_df; n_chains = n_chains)
println("\n############################################################")
println("# uFBA: FINAL STATUS                                       #")
println("############################################################")
display(status_df)
status_filename = joinpath("output", "ufba_sampling_status.csv")
CSV.write(status_filename, status_df)
println("Wrote $status_filename")
sampling_filename = joinpath("output", "ufba_sampling.csv")
CSV.write(sampling_filename, sampling_df)
println("Wrote $sampling_filename")
plot_all_histograms(sampling_df)

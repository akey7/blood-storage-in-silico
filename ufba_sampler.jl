using CSV
using ArgParse

include("src/UfbaSampler.jl")
using .UfbaSampler
metabolites_bounds_df = load_metabolite_bounds()

# println(query_metabolite_bounds(metabolites_bounds_df, "01-Ctrl AS3", "cys__L_c", 2))

s = ArgParseSettings()
@add_arg_table! s begin
    "--nchains"
    help = "Number of chains during sampling"
    arg_type = Int64
    default = 10
end
@add_arg_table s begin
    "--nmodels"
    help = "Number of uFBA models to analyze (-1 for all possible models)"
    arg_type = Int64
    default = -1
end
n_chains = parse_args(s)["nchains"]
n_models = parse_args(s)["nmodels"]

# _, standard_sampling_df = fba(model; n_chains = n_chains)
# standard_sampling_filename = joinpath("output", "standard_sampling.csv")
# if !isnothing(standard_sampling_df)
#     CSV.write(standard_sampling_filename, standard_sampling_df)
#     println("Wrote $standard_sampling_filename")
# else
#     println("Sampling failed, could not write")
#     exit(1)
# end

ufba_jobs = make_ufba_models_for_additives_and_times(metabolites_bounds_df, n_models)
all_sampling_df, status_df = execute_all_ufba_jobs(ufba_jobs, n_chains)
display(status_df)

# sampling_df, status_df, all_metabolite_status_df =
#     ufba_result = ufba_all_additives_all_times(
#         three_p_ufba_model,
#         metabolites_bounds_df;
#         n_chains = n_chains,
#     )

# @info "uFBA: Metabolite status"
# # display(first(all_metabolite_status_df, 10))
# all_metabolite_status_filename = joinpath("output", "all_metabolite_status.csv")
# CSV.write(all_metabolite_status_filename, all_metabolite_status_df)
# @info "uFBA: Final status"
# display(status_df)
# status_filename = joinpath("output", "ufba_sampling_status.csv")
# CSV.write(status_filename, status_df)
# println("Wrote $status_filename")
# sampling_filename = joinpath("output", "ufba_sampling.csv")
# CSV.write(sampling_filename, sampling_df)
# println("Wrote $sampling_filename")
# plot_all_histograms(sampling_df)

using CSV
using ArgParse

include("src/UfbaSampler.jl")
using .UfbaSampler

three_p_model = create_3p_model()
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

# _, standard_sampling_df = fba(model; n_chains = n_chains)
# standard_sampling_filename = joinpath("output", "standard_sampling.csv")
# if !isnothing(standard_sampling_df)
#     CSV.write(standard_sampling_filename, standard_sampling_df)
#     println("Wrote $standard_sampling_filename")
# else
#     println("Sampling failed, could not write")
#     exit(1)
# end

metabolite_status_df =
    find_metabolite_matches(three_p_model, metabolites_bounds_df, "01-Ctrl AS3", 2)
display(first(metabolite_status_df, 10))
add_sinks_for_unmatched_metabolites!(three_p_model, metabolite_status_df, "01-Ctrl AS3")
ct = case_3_constraint_tree!(three_p_model, metabolite_status_df, "01-Ctrl AS3")
optimize_case_3(ct, ct.objective.value)
# display_jump_results(branch_symbols, ct_symbols, jump_values, true)

# sampling_df, status_df, all_metabolite_status_df =
#     ufba_result =
#         ufba_all_additives_all_times(model, metabolites_bounds_df; n_chains = n_chains)

# println("\n############################################################")
# println("# uFBA: METABOLITE STATUS                                  #")
# println("############################################################")

# display(first(all_metabolite_status_df, 10))
# all_metabolite_status_filename = joinpath("output", "all_metabolite_status.csv")
# CSV.write(all_metabolite_status_filename, all_metabolite_status_df)

# println("\n############################################################")
# println("# uFBA: FINAL STATUS                                       #")
# println("############################################################")
# display(status_df)
# status_filename = joinpath("output", "ufba_sampling_status.csv")
# CSV.write(status_filename, status_df)
# println("Wrote $status_filename")
# sampling_filename = joinpath("output", "ufba_sampling.csv")
# CSV.write(sampling_filename, sampling_df)
# println("Wrote $sampling_filename")
# plot_all_histograms(sampling_df)

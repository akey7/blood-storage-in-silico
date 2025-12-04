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

three_p_ufba_model = create_3p_model(; add_exchanges = false)
metabolite_status_df =
    find_metabolite_matches(three_p_ufba_model, metabolites_bounds_df, "01-Ctrl AS3", 2)
display(first(metabolite_status_df, 10))
add_sinks_for_unmatched_metabolites!(
    three_p_ufba_model,
    metabolite_status_df,
    "01-Ctrl AS3",
    nothing,
)

println(">>>>>>>> TRANSPORTERS <<<<<<<<")
display(sort([rxn for (rxn, _) in three_p_ufba_model.reactions if contains(rxn, "t")]))
println(">>>>>>>> SINKS <<<<<<<<")
display(sort([rxn for (rxn, _) in three_p_ufba_model.reactions if contains(rxn, "SK")]))

ct = case_3_constraint_tree!(three_p_ufba_model, metabolite_status_df, "01-Ctrl AS3")
case_3_optimize_result_ct = optimize_case_3(ct, ct.objective.value)
if isnothing(case_3_optimize_result_ct)
    @info "Failed to optimize case 3, exiting"
    exit(1)
end
zero_case3_sinks, nonzero_case3_sinks = analyze_case_3(case_3_optimize_result_ct)
println(">>>>>>>>> ZERO CASE 3 SINKS <<<<<<<<<")
display(zero_case3_sinks)
println(">>>>>>>>> NON-ZERO CASE 3 SINKS <<<<<<<<<")
display(nonzero_case3_sinks)
three_p_ufba_model = create_3p_model(; add_exchanges = false)
add_sinks_for_unmatched_metabolites!(
    three_p_ufba_model,
    metabolite_status_df,
    "01-Ctrl AS3",
    string.(zero_case3_sinks),
)

sampling_df, status_df, all_metabolite_status_df =
    ufba_result = ufba_all_additives_all_times(
        three_p_ufba_model,
        metabolites_bounds_df;
        n_chains = n_chains,
    )

println("\n############################################################")
println("# uFBA: METABOLITE STATUS                                  #")
println("############################################################")

display(first(all_metabolite_status_df, 10))
all_metabolite_status_filename = joinpath("output", "all_metabolite_status.csv")
CSV.write(all_metabolite_status_filename, all_metabolite_status_df)

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
# plot_all_histograms(sampling_df)

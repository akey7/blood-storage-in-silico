using CSV
using ArgParse
using YAML

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

init_workers!()

three_p_model = create_fba_model(load_base_rbc_gem(); add_exchanges = false)
metabolite_status_df =
    find_metabolite_matches(three_p_model, metabolites_bounds_df, "01-Ctrl AS3", 2)
add_sinks_for_unmatched_metabolites!(
    three_p_model,
    metabolite_status_df,
    "01-Ctrl AS3",
    nothing,
)
rxn_ids_to_strings = map_reaction_ids_to_reaction_strings(three_p_model)
rxn_ids_to_strings_filename = joinpath("output", "rxn_ids_to_strings.yml")
YAML.write_file(rxn_ids_to_strings_filename, rxn_ids_to_strings)
@info "Wrote $rxn_ids_to_strings_filename"

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
case3_sinks_df = extract_case3_sinks(ufba_jobs)
sampling_df, status_df = execute_all_ufba_jobs(ufba_jobs, n_chains)

@info "uFBA: Final status"
display(status_df)
status_filename = joinpath("output", "ufba_sampling_status.csv")
CSV.write(status_filename, status_df)
println("Wrote $status_filename")
@info "Writing sampling results"
sampling_filename = joinpath("output", "ufba_sampling.csv")
CSV.write(sampling_filename, sampling_df)
println("Wrote $sampling_filename")
case3_sinks_filename = joinpath("output", "case3_sinks.csv")
CSV.write(case3_sinks_filename, case3_sinks_df)
println("Wrote $case3_sinks_filename")

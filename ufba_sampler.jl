using CSV
using ArgParse
using YAML

include("src/UfbaSampler.jl")
using .UfbaSampler
include("src/FbaModelBuilder.jl")
using .FbaModelBuilder
include("src/MetaboliteBounds.jl")
using .MetaboliteBounds

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

@info "Load flux bounds overrides"
flux_bounds_overrides_df = load_flux_bounds_overrides()

@info "Create reaction ids to strings mapping and save FBA model metabolites"
fba_model, fba_model_metabolites_df = create_fba_model(
    load_base_rbc_gem();
    exchanges = default_exchanges(),
    flux_bounds_overrides_df = flux_bounds_overrides_df,
)
mapping_additive = "01-Ctrl AS3"
metabolite_status_df =
    find_metabolite_matches(fba_model, metabolites_bounds_df, mapping_additive, 2)
mapping_sink_specifications = (
    metabolite_status_df = metabolite_status_df,
    additive = mapping_additive,
    prune_zero_sinks = nothing,
    sink_opt_outs = nothing,
)
add_sinks_for_unmatched_metabolites!(fba_model, mapping_sink_specifications)
rxn_ids_to_strings_dict, rxn_ids_to_strings_df =
    map_reaction_ids_to_reaction_strings(fba_model)
rxn_ids_to_strings_filename = joinpath("output", "rxn_ids_to_strings.yml")
YAML.write_file(rxn_ids_to_strings_filename, rxn_ids_to_strings_dict)
println("Wrote $rxn_ids_to_strings_filename")
fba_model_metabolites_filename = joinpath("output", "fba_model_metabolites.csv")
CSV.write(fba_model_metabolites_filename, fba_model_metabolites_df)
println("Wrote $fba_model_metabolites_filename")

@info "Run uFBA jobs"

ufba_jobs = make_ufba_models_for_additives_and_times(
    metabolites_bounds_df,
    n_models;
    exchanges = default_exchanges(),
    flux_bounds_overrides_df = flux_bounds_overrides_df,
)

# ufba_jobs = make_ufba_models_for_additives_and_times(
#     metabolites_bounds_df,
#     n_models;
#     exchanges = as3_exchanges(),
#     flux_bounds_overrides_df = flux_bounds_overrides_df,
# )

sink_status_df = extract_sinks(ufba_jobs)
added_sink_ids_df = extract_added_sink_ids(ufba_jobs)
ufba_jobs_result =
    execute_all_ufba_jobs(ufba_jobs, rxn_ids_to_strings_df; n_chains = n_chains)
sampling_df = ufba_jobs_result.sampling_df
status_df = ufba_jobs_result.status_df
status_counts_df = ufba_jobs_result.status_counts_df
joined_blocked_reaction_ids_df = ufba_jobs_result.joined_blocked_reaction_ids_df
prune_breaks_df = ufba_jobs_result.prune_breaks_df
fba_breaks_df = ufba_jobs_result.fba_breaks_df
@info "uFBA: Final status"
display(status_df)
status_filename = joinpath("output", "ufba_sampling_status.csv")
CSV.write(status_filename, status_df)
println("Wrote $status_filename")
display(status_counts_df)
@info "Writing sampling results"
sampling_filename = joinpath("output", "ufba_sampling.csv")
CSV.write(sampling_filename, sampling_df)
println("Wrote $sampling_filename")
sink_status_filename = joinpath("output", "sink_status.csv")
CSV.write(sink_status_filename, sink_status_df)
println("Wrote $sink_status_filename")
blocked_reactions_filename = joinpath("output", "ufba_blocked_reactions.csv")
CSV.write(blocked_reactions_filename, joined_blocked_reaction_ids_df)
println("Wrote $blocked_reactions_filename")
added_sink_ids_filename = joinpath("output", "ufba_added_sink_ids.csv")
CSV.write(added_sink_ids_filename, added_sink_ids_df)
println("Wrote $added_sink_ids_filename")
prune_breaks_filename = joinpath("output", "ufba_prune_breaks.csv")
CSV.write(prune_breaks_filename, prune_breaks_df)
println("Wrote $prune_breaks_filename")
fba_breaks_filename = joinpath("output", "ufba_fba_breaks.csv")
CSV.write(fba_breaks_filename, fba_breaks_df)
println("Wrote $fba_breaks_filename")

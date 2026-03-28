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
@add_arg_table s begin
    "--prune-method"
    help = "Prune method, case1 or case3"
    arg_type = Symbol
    default = :case3
end
n_chains = parse_args(s)["nchains"]
n_models = parse_args(s)["nmodels"]
prune_method = parse_args(s)["prune-method"]

init_workers!()

@info "Load flux bounds overrides"
flux_bounds_overrides_df = load_flux_bounds_overrides()

@info "Load metabolite measurement opt-outs (if available)"
metabolites_to_ignore = load_metabolite_measurement_opt_outs()

@info "Loading sink opt-ins (if available)"
sink_opt_ins = load_sink_opt_ins()

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
    sink_opt_ins = sink_opt_ins,
    metabolites_to_ignore = metabolites_to_ignore,
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
    metabolites_to_ignore = metabolites_to_ignore,
    prune_method = prune_method,
    relax_strategy = :tenth_minimum,
    sink_opt_ins = sink_opt_ins,
)

sink_overview_df = extract_sink_overview(ufba_jobs)
unmeasured_relaxations_df = extract_unmeasured_relaxations(ufba_jobs)
ufba_jobs_result =
    execute_all_ufba_jobs(ufba_jobs, rxn_ids_to_strings_df; n_chains = n_chains)
sampling_df = ufba_jobs_result.sampling_df
status_df = ufba_jobs_result.status_df
status_counts_df = ufba_jobs_result.status_counts_df
joined_blocked_reaction_ids_df = ufba_jobs_result.joined_blocked_reaction_ids_df
prune_breaks_df = ufba_jobs_result.prune_breaks_df
fba_breaks_df = ufba_jobs_result.fba_breaks_df
sinks_df = ufba_jobs_result.sinks_df

@info "uFBA: Final status"
display(status_df)
display(status_counts_df)

# Write all the files
@info "Writing sampling results"
status_filename = joinpath("output", "ufba_sampling_status.csv")
CSV.write(status_filename, status_df)
println("Wrote $status_filename")
sampling_filename = joinpath("output", "ufba_sampling.csv")
CSV.write(sampling_filename, sampling_df)
println("Wrote $sampling_filename")
sink_overview_filename = joinpath("output", "ufba_sink_overview.csv")
CSV.write(sink_overview_filename, sink_overview_df)
println("Wrote $sink_overview_filename")
blocked_reactions_filename = joinpath("output", "ufba_blocked_reactions.csv")
CSV.write(blocked_reactions_filename, joined_blocked_reaction_ids_df)
println("Wrote $blocked_reactions_filename")
prune_breaks_filename = joinpath("output", "ufba_prune_breaks.csv")
CSV.write(prune_breaks_filename, prune_breaks_df)
println("Wrote $prune_breaks_filename")
fba_breaks_filename = joinpath("output", "ufba_fba_breaks.csv")
CSV.write(fba_breaks_filename, fba_breaks_df)
println("Wrote $fba_breaks_filename")
unmeasured_relaxations_filename = joinpath("output", "ufba_unmeasured_relaxations.csv")
CSV.write(unmeasured_relaxations_filename, unmeasured_relaxations_df)
println("Wrote $unmeasured_relaxations_filename")
if !isnothing(sinks_df)
    sinks_filename = joinpath("output", "ufba_sinks_optimized.csv")
    CSV.write(sinks_filename, sinks_df)
    println("Wrote $sinks_filename")
else
    println("No sinks reported as optimized.")
end

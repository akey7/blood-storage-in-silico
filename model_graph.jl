using CSV
using DataFrames
using JSON3
using ProgressMeter
using COBREXA
using OrderedCollections
using YAML

include("src/ModelGraph.jl")
using .ModelGraph
include("src/UfbaSamplerAnalysisAndViz.jl")
using .UfbaSamplerAnalysisAndViz

# Common metabolites not to be traversed
common_metabolite_ids = [
    "M_pi_c",
    "M_h_c",
    "M_atp_c",
    "M_amp_c",
    "M_adp_c",
    "M_nad_c",
    "M_nadh_c",
    "M_nadp_c",
    "M_nadph_c",
    "M_ppi_c",
    # "M_admarg__L_c",
    "M_h2o_c",
    "M_h_e",
    "M_k_c",
    "M_k_e",
    "M_pi_e",
    "M_na1_c",
    "M_na1_e",
    "M_nh4_e",
    "M_nh4_c",
    "M_nh3_e",
    "M_nh3_c",
    "M_o2_c",
    # "M_o2_e",
    "M_co2_c",
    "M_co2_e",
]

@info "Loading uFBA models"
ufba_models = load_ufba_models()

@info "Loading uFBA sampling data"
sampling_filename = joinpath("output", "ufba_sampling.csv")
sampling_df = CSV.read(sampling_filename, DataFrame)

@info "Loading reaction ids to strings..."
rxn_ids_to_strings_filename = joinpath("output", "rxn_ids_to_strings.yml")
rxn_ids_to_strings =
    YAML.load_file(rxn_ids_to_strings_filename; dicttype = OrderedDict{String,Any})

@info "Creating uFBA model graphs"
ufba_model_graphs = make_graphs_for_ufba_models(ufba_models)

@info "Running graph search plan on all uFBA models"
dfs_plan_filename = joinpath("input", "dfs_plan.csv")
dfs_plan = CSV.read(dfs_plan_filename, DataFrame)
println("Read DFS plan from $dfs_plan_filename")

@info "Calculating median uFBA fluxes"
median_df = calc_median_flux_df(sampling_df)

@info "Analyzing visited metabolites and reactions"
visited_metabolite_filename = joinpath("output", "gem_dfs", "visited_metabolite.csv")
visited_reaction_filename = joinpath("output", "gem_dfs", "visited_reaction.csv")
plan_results = run_all_dfs_plans(ufba_model_graphs, dfs_plan, common_metabolite_ids)
enriched_visited_reactions_df =
    enrich_visited_reactions_df(plan_results.visited_reaction_df, rxn_ids_to_strings)
CSV.write(visited_metabolite_filename, plan_results.visited_metabolite_df)
println("Wrote $visited_metabolite_filename")
CSV.write(visited_reaction_filename, enriched_visited_reactions_df)
println("Wrote $visited_reaction_filename")

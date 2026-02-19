using CSV
using DataFrames
using JSON3
using ProgressMeter
import SBMLFBCModels as S
import AbstractFBCModels as A
using COBREXA

include("src/ModelGraph.jl")
using .ModelGraph
include("src/FbaModelBuilder.jl")
using .FbaModelBuilder

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
ufba_model_folder = joinpath("output", "ufba_models")
sbml_files = filter(
    f -> endswith(lowercase(f), ".xml") && isfile(joinpath(ufba_model_folder, f)),
    readdir(ufba_model_folder),
)
sbml_paths = joinpath.(ufba_model_folder, sbml_files)
n_sbml_paths = length(sbml_paths)
prog_sbml_paths = Progress(n_sbml_paths, "Loading uFBA models")
ufba_models = Dict()
for sbml_path in sbml_paths
    bn = replace(basename(sbml_path), ".xml" => "", "uFBA " => "")
    additive, final_time_str = split(bn, "_")
    final_time = parse(Int64, final_time_str)
    ufba_model = load_model(S.SBMLFBCModel, sbml_path, A.CanonicalModel.Model)
    ufba_models[(additive, final_time)] = ufba_model
    next!(prog_sbml_paths)
end

# @info "Loading uFBA sampling data"
# sampling_filename = joinpath("output", "ufba_sampling.csv")
# sampling_df = CSV.read(sampling_filename, DataFrame)

# @info "Placing uFBA results onto graph"
# graph_data = make_graph(fba_model; skip_exchanges = false)

# @info "Running graph search"
# dfs_plan_filename = joinpath("input", "dfs_plan.csv")
# dfs_plan = CSV.read(dfs_plan_filename, DataFrame)
# println("Read DFS plan from $dfs_plan_filename")

# dfs_plan_result = run_dfs_plan(dfs_plan, graph_data, common_metabolite_ids)
# dfs_plan_result_summary = summarize_dfs_plan_result(dfs_plan_result)

# df_plan_result_filename = joinpath("output", "gem_dfs", "dfs_plan_result.json")
# open(df_plan_result_filename, "w") do io
#     JSON3.pretty(io, dfs_plan_result)
# end

# dfs_plan_result_summary_filename =
#     joinpath("output", "gem_dfs", "dfs_plan_result_summary.json")
# open(dfs_plan_result_summary_filename, "w") do io
#     JSON3.pretty(io, dfs_plan_result_summary)
# end

# println("Wrote $df_plan_result_filename and $dfs_plan_result_summary_filename")

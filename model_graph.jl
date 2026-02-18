include("src/ModelGraph.jl")
using .ModelGraph

@info "Performing graph search"

graph_data = make_graph(model; skip_exchanges = true)
println("metabolite_ids placed into a graph.")

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
    "M_admarg__L_c",
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
    "M_o2_e",
    "M_co2_c",
    "M_co2_e",
]

dfs_plan_filename = joinpath("input", "gem_dfs", "dfs_plan.csv")
dfs_plan = CSV.read(dfs_plan_filename, DataFrame)
println("Read DFS plan from $dfs_plan_filename")

dfs_plan_result = run_dfs_plan(dfs_plan, graph_data, common_metabolite_ids)
dfs_plan_result_summary = summarize_dfs_plan_result(dfs_plan_result)

df_plan_result_filename = joinpath("output", "gem_dfs", "dfs_plan_result.json")
open(df_plan_result_filename, "w") do io
    JSON3.pretty(io, dfs_plan_result)
end

dfs_plan_result_summary_filename = joinpath("output", "gem_dfs", "dfs_plan_result_summary.json")
open(dfs_plan_result_summary_filename, "w") do io
    JSON3.pretty(io, dfs_plan_result_summary)
end

println("Wrote $df_plan_result_filename and $dfs_plan_result_summary_filename")

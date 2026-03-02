using CSV
using Random
using CairoMakie

include("src/AbsoluteQuant.jl")
using .AbsoluteQuant

num_threads = Threads.nthreads()
println("Num threads $num_threads")

Random.seed!(123)

println(">" ^ 10, " COMBINING RELATIVE AND ABSOLUTE QUANT DATA ", "<" ^ 10)
absolute_quant_df, absolute_quant_medians_df = load_absolute_quant()
absolute_extracellular_quant_df = load_extracellular_absolute_quant()
absolute_quant_filename = joinpath("output", "absolute_quant.csv")
CSV.write(absolute_quant_filename, absolute_quant_df)
println(first(absolute_quant_medians_df, 10))
fold_changes_df = load_relative_quant()
qc_fold_changes_df, qc_fold_change_zeros_df = qc(fold_changes_df)
qc_fold_changes_filename = joinpath("output", "qc_fold_changes.csv")
CSV.write(qc_fold_changes_filename, qc_fold_changes_df)
qc_fold_change_zeros_filename = joinpath("output", "qc_fold_change_zeros.csv")
CSV.write(qc_fold_change_zeros_filename, qc_fold_change_zeros_df)
long_df, wide_df =
    combine_relative_and_absolute_quant_c(fold_changes_df, absolute_quant_medians_df)
long_df_filename = joinpath("output", "absolute_quant_c_long_df.csv")
CSV.write(long_df_filename, long_df)
println("Wrote $long_df_filename")
relative_absolute_quant_filename = joinpath("output", "absolute_quant_c_wide_df.csv")
CSV.write(relative_absolute_quant_filename, wide_df)
println("Wrote $relative_absolute_quant_filename")
absolute_quant_e_long_df =
    combine_relative_and_absolute_quant_e(fold_changes_df, absolute_extracellular_quant_df)
# absolute_quant_e_long_filename = joinpath("output", "absolute_quant_e_long.csv")
# absolute_quant_e_wide_filename = joinpath("output", "absolute_quant_e_wide.csv")
# CSV.write(absolute_quant_e_long_filename, absolute_quant_e_long_df)
# CSV.write(absolute_quant_e_wide_filename, absolute_quant_e_wide_df)
# println("Wrote $absolute_quant_e_long_filename")
# println("Wrote $absolute_quant_e_wide_filename")

# println(">" ^ 10, " TIMESERIES PLOTS ", "<" ^ 10)
# plot_all_mM_timeseries(long_df)

# println(">" ^ 10, " C-MEANS CLUSTERING ", "<" ^ 10)
# all_memberships_dfs, fuzzy_objectives_df =
#     cluster_all_additives_all_n_clusters(long_df; max_clusters = 7)
# println(first(fuzzy_objectives_df, 10))

# println(">" ^ 10, " MAKING C-MEANS PLOTS ", "<" ^ 10)
# plot_elbows(fuzzy_objectives_df)
# all_primary_cluster_df = plot_c_means_all_additives(long_df, all_memberships_dfs, 6)
# all_primary_cluster_df_filename = joinpath("output", "c_means_primary_clusters.csv")
# CSV.write(all_primary_cluster_df_filename, all_primary_cluster_df)
# println("Wrote $all_primary_cluster_df_filename")

# println(">" ^ 10, " RATE REGRESSIONS ", "<" ^ 10)
# rate_df = regress_concentration_vs_time(long_df)
# rate_filename = joinpath("output", "concentration_rates.csv")
# CSV.write(rate_filename, rate_df)
# fig_conc_vs_time = plot_all_regressions(long_df)

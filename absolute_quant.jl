using CSV
using Random
using CairoMakie

include("src/AbsoluteQuant.jl")
using .AbsoluteQuant

num_threads = Threads.nthreads()
println("Num threads $num_threads")

Random.seed!(123)

println(">" ^ 10, " WRANGLING DATA ", "<" ^ 10)
absolute_quant_df, absolute_quant_medians_df = load_absolute_quant()
absolute_quant_filename = joinpath("output", "absolute_quant.csv")
CSV.write(absolute_quant_filename, absolute_quant_df)
println(first(absolute_quant_medians_df, 10))
fold_changes_df = load_relative_quant()
println(first(fold_changes_df, 10))
long_df, wide_df =
    combine_relative_and_absolute_quant(fold_changes_df, absolute_quant_medians_df)
println(first(long_df, 10))
long_df_filename = joinpath("output", "long_df.csv")
CSV.write(long_df_filename, long_df)
relative_absolute_quant_filename = joinpath("output", "relative_absolute_quant.csv")
CSV.write(relative_absolute_quant_filename, wide_df)
println("Wrote $relative_absolute_quant_filename")

# println(">" ^ 10, " TIMESERIES PLOTS ", "<" ^ 10)
# plot_all_mM_timeseries(long_df)

# println(">" ^ 10, " DIFFING ", "<" ^ 10)
# diffed_df = diff_mM(long_df)
# diff_filename = joinpath("output", "diff_mM.csv")
# CSV.write(diff_filename, diffed_df)

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

println(">" ^ 10, " PCA ANALYSIS ", "<" ^ 10)
plot_pca_all_additives(long_df)

# println(">" ^ 10, " RATE REGRESSION ", "<" ^ 10)
# rate_df = regress_concentration_dxdt(long_df, 1000)
# display(first(rate_df, 20))
# rate_filename = joinpath("output", "concentration_rates.csv")
# CSV.write(rate_filename, rate_df)

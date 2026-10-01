using CSV
using Random
using CairoMakie
using DataFrames

include("src/AbsoluteQuant.jl")
using .AbsoluteQuant

num_threads = Threads.nthreads()
println("Num threads $num_threads")

Random.seed!(123)

@info "Combining relative and absolute quant"
absolute_quant_df, absolute_quant_medians_df = load_absolute_quant()
absolute_extracellular_quant_df = load_extracellular_absolute_quant()
fold_changes_df = load_relative_quant()
fold_change_filename = joinpath("output", "relative_quant_1", "fold_changes.csv")
CSV.write(fold_change_filename, fold_changes_df)
display(first(fold_changes_df, 10))
println("Wrote $fold_change_filename")
qc_fold_changes_df, qc_fold_change_zeros_df = qc(fold_changes_df)
qc_fold_changes_filename = joinpath("output", "qc_fold_changes.csv")
CSV.write(qc_fold_changes_filename, qc_fold_changes_df)
qc_fold_change_zeros_filename = joinpath("output", "qc_fold_change_zeros.csv")
CSV.write(qc_fold_change_zeros_filename, qc_fold_change_zeros_df)
absolute_quant_c_long_df =
    combine_relative_and_absolute_quant_c(fold_changes_df, absolute_quant_medians_df)
absolute_quant_e_long_df =
    combine_relative_and_absolute_quant_e(fold_changes_df, absolute_extracellular_quant_df)
absolute_quant_long_df, absolute_quant_wide_df = union_and_pivot_wide(
    absolute_quant_c_long_df,
    absolute_quant_e_long_df;
    include_extracellular = true,
)
absolute_quant_long_filename = joinpath("output", "absolute_quant_long.csv")
CSV.write(absolute_quant_long_filename, absolute_quant_long_df)
println("Write $absolute_quant_long_filename")
absolute_quant_wide_filename = joinpath("output", "absolute_quant_wide.csv")
CSV.write(absolute_quant_wide_filename, absolute_quant_wide_df)
println("Wrote $absolute_quant_wide_filename")

# @info "Relative quant, SECOND dataset"
# fold_changes_2_result = load_relative_quant_2()
# fold_changes_df = fold_changes_2_result.fold_changes_df
# metabolite_cleaning_df = fold_changes_2_result.metabolite_cleaning_df
# clean_metabolite_cleaning_df = copy(metabolite_cleaning_df)
# clean_metabolite_cleaning_df[!, :col_name] =
#     [replace(string(s), "’" => "'") for s in clean_metabolite_cleaning_df[!, :col_name]]
# n_samples_per_condition_df = fold_changes_2_result.n_samples_per_condition_df
# fold_changes_filename = joinpath("output", "relative_quant_2", "fold_changes.csv")
# metabolite_cleaning_filename =
#     joinpath("output", "relative_quant_2", "metabolite_cleaning.csv")
# n_samples_per_condition_filename =
#     joinpath("output", "relative_quant_2", "n_samples_per_condition.csv")
# CSV.write(fold_changes_filename, fold_changes_df)
# CSV.write(metabolite_cleaning_filename, clean_metabolite_cleaning_df)
# CSV.write(n_samples_per_condition_filename, n_samples_per_condition_df)
# println("Wrote $fold_changes_filename")
# println("Wrote $metabolite_cleaning_filename")
# println("Wrote $n_samples_per_condition_filename")
# absolute_quant_df, absolute_quant_medians_df = load_absolute_quant()
# absolute_extracellular_quant_df = load_extracellular_absolute_quant()
# absolute_quant_c_long_df =
#     combine_relative_and_absolute_quant_c(fold_changes_df, absolute_quant_medians_df)
# absolute_quant_e_long_df =
#     combine_relative_and_absolute_quant_e(fold_changes_df, absolute_extracellular_quant_df)
# absolute_quant_long_df, absolute_quant_wide_df = union_and_pivot_wide(
#     absolute_quant_c_long_df,
#     absolute_quant_e_long_df;
#     include_extracellular = true,
# )
# absolute_quant_long_filename = joinpath("output", "absolute_quant_long.csv")
# CSV.write(absolute_quant_long_filename, absolute_quant_long_df)
# println("Write $absolute_quant_long_filename")
# absolute_quant_wide_filename = joinpath("output", "absolute_quant_wide.csv")
# CSV.write(absolute_quant_wide_filename, absolute_quant_wide_df)
# println("Wrote $absolute_quant_wide_filename")

n_plots = 100

@info "Rate regressions"
rate_df = regress_concentration_vs_time(absolute_quant_long_df; remove_zero_rates = true)
rate_filename = joinpath("output", "concentration_rates.csv")
CSV.write(rate_filename, rate_df)
println("Wrote $rate_filename")

@info "Regression plots (limiting to first $n_plots)"
plot_all_regressions(absolute_quant_long_df; n_plots = n_plots)

# For the second dataset there are too many additives to plot in this way using given color palette
# @info "Timeseries plots"
# plot_all_mM_timeseries(absolute_quant_long_df)

# For the second relative quant dataset, ways to visualize and save c-means clusters needs to be scaled to higher numbers of additives.
# max_clusters = 7
# @info "C-Means clustering (max clusters: $max_clusters)"
# all_memberships_dfs, fuzzy_objectives_df = cluster_all_additives_all_n_clusters(
#     absolute_quant_long_df;
#     max_clusters = max_clusters,
# )
# println(first(fuzzy_objectives_df, 10))

# @info "Making c-means plots"
# plot_elbows(fuzzy_objectives_df)
# all_primary_cluster_df =
#     plot_c_means_all_additives(absolute_quant_long_df, all_memberships_dfs, 6)
# all_primary_cluster_df_filename = joinpath("output", "c_means_primary_clusters.csv")
# CSV.write(all_primary_cluster_df_filename, all_primary_cluster_df)
# println("Wrote $all_primary_cluster_df_filename")

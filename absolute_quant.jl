using CSV
using Random
using CairoMakie

include("src/AbsoluteQuant.jl")
using .AbsoluteQuant

num_threads = Threads.nthreads()
println("Num threads $num_threads")

Random.seed!(123)

# @info "Combining relative and absolute quant, FIRST dataset"
# absolute_quant_df, absolute_quant_medians_df = load_absolute_quant()
# absolute_extracellular_quant_df = load_extracellular_absolute_quant()
# fold_changes_df = load_relative_quant()
# fold_change_filename = joinpath("output", "relative_quant_1", "fold_changes.csv")
# CSV.write(fold_change_filename, fold_changes_df)
# display(first(fold_changes_df, 10))
# println("Wrote $fold_change_filename")
# qc_fold_changes_df, qc_fold_change_zeros_df = qc(fold_changes_df)
# qc_fold_changes_filename = joinpath("output", "qc_fold_changes.csv")
# CSV.write(qc_fold_changes_filename, qc_fold_changes_df)
# qc_fold_change_zeros_filename = joinpath("output", "qc_fold_change_zeros.csv")
# CSV.write(qc_fold_change_zeros_filename, qc_fold_change_zeros_df)
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

@info "Relative quant, SECOND dataset"
fold_changes_2_result = load_relative_quant_2()
valid_measurements_df = fold_changes_2_result.valid_measurements_df
as3_fold_change_df = fold_changes_2_result.as3_fold_change_df
conditions_fold_change_df = fold_changes_2_result.conditions_fold_change_df
all_conditions_long_df = fold_changes_2_result.all_conditions_long_df
fold_changes_df_2 = fold_changes_2_result.fold_changes_df
cols_to_keep = fold_changes_2_result.cols_to_keep
metabolite_cleaning_df = fold_changes_2_result.metabolite_cleaning_df
# display(first(fold_changes_df_2, 10))
valid_measurements_filename =
    joinpath("output", "relative_quant_2", "valid_measurements.csv")
as3_fold_change_filename = joinpath("output", "relative_quant_2", "as3_fold_change.csv")
conditions_fold_change_filename = joinpath("output", "relative_quant_2", "conditions_fold_change.csv")
all_conditions_long_filename = joinpath("output", "relative_quant_2", "all_conditions_long_df.csv")
fold_changes_filename_2 = joinpath("output", "relative_quant_2", "fold_changes_2.csv")
CSV.write(valid_measurements_filename, valid_measurements_df)
CSV.write(as3_fold_change_filename, as3_fold_change_df)
CSV.write(conditions_fold_change_filename, conditions_fold_change_df)
CSV.write(all_conditions_long_filename, all_conditions_long_df)
CSV.write(fold_changes_filename_2, fold_changes_df_2)
println("Wrote $valid_measurements_filename")
println("Wrote $as3_fold_change_filename")
println("Wrote $conditions_fold_change_filename")
println("Wrote $all_conditions_long_filename")
println("Wrote $fold_changes_filename_2")
absolute_quant_df, absolute_quant_medians_df = load_absolute_quant()
absolute_extracellular_quant_df = load_extracellular_absolute_quant()
absolute_quant_c_long_df =
    combine_relative_and_absolute_quant_c(fold_changes_df_2, absolute_quant_medians_df)
absolute_quant_e_long_df =
    combine_relative_and_absolute_quant_e(fold_changes_df_2, absolute_extracellular_quant_df)
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
display(first(metabolite_cleaning_df, 10))

# For the second dataset there are too many additives to plot in this way
# @info "Timeseries plots"
# plot_all_mM_timeseries(absolute_quant_long_df)

# @info "C-Means clustering"
# all_memberships_dfs, fuzzy_objectives_df =
#     cluster_all_additives_all_n_clusters(absolute_quant_long_df; max_clusters = 7)
# println(first(fuzzy_objectives_df, 10))

# @info "Making c-means plots"
# plot_elbows(fuzzy_objectives_df)
# all_primary_cluster_df =
#     plot_c_means_all_additives(absolute_quant_long_df, all_memberships_dfs, 6)
# all_primary_cluster_df_filename = joinpath("output", "c_means_primary_clusters.csv")
# CSV.write(all_primary_cluster_df_filename, all_primary_cluster_df)
# println("Wrote $all_primary_cluster_df_filename")

# @info "Rate regressions"
# rate_df = regress_concentration_vs_time(absolute_quant_long_df)
# rate_filename = joinpath("output", "concentration_rates.csv")
# CSV.write(rate_filename, rate_df)
# println("Wrote $rate_filename")
# plot_all_regressions(absolute_quant_long_df)

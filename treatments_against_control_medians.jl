using Base.Threads
using Random
using CSV

include("src/TreatmentsAgainstControlMedians.jl")
using .TreatmentsAgainstControlMedians

num_threads = Threads.nthreads()
println("Num threads $num_threads")

Random.seed!(123)

df = load_and_clean_2()

# df_filename = joinpath("output", "control_median_normalized_intensity.csv")
# CSV.write(df_filename, df)
# println("Wrote $df_filename")
# plot_loess_for_all_metabolites(df)
# results_df = find_significant_metabolites_additives(df)
# results_filename = joinpath("output", "signficant_pairs.csv")
# CSV.write(results_filename, results_df)
# println("Wrote $results_filename")

c_means_dict, wide_timeseries_dict, fuzzy_objectives_df =
    c_means_metabolite_trajectories(df)
c_means_filename = joinpath("output", "c_means_clusters.csv")
CSV.write(c_means_filename, c_means_dict[5])
println("Wrote $c_means_filename")
plot_c_means_for_all_additives(c_means_dict[5], wide_timeseries_dict[5])
println(fuzzy_objectives_df)
plot_fuzzy_objectives_elbow(fuzzy_objectives_df)

using Base.Threads
using Random
using CSV
using DataFrames

include("src/TreatmentsAgainstControlMedians.jl")
using .TreatmentsAgainstControlMedians

num_threads = Threads.nthreads()
println("Num threads $num_threads")

Random.seed!(123)

normalized_intensity_df = load_and_clean_2()
df_filename = joinpath("output", "control_median_normalized_intensity.csv")
CSV.write(df_filename, normalized_intensity_df)
println("Wrote $df_filename")

# all_c_means_df, all_wide_timeseries_df, fuzzy_objectives_df =
#     c_means_metabolite_trajectories(normalized_intensity_df, 10)
# c_means_filename = joinpath("output", "c_means_clusters.csv")
# CSV.write(c_means_filename, all_c_means_df)
# println("Wrote $c_means_filename")
# plot_c_means_for_all_additives(7, all_c_means_df, all_wide_timeseries_df)
# println(fuzzy_objectives_df)
# plot_fuzzy_objectives_elbow(fuzzy_objectives_df)

# pathways_filename = joinpath("input", "Pathways 1.csv")
# pathways_df = CSV.read(pathways_filename, DataFrame)
# cluster_enrichment_df = cluster_enrichment_analysis(7, all_c_means_df, pathways_df)
# cluster_enrichment_filename = joinpath("output", "cluster_enrichment.csv")
# CSV.write(cluster_enrichment_filename, cluster_enrichment_df)
# println("Wrote $cluster_enrichment_filename")

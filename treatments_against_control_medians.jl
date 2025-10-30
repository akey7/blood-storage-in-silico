using Base.Threads
using Random
using CSV
using DataFrames

include("src/TreatmentsAgainstControlMedians.jl")
using .TreatmentsAgainstControlMedians

num_threads = Threads.nthreads()
println("Num threads $num_threads")

Random.seed!(123)

println(">" ^ 10, " LOADING METABOLOMICS DATA ", "<" ^ 10)
normalized_intensity_df = load_and_clean_2()
df_filename = joinpath("output", "control_median_normalized_intensity.csv")
CSV.write(df_filename, normalized_intensity_df)
println("Wrote $df_filename")

println(">" ^ 10, " CALCULATING C-MEANS CLUSTERS ", "<" ^ 10)
all_c_means_df, all_wide_timeseries_df, fuzzy_objectives_df =
    c_means_metabolite_trajectories(normalized_intensity_df, 10)
c_means_filename = joinpath("output", "c_means_clusters.csv")
CSV.write(c_means_filename, all_c_means_df)
println("Wrote $c_means_filename")
plot_c_means_for_all_additives(7, all_c_means_df, all_wide_timeseries_df)
println(fuzzy_objectives_df)
plot_fuzzy_objectives_elbow(fuzzy_objectives_df)

println(">" ^ 10, " LOADING GEMS AND MAPPING SUBSYSTEMS ", "<" ^ 10)
gem_reactions_df, gem_metabolites_df = load_gem_and_subsystems()
println(first(gem_reactions_df, 10))
println(first(gem_metabolites_df, 10))

println(">" ^ 10, " REACTION CATEGORY ENRICHMENT ANALYSIS ", "<" ^ 10)
enrichment_df, metabolites_subsystems_df, top3_df, category_colors =
    cluster_enrichment_analysis(7, all_c_means_df, gem_reactions_df, gem_metabolites_df)
println(first(enrichment_df, 10))
println(first(metabolites_subsystems_df, 10))
println(first(top3_df, 10))
println(category_colors)
plot_bars_for_all_additives(top3_df, category_colors)

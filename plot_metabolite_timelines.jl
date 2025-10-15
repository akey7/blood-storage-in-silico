using Base.Threads
using CSV

include("src/MetaboliteTimelines.jl")
using .MetaboliteTimelines

num_threads = Threads.nthreads()
println("Num threads $num_threads")

df_01 = load_and_clean_01()
df_02 = load_and_clean_02()

normalized_abundance_correlations_df = normalized_abundance_correlations(df_01)
println(first(normalized_abundance_correlations_df, 100))
normalized_abundance_correlations_filename =
    joinpath("output", "normalized_abundance_correlations.csv")
CSV.write(normalized_abundance_correlations_filename, normalized_abundance_correlations_df)
plot_aggregations_for_all_metabolites(df_02)

time_trend_df = normalized_abundance_correlations_by_additive(df_01)
println(first(time_trend_df, 100))
time_trend_filename = joinpath("output", "time_trends.csv")
CSV.write(time_trend_filename, time_trend_df)

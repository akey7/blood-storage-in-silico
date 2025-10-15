using PlotlyBase
using Base.Threads
using CSV

include("src/MetaboliteTimelines.jl")
using .MetaboliteTimelines

num_threads = Threads.nthreads()
println("Num threads $num_threads")

df_01 = load_and_clean_01()
df_02 = load_and_clean_02()

normalized_abundance_correlations_df = normalized_abundance_correlations(df_01)
println(first(normalized_abundance_correlations_df, 200))
normalized_abundance_correlations_filename = joinpath("output", "normalized_abundance_correlations.csv")
CSV.write(normalized_abundance_correlations_filename, normalized_abundance_correlations_df)
plot_aggregations_for_all_metabolites(df_02)

additive_and_time_df = normalized_abundance_correlations_by_additive_and_time(df_01)
additive_and_time_filename = joinpath("output", "additive_and_time.csv")
CSV.write(additive_and_time_filename, additive_and_time_df)

using PlotlyBase
using Base.Threads

include("src/MetaboliteTimelines.jl")
using .MetaboliteTimelines

num_threads = Threads.nthreads()
println("Num threads $num_threads")

df_01 = load_and_clean_01()
normalized_abundance_correlations_df = normalized_abundance_correlations(df_01)
println(first(normalized_abundance_correlations_df, 200))

df_02 = load_and_clean_02()
plot_aggregations_for_all_metabolites(df_02)

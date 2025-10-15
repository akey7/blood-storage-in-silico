using Base.Threads
using CSV

include("src/MetaboliteTimelines.jl")
using .MetaboliteTimelines

num_threads = Threads.nthreads()
println("Num threads $num_threads")

df = load_and_clean()

normalized_abundance_correlations_df = normalized_abundance_correlations(df)
println(first(normalized_abundance_correlations_df, 100))
normalized_abundance_correlations_filename =
    joinpath("output", "normalized_abundance_correlations.csv")
CSV.write(normalized_abundance_correlations_filename, normalized_abundance_correlations_df)
plot_aggregations_for_all_metabolites(df)

using CSV

include("src/AbsoluteQuant.jl")
using .AbsoluteQuant

num_threads = Threads.nthreads()
println("Num threads $num_threads")

println(">" ^ 10, " WRANGLING DATA ", "<" ^ 10)
absolute_quant_df, absolute_quant_medians_df = load_absolute_quant()
absolute_quant_filename = joinpath("output", "absolute_quant.csv")
CSV.write(absolute_quant_filename, absolute_quant_df)
println(first(absolute_quant_medians_df, 10))
fold_changes_df = load_relative_quant()
println(first(fold_changes_df, 10))
long_df, wide_df = combine_relative_and_absolute_quant(fold_changes_df, absolute_quant_medians_df)
println(first(long_df, 10))
relative_absolute_quant_filename = joinpath("output", "relative_absolute_quant.csv")
CSV.write(relative_absolute_quant_filename, wide_df)
println("Wrote $relative_absolute_quant_filename")

println(">" ^ 10, " C-MEANS CLUSTERING ", "<" ^ 10)
all_memberships_dfs, fuzzy_objectives_df = cluster_all_additives_all_n_clusters(long_df)
println(fuzzy_objectives_df)

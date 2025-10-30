include("src/AbsoluteQuant.jl")
using .AbsoluteQuant

num_threads = Threads.nthreads()
println("Num threads $num_threads")

absolute_quant_df, normalization_df = load_absolute_quant()
println(first(normalization_df, 10))

relative_fold_change_df = load_relative_quant()
println(first(relative_fold_change_df, 10))

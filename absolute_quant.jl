include("src/AbsoluteQuant.jl")
using .AbsoluteQuant

num_threads = Threads.nthreads()
println("Num threads $num_threads")

absolute_quant_df = load_and_clean_3()
println(first(absolute_quant_df, 10))

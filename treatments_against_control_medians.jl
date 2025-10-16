using Base.Threads
using CSV

include("src/TreatmentsAgainstControlMedians.jl")
using .TreatmentsAgainstControlMedians

num_threads = Threads.nthreads()
println("Num threads $num_threads")

df = load_and_clean_2()
df_filename = joinpath("output", "control_median_normalized_intensity.csv")
CSV.write(df_filename, df)
println("Wrote $df_filename")
plot_loess_for_all_metabolites(df)

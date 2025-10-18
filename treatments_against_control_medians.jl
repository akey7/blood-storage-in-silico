using Base.Threads
using CSV

include("src/TreatmentsAgainstControlMedians.jl")
using .TreatmentsAgainstControlMedians

num_threads = Threads.nthreads()
println("Num threads $num_threads")

df = load_and_clean_2()

# df_filename = joinpath("output", "control_median_normalized_intensity.csv")
# CSV.write(df_filename, df)
# println("Wrote $df_filename")
# plot_loess_for_all_metabolites(df)
# results_df = find_significant_metabolites_additives(df)
# results_filename = joinpath("output", "signficant_pairs.csv")
# CSV.write(results_filename, results_df)

# c_means_df = c_means_metabolite_trajectories_in_additive(df, "01-Ctrl AS3")
c_means_df = c_means_metabolite_trajectories_in_additive(df, "02-Adenosine")
println(first(c_means_df, 10))

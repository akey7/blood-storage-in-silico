using CSV
using DataFrames

include("src/UfbaSampler.jl")
using .UfbaSampler

model = create_3p_model()
metabolites_bounds_df = load_metabolite_bounds()

# println(query_metabolite_bounds(metabolites_bounds_df, "01-Ctrl AS3", "cys__L_c", 2))

sampling_df, status_df =
    ufba_result = ufba_all_additive_all_times(model, metabolites_bounds_df)
println("\n############################################################")
println("# uFBA: FINAL STATUS                                       #")
println("############################################################")
display(status_df)
status_filename = joinpath("output", "ufba_sampling_status.csv")
CSV.write(status_filename, status_df)
println("Wrote $status_filename")
sampling_filename = joinpath("output", "ufba_sampling.csv")
CSV.write(sampling_filename, sampling_df)
println("Wrote $sampling_filename")
fig = histograms_for_reaction_in_additive(sampling_df, "01-Ctrl AS3", "R_HXPRT")
display(fig)
println("Press enter to exit...")
readline()

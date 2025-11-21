using CSV
using DataFrames

include("src/UfbaSampler.jl")
using .UfbaSampler

model = create_3p_model()
# fluxes_df = sample_fluxes(model)
# display(first(fluxes_df[!, :R_HEX1], 10))
metabolites_bounds_df = load_metabolite_bounds()
# println(query_metabolite_bounds(metabolites_bounds_df, "01-Ctrl AS3", "23dpg_c", 2))
# println(query_metabolite_bounds(metabolites_bounds_df, "01-Ctrl AS3", "nonexistent_c", 2))
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

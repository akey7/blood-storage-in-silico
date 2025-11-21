using DataFrames

include("src/UfbaSampler.jl")
using .UfbaSampler

model = create_3p_model()
# fluxes_df = sample_fluxes(model)
# display(first(fluxes_df[!, :R_HEX1], 10))
metabolites_bounds_df = load_metabolite_bounds()
# println(query_metabolite_bounds(metabolites_bounds_df, "01-Ctrl AS3", "23dpg_c", 2))
# println(query_metabolite_bounds(metabolites_bounds_df, "01-Ctrl AS3", "nonexistent_c", 2))
ufba_result = ufba(model, metabolites_bounds_df, "01-Ctrl AS3", 2)
if !isnothing(ufba_result)
    println("Generated $(nrow(ufba_result)) samples")
end

include("src/UfbaSampler.jl")
using .UfbaSampler

model = create_3p_model()
# constraints_explorer(model)
fluxes_df = sample_fluxes(model)
# display(first(fluxes_df[!, :R_HEX1], 10))
# convert_to_jump(model)
ufba(model)

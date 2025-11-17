include("src/UfbaSampler.jl")
using .UfbaSampler

model = create_3p_model()
fluxes_df = sample_fluxes(model)
display(first(fluxes_df, 10))

using HiGHS

include("src/UfbaSampler.jl")
using .UfbaSampler

model = create_3p_model()
# fluxes_df = sample_fluxes(model)
# display(first(fluxes_df[!, :R_HEX1], 10))
ufba_result = ufba(model; optimizer = HiGHS.Optimizer)
display(ufba_result)

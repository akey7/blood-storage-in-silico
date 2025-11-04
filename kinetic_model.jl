include("src/KineticModel.jl")
using .KineticModel

glycolysis_network, sol = run_glycolysis()

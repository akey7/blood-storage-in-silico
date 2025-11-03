include("src/KineticModel.jl")
using .KineticModel

println(glycolysis())
@info "Press enter to exit..."
readline()

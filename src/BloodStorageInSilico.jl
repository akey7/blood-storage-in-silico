module BloodStorageInSilico

include("MetaboliteTimelines.jl")
include("TreatmentsAgainstControlMedians.jl")
include("AbsoluteQuant.jl")
include("KineticModel.jl")

export hello_world,
    MetaboliteTimelines, TreatmentsAgainstControlMedians, AbsoluteQuant, KineticModel

hello_world() = println("Hello world!")

end

module BloodStorageInSilico

include("MetaboliteTimelines.jl")
include("TreatmentsAgainstControlMedians.jl")
include("AbsoluteQuant.jl")

export hello_world, MetaboliteTimelines, TreatmentsAgainstControlMedians, AbsoluteQuant

hello_world() = println("Hello world!")

end

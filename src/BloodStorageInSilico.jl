module BloodStorageInSilico

include("MetaboliteTimelines.jl")
include("TreatmentsAgainstControlMedians.jl")

export hello_world, MetaboliteTimelines, TreatmentsAgainstControlMedians

hello_world() = println("Hello world!")

end

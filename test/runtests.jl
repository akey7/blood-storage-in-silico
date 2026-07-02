using Test

@testset "BloodStorageInSilico Core Pipeline Initialization" begin
    @test isfile(joinpath("input", "RBC-GEM.xml"))
end

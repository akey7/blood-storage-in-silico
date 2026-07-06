using Test

@testset "BloodStorageInSilico Core Pipeline Initialization" begin
    repo_root = normpath(joinpath(@__DIR__, ".."))
    model_path = joinpath(repo_root, "input", "RBC-GEM.xml")
    @test isfile(model_path)
end

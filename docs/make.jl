using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.instantiate()

using Documenter
using BloodStorageInSilico

makedocs(
    sitename = "BloodStorageInSilico.jl",
    modules = [
        BloodStorageInSilico.RawRelativeIntensities,
        BloodStorageInSilico.UfbaSampler,
        BloodStorageInSilico.UfbaSamplerViz,
    ],
)

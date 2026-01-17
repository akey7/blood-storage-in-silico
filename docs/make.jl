using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.instantiate()

using Documenter
using BloodStorageInSilico

ENV["JULIA_DOCUMENTER_BUILD"] = "true"

makedocs(
    sitename = "UfbaSampler.jl",
    modules = [BloodStorageInSilico.UfbaSampler],
)

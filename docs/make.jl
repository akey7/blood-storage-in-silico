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
        BloodStorageInSilico.AbsoluteQuant,
        BloodStorageInSilico.UfbaSampler,
        BloodStorageInSilico.UfbaSamplerAnalysisAndViz,
        BloodStorageInSilico.MetaboliteTimelines,
        BloodStorageInSilico.ModelGraph,
        BloodStorageInSilico.UfbaSampler.FbaModelBuilder,
        BloodStorageInSilico.UfbaSampler.PruningOptimizations,
        BloodStorageInSilico.UfbaSampler.MetaboliteBounds,
    ],
    format = Documenter.HTML(
        prettyurls = false,
    ),
    pages = [
        "Home" => "index.md",
        "Quantification Workflow" => "quantification_workflow.md",
        "uFBA Workflow" => "ufba_workflow.md",
        "uFBA Analysis Workflow" => "ufba_analysis_workflow.md",
    ],
)

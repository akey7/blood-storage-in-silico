using Makie
using CairoMakie
using AlgebraOfGraphics

CairoMakie.activate!()

include("src/MetaboliteTimelines.jl")
using .MetaboliteTimelines

df = load_and_clean_metabolite_timelines()
plt = plot_timeline_for_metabolite(df, "L-alanine")
fig = draw(
    plt;
    axis=(; title="L-alanine"),
    legend=(;
        position=:right,
        orientation=:vertical,
        title="Additive",
        labelsize=14,
        titlesize=16,
    )
)
plt_filename = joinpath("output", "L-alanine.png")
save(plt_filename, fig)
println(first(df, 10))

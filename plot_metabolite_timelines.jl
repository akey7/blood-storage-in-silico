using CairoMakie
using AlgebraOfGraphics

CairoMakie.activate!()

include("src/MetaboliteTimelines.jl")
using .MetaboliteTimelines

df = load_and_clean_metabolite_timelines()
plt = plot_timeline_for_metabolite(df, "L-alanine")
fig = Figure(; size = (1000, 1000))
ax = Axis(fig[1, 1]; xlabel = "Time", ylabel = "Abundance", title = "L-alanine")
draw!(ax, plt)
plt_filename = joinpath("output", "L-alanine.png")
save(plt_filename, fig)
println(first(df, 10))

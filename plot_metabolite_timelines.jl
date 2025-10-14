using PlotlyBase

include("src/MetaboliteTimelines.jl")
using .MetaboliteTimelines

df = load_and_clean_all_normalized_abundances()
println(first(df, 10))
plt = plot_scatter_all_normalized_abundances(df, "L-alanine")
plt_save_path = joinpath("output", "L-alanine.html")
open(plt_save_path, "w") do io
    PlotlyBase.to_html(io, plt.plot, include_plotlyjs = "cdn")
end
println("Wrote $plt_save_path")
display(plt)
println("Press enter to exit...")
readline()

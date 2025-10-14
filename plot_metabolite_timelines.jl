include("src/MetaboliteTimelines.jl")
using .MetaboliteTimelines

println(first(load_and_clean_metabolite_timelines(), 10))

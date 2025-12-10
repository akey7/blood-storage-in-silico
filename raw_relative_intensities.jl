include("src/RawRelativeIntensities.jl")
using .RawRelativeIntensities

relative_intensities_df = load_relative_intensities()
display(first(relative_intensities_df, 100))
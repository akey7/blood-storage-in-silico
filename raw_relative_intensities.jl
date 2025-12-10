include("src/RawRelativeIntensities.jl")
using .RawRelativeIntensities

relative_intensities_df = load_relative_intensities()
pca_result = pca_relative_intensities(relative_intensities_df, "01-Ctrl AS3")
display(pca_result.model)

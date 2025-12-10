module RawRelativeIntensities

using CSV
using DataFrames
using DataFramesMeta

export load_relative_intensities

function load_relative_intensities()
    relative_filename = joinpath("input", "Data Sheet 1.CSV")
    wide_df = CSV.read(relative_filename, DataFrame)
    long_df = stack(
        wide_df,
        Not([:Sample, :Time, :Additive]),
        variable_name = :MixedName,
        value_name = :Intensity,
    )
    return long_df
end

end

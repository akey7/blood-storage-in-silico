module MetaboliteTimelines

using CSV
using DataFrames

export load_and_clean_metabolite_timelines

function load_and_clean_metabolite_timelines()
    filename = joinpath("input", "Data Sheet 1.CSV")
    df_1 = CSV.read(filename, DataFrame)
    df_2 = stack(
        df_1,
        Not([:Sample, :Time, :Additive]),
        variable_name = :Metabolite,
        value_name = :PeakAreaTop,
    )
    return df_2
end

end

module MetaboliteTimelines

using CSV
using DataFrames
using DataFramesMeta
using Statistics

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
    df_3 = @combine(
        groupby(df_2, [:Time, :Additive, :Metabolite]),
        :MedianPeakAreaTop = median(skipmissing(:PeakAreaTop))
    )
    df_4 = select(df_2, Not(:Sample))
    df_5 = innerjoin(df_4, df_3, on = [:Time, :Additive, :Metabolite])
    df_6 = transform(df_5, [:PeakAreaTop, :MedianPeakAreaTop] => ByRow((x, y) -> x / y) => :Abundance)
    df_7 = select(df_6, [:Time, :Additive, :Metabolite, :Abundance])
    return df_7
end

end

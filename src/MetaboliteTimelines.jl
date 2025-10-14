module MetaboliteTimelines

using CSV
using DataFrames
using DataFramesMeta
using Statistics
using PlotlyJS

export load_and_clean_all_normalized_abundances,
    makie_plot_timeline_for_metabolite, plot_scatter_all_normalized_abundances

function load_and_clean_all_normalized_abundances()
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
    df_6 = transform(
        df_5,
        [:PeakAreaTop, :MedianPeakAreaTop] => ByRow((x, y) -> x / y) => :Abundance,
    )
    all_abundances_df = select(df_6, [:Time, :Additive, :Metabolite, :Abundance])
    return all_abundances_df
end

function plot_scatter_all_normalized_abundances(everything_df, metabolite)
    df = subset(everything_df, :Metabolite => x -> x .== metabolite)
    traces = [
        scatter(
            x = df[df.Additive .== additive, :Time],
            y = df[df.Additive .== additive, :Abundance],
            mode = "markers",
            name = string(additive),
            marker = attr(size = 10),
        ) for additive in unique(df.Additive)
    ]
    layout = Layout(
        title = "Abundance vs. Time by Additive",
        xaxis = attr(
            title = "Time",
            tickvals = unique(df.Time),
            ticktext = string.(unique(df.Time)),
        ),
        yaxis = attr(title = "Abundance"),
        legend = attr(title = "Additive"),
    )
    plt = plot(traces, layout)
    return plt
end

end

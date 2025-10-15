module MetaboliteTimelines

using CSV
using DataFrames
using DataFramesMeta
using Statistics
using Distributions
using PlotlyJS

export makie_plot_timeline_for_metabolite,
    plot_scatter_all_normalized_abundances,
    load_and_clean_02,
    plot_means_for_metabolite,
    load_and_clean_01

function load_and_clean_01()
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
        [:PeakAreaTop, :MedianPeakAreaTop] =>
            ByRow((x, y) -> x / y) => :NormalizedAbundance,
    )
    df_7 = select(df_6, [:Time, :Additive, :Metabolite, :NormalizedAbundance])
    return df_7
end

function load_and_clean_02()
    filename = joinpath("input", "Data Sheet 1.CSV")
    df_1 = CSV.read(filename, DataFrame)
    df_2 = stack(
        df_1,
        Not([:Sample, :Time, :Additive]),
        variable_name = :Metabolite,
        value_name = :Intensity,
    )
    alpha = 0.05
    abundances_df = @combine(
        groupby(df_2, [:Time, :Additive, :Metabolite]),
        :MeanIntensity = mean(skipmissing(:Intensity)),
        :SEM_min =
            mean(skipmissing(:Intensity)) -
            quantile(TDist(length(:Intensity)-1), 1-alpha/2) * std(:Intensity) /
            sqrt(length(:Intensity)),
        :SEM_max =
            mean(skipmissing(:Intensity)) +
            quantile(TDist(length(:Intensity)-1), 1-alpha/2) * std(:Intensity) /
            sqrt(length(:Intensity)),
        :MedianIntensity = median(skipmissing(:Intensity))
    )
    return abundances_df
end

function plot_means_for_metabolite(everything_df, metabolite)
    metabolite_df = subset(everything_df, :Metabolite => x -> x .== metabolite)
    traces::Vector{GenericTrace} = []
    for additive in unique(metabolite_df.Additive)
        additive_df = subset(metabolite_df, :Additive => x -> x .== additive)
        trace = scatter(
            x = additive_df.Time,
            y = additive_df.MedianIntensity,
            mode = "lines+markers",
            name = additive,
            marker = attr(size = 10),
            # error_y = attr(
            #     type = "data",
            #     array = additive_df.SEM_max,
            #     arraymin = additive_df.SEM_min,
            #     visible = true,
            # )
        )
        push!(traces, trace)
    end
    layout = Layout(
        title = "Mean $metabolite Abundance vs. Time by Additive",
        xaxis = attr(
            title = "Time",
            tickvals = unique(everything_df.Time),
            ticktext = string.(unique(everything_df.Time)),
        ),
        yaxis = attr(title = "Abundance"),
        legend = attr(title = "Additive"),
    )
    plt = plot(traces, layout)
    return plt
end

end

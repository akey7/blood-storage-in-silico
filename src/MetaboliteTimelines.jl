module MetaboliteTimelines

using CSV
using DataFrames
using DataFramesMeta
using Statistics
using PlotlyJS

export load_and_clean_all_normalized_abundances,
    makie_plot_timeline_for_metabolite,
    plot_scatter_all_normalized_abundances,
    load_and_clean_means_only,
    plot_means_for_metabolite

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

function load_and_clean_means_only()
    filename = joinpath("input", "Data Sheet 1.CSV")
    df_1 = CSV.read(filename, DataFrame)
    df_2 = stack(
        df_1,
        Not([:Sample, :Time, :Additive]),
        variable_name = :Metabolite,
        value_name = :Intensity,
    )

    # TODO: Calculate the critical value of the t distribution based on
    # degrees of freedom. Do not hard code 1.96!

    median_abundances_df = @combine(
        groupby(df_2, [:Time, :Additive, :Metabolite]),
        :MeanIntensity = mean(skipmissing(:Intensity)),
        :SEM_min =
            mean(skipmissing(:Intensity)) -
            1.96 * std(:Intensity) / sqrt(length(:Intensity)),
        :SEM_max =
            mean(skipmissing(:Intensity)) +
            1.96 * std(:Intensity) / sqrt(length(:Intensity))
    )
    return median_abundances_df
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

function plot_means_for_metabolite(everything_df, metabolite)
    metabolite_df = subset(everything_df, :Metabolite => x -> x .== metabolite)
    traces::Vector{GenericTrace} = []
    for additive in unique(metabolite_df.Additive)
        additive_df = subset(metabolite_df, :Additive => x -> x .== additive)
        trace = scatter(
            x = additive_df.Time,
            y = additive_df.MeanIntensity,
            mode = "lines+markers",
            name = additive,
            marker = attr(size = 10),
            error_y = attr(
                type = "data",
                array = additive_df.SEM_max,
                arraymin = additive_df.SEM_min,
                visible = true,
            )
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

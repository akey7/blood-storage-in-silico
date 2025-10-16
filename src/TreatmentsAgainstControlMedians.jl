module TreatmentsAgainstControlMedians

using CSV
using DataFrames
using DataFramesMeta
using Statistics
using AlgebraOfGraphics
using CairoMakie
using Makie

export load_and_clean_2, plot_loess_for_all_metabolites

function load_and_clean_2()
    filename = joinpath("input", "Data Sheet 1.CSV")
    df1 = CSV.read(filename, DataFrame)
    df2 = stack(
        df1,
        Not([:Sample, :Time, :Additive]),
        variable_name = :Metabolite,
        value_name = :Intensity,
    )
    df3 = subset(df2, :Additive => x -> x .== "01-Ctrl AS3")
    df4 = @combine(
        groupby(df3, [:Metabolite, :Time]),
        :ControlMedianIntensity = median(skipmissing(:Intensity))
    )
    df5 = innerjoin(df2, df4, on = [:Metabolite, :Time])
    df6 = transform(
        df5,
        [:Intensity, :ControlMedianIntensity] =>
            ByRow((x, y) -> x / y) => :ControlMedianNormalizedIntensity,
    )
    df7 = select(
        df6,
        [:Sample, :Time, :Additive, :Metabolite, :ControlMedianNormalizedIntensity],
    )
    return df7
end

function plot_loess_for_metabolite(everything_df, metabolite)
    df = subset(everything_df, :Metabolite => x -> x .== metabolite)
    time_points = unique(df.Time)
    plt =
        data(df) *
        mapping(
            :Time => "Time",
            :ControlMedianNormalizedIntensity => "Control Median Normalized Intensity",
            color = :Additive => "Additive",
        ) *
        (visual(Scatter; markersize = 10) + visual(Lines) + linear())
    fig = draw(
        plt;
        figure = (; size = (750, 500)),
        axis = (; title = metabolite, xticks = time_points),
    )
    return fig
end

function plot_loess_for_all_metabolites(df)
    metabolites = unique(df.Metabolite)
    for metabolite in metabolites
        fig = plot_loess_for_metabolite(df, metabolite)
        clean_metabolite = replace(metabolite, r"[^A-Za-z0-9]" => "_")
        filename = joinpath(
            "output",
            "control_median_normalized_intensity_plots",
            "$(clean_metabolite).png",
        )
        save(filename, fig)
        println("Wrote $filename")
    end
end

end

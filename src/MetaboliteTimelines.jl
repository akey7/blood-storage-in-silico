module MetaboliteTimelines

using Base.Iterators
using CSV
using DataFrames
using DataFramesMeta
using Statistics
using Distributions
using HypothesisTests
using HypothesisTests: pvalue
using StatsBase
using MultipleTesting
using Combinatorics
using ThreadsX
using AlgebraOfGraphics
using CairoMakie
using Makie

export load_and_clean,
    plot_aggregations_for_all_metabolites, normalized_abundance_correlations

function load_and_clean()
    filename = joinpath("input", "Data Sheet 1.CSV")
    df1 = CSV.read(filename, DataFrame)
    df2 = stack(
        df1,
        Not([:Sample, :Time, :Additive]),
        variable_name = :Metabolite,
        value_name = :Intensity,
    )
    df3 = @combine(
        groupby(df2, :Metabolite),
        :MedianIntensity = median(skipmissing(:Intensity))
    )
    df4 = innerjoin(df2, df3, on = :Metabolite)
    df5 = transform(
        df4,
        [:Intensity, :MedianIntensity] =>
            ByRow((x, y) -> x / y) => :MedianNormalizedIntensity,
    )
    df6 = select(df5, [:Sample, :Time, :Additive, :Metabolite, :MedianNormalizedIntensity])
    return df6
end

function aggregate_metabolite_additive(everything_df, metabolite, additive)
    additive_df = subset(
        everything_df,
        :Additive => x -> x .== additive,
        :Metabolite => x -> x .== metabolite,
    )
    aggregated_df = @combine(
        groupby(additive_df, :Time),
        :Aggregated = mean(skipmissing(:MedianNormalizedIntensity))
    )
    return aggregated_df
end

function plot_aggregations_for_metabolite(everything_df, metabolite)
    metabolite_df = subset(everything_df, :Metabolite => x -> x .== metabolite)
    aggregated_df = @combine(
        groupby(metabolite_df, [:Additive, :Time]),
        :Aggregated = mean(skipmissing(:MedianNormalizedIntensity))
    )
    time_points = unique(everything_df.Time)
    plt =
        data(aggregated_df) *
        mapping(:Time, :Aggregated, color = :Additive) *
        (visual(Lines) + visual(Scatter; markersize = 10))
    fig = draw(
        plt;
        figure = (; size = (750, 500)),
        axis = (;
            title = metabolite,
            xlabel = "Time",
            ylabel = "Normalized Abundance",
            xticks = time_points,
        ),
    )
    return fig
end

function plot_aggregations_for_all_metabolites(df)
    metabolites = unique(df.Metabolite)
    for metabolite in metabolites
        fig = plot_aggregations_for_metabolite(df, metabolite)
        clean_metabolite = replace(metabolite, r"[^A-Za-z0-9]" => "_")
        filename = joinpath("output", "plots", "$(clean_metabolite).png")
        save(filename, fig)
        println("Wrote $filename")
    end
end

function normalized_abundance_correlations(df)
    println("Calculating MedianNormalizedIntensity correlations")
    metabolites = unique(df.Metabolite)
    unique_pairs = collect(combinations(metabolites, 2))
    rows = ThreadsX.map(unique_pairs) do unique_pair
        m1, m2 = unique_pair
        m1_df = subset(df, :Metabolite => x -> x .== m1)
        m2_df = subset(df, :Metabolite => x -> x .== m2)
        xvs = tiedrank(m1_df.MedianNormalizedIntensity)
        yvs = tiedrank(m2_df.MedianNormalizedIntensity)
        spearman = CorrelationTest(xvs, yvs)
        p_value = pvalue(spearman)
        rho = spearman.r
        return (m1 = m1, m2 = m2, rho = rho, p_value = p_value)
    end
    df = DataFrame(rows)
    fdr_threshold = 0.05
    df.adj_p_value = adjust(df.p_value, BenjaminiHochberg())
    df.significant = df.adj_p_value .< fdr_threshold
    final_df = sort(df, :adj_p_value)
    return final_df
end

end

module MetaboliteTimelines

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
using PlotlyJS
using PlotlyBase

# replace(s, r"[^A-Za-z0-9]" => "_")

export makie_plot_timeline_for_metabolite,
    plot_scatter_all_normalized_abundances,
    load_and_clean_02,
    plot_aggregations_for_metabolite,
    load_and_clean_01,
    normalized_abundance_correlations,
    save_aggregation_plot,
    plot_aggregations_for_all_metabolites

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

function plot_aggregations_for_metabolite(everything_df, metabolite)
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

function save_aggregation_plot(plt, filename)
    open(filename, "w") do io
        PlotlyBase.to_html(io, plt.plot, include_plotlyjs = "cdn")
    end
    println("Wrote $filename")
end

function plot_aggregations_for_all_metabolites(df)
    metabolites = unique(df.Metabolite)
    for metabolite in metabolites[1:10]
        plt = plot_aggregations_for_metabolite(df, metabolite)
        filename = joinpath("output", "plots", "$(metabolite).html")
        save_aggregation_plot(plt, filename)
    end
end

function normalized_abundance_correlations(df)
    metabolites = unique(df.Metabolite)
    unique_pairs = collect(combinations(metabolites, 2))
    rows = ThreadsX.map(unique_pairs) do unique_pair
        m1, m2 = unique_pair
        m1_df = subset(df, :Metabolite => x -> x .== m1)
        m2_df = subset(df, :Metabolite => x -> x .== m2)
        xvs = tiedrank(m1_df.NormalizedAbundance)
        yvs = tiedrank(m2_df.NormalizedAbundance)
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

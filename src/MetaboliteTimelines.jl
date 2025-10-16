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
using PlotlyJS
using PlotlyBase

export makie_plot_timeline_for_metabolite,
    plot_scatter_all_normalized_abundances,
    plot_aggregations_for_metabolite,
    load_and_clean,
    normalized_abundance_correlations,
    save_aggregation_plot,
    plot_aggregations_for_all_metabolites,
    normalized_abundance_correlations_by_additive

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

function plot_aggregations_for_metabolite(everything_df, metabolite)
    println("Plotting $metabolite")
    traces::Vector{GenericTrace} = []
    additives = unique(everything_df.Additive)
    for additive in additives
        additive_df = subset(
            everything_df,
            :Additive => x -> x .== additive,
            :Metabolite => x -> x .== metabolite,
        )
        aggregated_df = @combine(
            groupby(additive_df, :Time),
            :Aggregated = mean(skipmissing(:MedianNormalizedIntensity))
        )
        trace = scatter(
            x = aggregated_df.Time,
            y = aggregated_df.Aggregated,
            mode = "lines+markers",
            name = additive,
            marker = attr(size = 10),
        )
        push!(traces, trace)
    end
    layout = Layout(
        title = "<b>$metabolite</b>",
        xaxis = attr(
            title = "Time",
            tickvals = unique(everything_df.Time),
            ticktext = string.(unique(everything_df.Time)),
        ),
        yaxis = attr(title = "Median Normalized Abundance"),
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
    for metabolite in metabolites
        plt = plot_aggregations_for_metabolite(df, metabolite)
        clean_metabolite = replace(metabolite, r"[^A-Za-z0-9]" => "_")
        filename = joinpath("output", "plots", "$(clean_metabolite).html")
        save_aggregation_plot(plt, filename)
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

function mann_kendall_no_ties(xs)
    vals = string(xs)
    s = 0.0
    for k in eachindex(xs)[1:(end-1)]
        for j in eachindex(xs)[2:end]
            s += sign(xs[j] - xs[k])
        end
    end
    n = length(xs)
    var_s = n*(n-1)*(2*n+5)/18.0
    if s > 0.0
        return vals, n, s, var_s, abs((s-1) / sqrt(var_s))
    elseif s < 0.0
        return vals, n, s, var_s, abs((s+1) / sqrt(var_s))
    else
        return vals, n, s, var_s, 0.0
    end
end

end

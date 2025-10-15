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
    load_and_clean_02,
    plot_aggregations_for_metabolite,
    load_and_clean_01,
    normalized_abundance_correlations,
    save_aggregation_plot,
    plot_aggregations_for_all_metabolites,
    normalized_abundance_correlations_by_additive

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
    println("Plotting $metabolite")
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
    ThreadsX.map(metabolites) do metabolite
        plt = plot_aggregations_for_metabolite(df, metabolite)
        clean_metabolite = replace(metabolite, r"[^A-Za-z0-9]" => "_")
        filename = joinpath("output", "plots", "$(clean_metabolite).html")
        save_aggregation_plot(plt, filename)
    end
end

function normalized_abundance_correlations(df)
    println("Calculating normalized abundance correlations")
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

function mann_kendall_no_ties(xs)
    s = 0.0
    for k in eachindex(xs)[1:(end-1)]
        for j in eachindex(xs)[2:end]
            s += sign(xs[j] - xs[k])
        end
    end
    n = length(xs)
    var_s = n*(n-1)*(2*n+5)/18.0
    if s > 0.0
        return n, s, var_s, abs((s-1) / sqrt(var_s))
    elseif s < 0.0
        return n, s, var_s, abs((s+1) / sqrt(var_s))
    else
        return n, s, var_s, 0.0
    end
end

function normalized_abundance_correlations_by_additive(input_df)
    println("Calculating normalized abundance by additive time trends")
    additives = unique(input_df.Additive)
    metabolites = unique(input_df.Metabolite)
    jobs = vec(collect(product(additives, metabolites)))
    rows = ThreadsX.map(jobs) do job
        additive, metabolite = job
        df1 = subset(
            input_df,
            :Metabolite => x -> x .== metabolite,
            :Additive => x -> x .== additive,
        )
        df2 = @combine(
            groupby(df1, :Time),
            :MedianNormalizedAbundance = median(skipmissing(:NormalizedAbundance))
        )
        df3 = sort(df2, :Time)
        xs = df3.MedianNormalizedAbundance
        n, s, var_s, mann_kendall_z = mann_kendall_no_ties(xs)
        is_significant = mann_kendall_z > 1.96
        return (
            additive = additive,
            metabolite = metabolite,
            n = n,
            s = s,
            var_s = var_s,
            mann_kendall_z = mann_kendall_z,
            is_significant = is_significant,
        )
    end
    rows_df = DataFrame(rows)
    final_df = sort(rows_df, [:additive, :metabolite, :is_significant])
    return final_df
end

end

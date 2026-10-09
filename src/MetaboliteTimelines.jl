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
using ProgressMeter

export load_and_clean,
    plot_aggregations_for_all_metabolites, normalized_abundance_correlations

"""
    load_and_clean()

Loads, cleans, and preprocesses relative quant metabolomics data from `input/Data Sheet 1.CSV`. This prepares the data to be used by other functions in this module.

# Returns
`DataFrame`

Returns a DataFrame, pivoted long, with the following columns:
1. `:Sample`: The sample id
2. `:Time`: Time of the measurement in weeks.
3. `:Additive`: The additive the measurement was taken in.
4. `:Metabolite`: The name of metabolite (from the original sheet). Note that this is different from other modules that translates and proportionates the metabolite names into RBC-GEM identifiers.
5. `:MedianNormalizedIntensity`: The intensity of the metabolite normalized by the median value for that metabolite at that time point.
"""
function load_and_clean()
    filename = joinpath("input", "Data Sheet 1.CSV")
    df1 = CSV.read(filename, DataFrame; ntasks = 1)
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

"""
    plot_aggregations_for_metabolite(everything_df, metabolite)

Plot the mean `:MedianNormalizedIntensity` for the given metabolite.

# Arguments
1. `everything_df`: The DataFrame returned by [`load_and_clean`](@ref BloodStorageInSilico.MetaboliteTimelines.plot_aggregations_for_metabolite)
2. `metabolite`: Name of the metabolite to plot.

# Returns
`Figure`

Returns a Makie figure that can be displayed or saved.
"""
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

"""
    plot_aggregations_for_all_metabolites(df)

Plots aggregations for all metabolites with [`plot_aggregations_for_metabolite`](@ref BloodStorageInSilico.MetaboliteTimelines.plot_aggregations_for_metabolite). Saves each file to `output/plots`, with the metabolite name "cleaned" to make a well-behaved filename. Displays a nifty status bar while generating the plots.

# Arguments
1. `df`: The DataFrame returned by [`load_and_clean`](@ref BloodStorageInSilico.MetaboliteTimelines.plot_aggregations_for_metabolite).
"""
function plot_aggregations_for_all_metabolites(df)
    metabolites = unique(df.Metabolite)
    n_metabolites = length(metabolites)
    prog = Progress(n_metabolites, "Writing metabolite timelines")
    for metabolite in metabolites
        fig = plot_aggregations_for_metabolite(df, metabolite)
        clean_metabolite = replace(metabolite, r"[^A-Za-z0-9]" => "_")
        filename = joinpath("output", "plots", "$(clean_metabolite).png")
        save(filename, fig)
        next!(prog)
    end
end

"""
    normalized_abundance_correlations(df)

Calculates the correlations and adjusted p-values of correlations of abundances for metabolites. FDR threshold 0.05.

# Arguments
1. `df`: The DataFrame returned by [`load_and_clean`](@ref BloodStorageInSilico.MetaboliteTimelines.plot_aggregations_for_metabolite). Uses ThreadsX to compute on multiple threads.

# Returns
`DataFrame`

Returns a DataFrame with the following columns:

1. `:m1`: The first metabolite name
2. `:m2`: The second metabolite name
3. `:rho`: Spearman correlation coefficient
4. `:p_value`: Unadjusted p-value
5. `:adj_p_value`: Benjamini-Hochberg adjusted p-value.
6. `:signifcant`: `true` if the FDR is significant
"""
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

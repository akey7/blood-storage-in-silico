module UfbaSamplerAnalysisAndViz

using Base.Iterators
using Statistics
using StatsBase
using CairoMakie
using AlgebraOfGraphics
using ColorSchemes
using DataFrames
using DataFramesMeta
using ProgressMeter
using HypothesisTests
using MultipleTesting

export histograms_for_reaction_v2,
    plot_all_histograms_for_reactions,
    diagnose_flux_stats,
    pivot_sampling_df_long,
    interesting_reactions_and_times

"""
    histograms_for_reaction_v2(long_sampling_df, reaction_id, reaction_string)

Plots histograms for a single reaction, with time points as separate panels and additives layered on top of each other in different colors. Draws a thick black dashed vertical line at the 0 point on all rows.

# Arguments
1. `long_sampling_df`: Sampling DataFrame, pivoted long
2. `reaction_id`: The reaction id for which the samples are being plotted.
3. `reaction_string`: The human-readable reaction string to place as a subtitle on the plot.

# Returns
`Figure`

Returns a Makie `Figure` to display or save.
"""
function histograms_for_reaction_v2(long_sampling_df, reaction_id, reaction_string)
    plt_df = @chain long_sampling_df begin
        @rsubset(:reaction_id == reaction_id)
        @rtransform(:time_span = "Week $(:final_time - 1) to $(:final_time)")
    end
    title = "$reaction_id\n$reaction_string"
    additive_palette = [
        "01-Ctrl AS3" => :dodgerblue,
        "02-Adenosine" => :orange,
        "03-Glutamine" => :blueviolet,
        "04-Methionine" => :crimson,
        "07-NAC" => :brown,
        "08-Taurine" => :magenta,
    ]
    hist_layer =
        data(plt_df) *
        mapping(:flux; color = :additive, row = :time_span => nonnumeric) *
        histogram(bins = 20) *
        visual(alpha = 0.5)
    zero_line_layer =
        data((flux = [0],)) *
        mapping(:flux) *
        visual(VLines; color = :black, linestyle = :dash, linewidth = 3)
    plt = hist_layer + zero_line_layer
    return draw(
        plt,
        scales(Color = (; palette = additive_palette));
        facet = (; linkxaxes = :all, linkyaxes = :all),
        figure = (; title = title, size = (700, 700)),
    )
end

"""
    pivot_sampling_df_long(sampling_df)

Pivots the sampling_df longer.

# Arguments
1. `sampling_df`: The sampling DataFrame in long format. The DataFrame should have a column for each reaction sampled, along with `:additive` and `:final_time` columns.

# Returns
`DataFrame`

Returns a DataFrame pivoted to long with the following columns:
1. `:additive`: The additive
2. `:final_time`: Final time
3. `:reaction_id`: The reaction id
4. `:flux`: THe flux through that reaction at that sample.
"""
function pivot_sampling_df_long(sampling_df)
    long_sampling_df = stack(
        sampling_df,
        Not([:additive, :final_time]),
        variable_name = :reaction_id,
        value_name = :flux,
    )
    return long_sampling_df
end

"""
    plot_all_histograms_for_reactions(sampling_df, rxn_ids_to_strings)

Plots version 2 of all histograms (with time points for all additives on the same figure). This function saves each figure as they are made to the `output/uFBA_histograms_v2` folder. Displays a progress meter as the plots are made.

# Arguments
1. `sampling_df`: Wide DataFrame of uFBA sampling results.

2. `rxn_ids_to_strings`: Dictionary mapping reaction ids to human readable strings for plot subtitles.
"""
function plot_all_histograms_for_reactions(sampling_df, rxn_ids_to_strings)
    if nrow(sampling_df) == 0
        @info "uFBA: Nothing to plot"
    else
        @info "uFBA: Plotting histograms, version 2"
        long_sampling_df = pivot_sampling_df_long(sampling_df)
        reaction_ids = unique(long_sampling_df.reaction_id)
        n_reaction_ids = length(reaction_ids)
        prog = Progress(n_reaction_ids, desc = "Writing histograms, version 2")
        for reaction_id in reaction_ids
            reaction_string = rxn_ids_to_strings[reaction_id]
            fig = histograms_for_reaction_v2(long_sampling_df, reaction_id, reaction_string)
            filename = joinpath("output", "uFBA_histograms_v2", "$reaction_id.png")
            save(filename, fig)
            next!(prog)
        end
        finish!(prog)
    end
end

"""
    diagnose_flux_stats(sampling_df)

Calculates diagnostic statistics for uFBA models. Expecially useful for finding fluxes that average zero flux.

# Arugments
1. `sampling_df`: Wide format DataFrame of sampling results.

# Returns
`DataFrame`

Returns a DataFrame with the following columns
1. `:additive`: Additive of the model
2. `:final_time`: Final time of the regression of rates of metabolite concentration change.
3. `:reaction_id`: Id of the flux.
4. `:mean_is_approx_zero`: True if the mean is approximately zero, false if the mean is non-zero.
5. `:mean`: Mean of the flux distribution
6. `:std`: Sample standard deviation of the flux distribution.
"""
function diagnose_flux_stats(sampling_df)
    @info "Calculating flux descriptive statistics"
    long_sampling_df = pivot_sampling_df_long(sampling_df)
    descriptions_df = @chain long_sampling_df begin
        @groupby(:additive, :final_time, :reaction_id)
        @combine(
            :mean_is_approx_zero = isapprox(mean(:flux), 0.0, atol = 1.0e-10),
            :mean = mean(:flux),
            :std = std(:flux),
        )
    end
    return descriptions_df
end

"""
    interesting_reactions_and_times(sampling_df)

Looks at all pairs of reactions and timepoints to determine if there are significant differences among the flux distributions using a k-sample Anderson-Darling test.

# Arguments
1. `sampling_df`: The wide sampling DataFrame.

# Returns
`DataFrame`

Returns a DataFrame with the following columns:
1. `:reaction_id`: The reaction id
2. `:final_time`: Final time point
3. `:p_value`: Unadjsuted p-value
4. `:adj_p_value`: Benjamini-Hochberg adjusted p-value

The resulting DataFrame is sortedin ascending p-value order.
"""
function interesting_reactions_and_times(sampling_df; n_samples = 10)
    @info "Finding interesting reactions and time points"
    long_df = pivot_sampling_df_long(sampling_df)
    reaction_ids = sort(unique(long_df.reaction_id))
    final_times = sort(unique(long_df.final_time))
    pairs = product(reaction_ids, final_times)
    n_pairs = length(pairs)
    pair_results = []
    prog = Progress(n_pairs, desc = "Evaluating reactions and time points")
    for (reaction_id, final_time) in pairs
        df = @chain long_df begin
            @rsubset(:reaction_id == reaction_id, :final_time == final_time)
            @groupby(:additive)
            transform(eachindex => :sample)
            unstack(:sample, :additive, :flux)
            select(Not(:sample))
        end
        xs = [
            sample(Float64[coalesce(x, 0.0) for x in col], n_samples; replace = false)
            for col in eachcol(df)
        ]
        ad_test = KSampleADTest(xs...)
        pv = pvalue(ad_test)
        result = (reaction_id = reaction_id, final_time = final_time, p_value = pv)
        push!(pair_results, result)
        next!(prog)
    end
    pair_results_df = DataFrame(pair_results)
    pair_results_df.adj_p_value = adjust(pair_results_df.p_value, BenjaminiHochberg())
    final_df = sort(pair_results_df, :adj_p_value)
    return final_df
end

end

module UfbaSamplerAnalysisAndViz

using Base.Iterators
using Statistics
using CairoMakie
using AlgebraOfGraphics
using ColorSchemes
using DataFrames
using DataFramesMeta
using ProgressMeter
using HypothesisTests

export histograms_for_reaction_in_additive_v1,
    histograms_for_reaction_v2,
    plot_all_histograms_v1,
    plot_all_histograms_for_reactions,
    diagnose_flux_stats,
    pivot_sampling_df_long,
    interesting_reactions_and_times

"""
    histograms_for_reaction_in_additive_v1(long_sampling_df, additive, reaction_id, reaction_string)

Plots histograms for uFBA results at all time points in a SINGLE additive on one Makie plot.

# Arguments
1. `long_sampling_df`: Sampling DataFrame, pivoted long
2. `additive`: Additive to make plots for.
3. `reaction_id`: The reaction id for which the samples are being plotted.
4. `reaction_string`: The human-readable reaction string to place as a subtitle on the plot.

# Returns
`Figure`

Returns a Makie `Figure` that can be displayed or saved.
"""
function histograms_for_reaction_in_additive_v1(
    long_sampling_df,
    additive,
    reaction_id,
    reaction_string,
)
    plt_df = @chain long_sampling_df begin
        @rsubset(:additive == additive, :reaction_id == reaction_id)
        select(:final_time, :flux)
    end
    title = "$reaction_id in $additive\n$reaction_string"
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Flux (mM/week)", ylabel = "Density", title = title)
    final_times = sort(unique(plt_df.final_time))
    colors = [:dodgerblue, :orange, :blueviolet, :crimson, :deeppink]
    for (final_time, color) in zip(final_times, colors)
        hist_df = @rsubset(plt_df, :final_time == final_time)
        hist!(
            ax,
            hist_df.flux;
            bins = 50,
            color = (color, 0.33),
            label = string(final_time),
        )
    end
    axislegend(ax)
    return fig
end

"""
    histograms_for_reaction_v2(long_sampling_df, reaction_id, reaction_string)

Plots histograms for a single reaction, with time points as separate panels and additives layered on top of each other in different colors.

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
    plt =
        data(plt_df) *
        mapping(:flux; color = :additive, row = :time_span => nonnumeric) *
        histogram(bins = 20) *
        visual(alpha = 0.5)
    return draw(
        plt,
        scales(Color = (; palette = additive_palette));
        facet = (; linkxaxes = :all, linkyaxes = :all),
        figure = (; title = title, size = (700, 700)),
    )
end

"""
    plot_all_histograms_v1(sampling_df, rxn_ids_to_strings)

Plots version 1 histograms for all reactions in all addititves (with separate figures for each additive). This function saves each figure as they are made to the `output/uFBA_histograms_v1` folder. Displays a progress meter as the plots are made.

# Arguments
1. `sampling_df`: Wide DataFrame of uFBA sampling results.

2. `rxn_ids_to_strings`: Dictionary mapping reaction ids to human readable strings for plot subtitles.
"""
function plot_all_histograms_v1(sampling_df, rxn_ids_to_strings)
    if nrow(sampling_df) == 0
        @info "uFBA: Nothing to plot"
    else
        @info "uFBA: Plotting histograms, version 1"
        long_sampling_df = stack(
            sampling_df,
            Not([:additive, :final_time]),
            variable_name = :reaction_id,
            value_name = :flux,
        )
        additives = unique(long_sampling_df.additive)
        reaction_ids = unique(long_sampling_df.reaction_id)
        pairs = product(additives, reaction_ids)
        n_pairs = length(pairs)
        prog = Progress(n_pairs, desc = "Writing histograms, version 1...")
        for (additive, reaction_id) in pairs
            reaction_string = rxn_ids_to_strings[reaction_id]
            fig = histograms_for_reaction_in_additive_v1(
                long_sampling_df,
                additive,
                reaction_id,
                reaction_string,
            )
            filename =
                joinpath("output", "uFBA_histograms_v1", "$additive $(reaction_id).png")
            save(filename, fig)
            next!(prog)
        end
        finish!(prog)
    end
end

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

function interesting_reactions_and_times(sampling_df)
    @info "Finding interesting reactions and time points"
    long_df = pivot_sampling_df_long(sampling_df)
    reaction_ids = sort(unique(long_df.reaction_id))
    final_times = sort(unique(long_df.final_time))
    pairs = collect(product(reaction_ids, final_times))[1:2]
    n_pairs = length(pairs)
    # prog = Progress(n_pairs, desc = "Evaluating reactions and time points")
    for (reaction_id, final_time) in pairs
        df = @chain long_df begin
            @rsubset(:reaction_id == reaction_id, :final_time == final_time)
            @groupby(:additive)
            transform(eachindex => :sample)
            unstack(:sample, :additive, :flux)
            select(Not(:sample))
        end
        xs = [Float64[coalesce(x, 0.0) for x in col] for col in eachcol(df)]
        result = KSampleADTest(xs...)
        display(result)
        # next!(prog)
    end
end

end

module UfbaSamplerViz

using Base.Iterators
using CairoMakie
using AlgebraOfGraphics
using ColorSchemes
using DataFrames
using DataFramesMeta
using ProgressMeter

export plot_all_histograms_v1, plot_all_histograms_for_reactions

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

function plot_all_histograms_for_reactions(sampling_df, rxn_ids_to_strings)
    if nrow(sampling_df) == 0
        @info "uFBA: Nothing to plot"
    else
        @info "uFBA: Plotting histograms, version 2"
        long_sampling_df = stack(
            sampling_df,
            Not([:additive, :final_time]),
            variable_name = :reaction_id,
            value_name = :flux,
        )
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

end

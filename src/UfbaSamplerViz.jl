module UfbaSamplerViz

using Base.Iterators
using CairoMakie
using DataFrames
using DataFramesMeta
using ThreadsX

export plot_all_histograms

function histograms_for_reaction_in_additive(long_sampling_df, additive, reaction_id)
    plt_df = @chain long_sampling_df begin
        @rsubset(:additive == additive, :reaction_id == reaction_id)
        select(:final_time, :flux)
    end
    title = "$additive $reaction_id"
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Flux (mM/week)", ylabel = "Density", title = title)
    final_times = [2, 3, 4, 5, 6]
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

function plot_all_histograms(sampling_df)
    if nrow(sampling_df) == 0
        @info "uFBA: Nothing to plot"
    else
        @info "uFBA: Plotting histograms"
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
        for (i, (additive, reaction_id)) in enumerate(pairs)
            fig =
                histograms_for_reaction_in_additive(long_sampling_df, additive, reaction_id)
            filename = joinpath("output", "uFBA_histograms", "$additive $(reaction_id).png")
            save(filename, fig)
            println("Wrote $i of $n_pairs: $filename")
        end
    end
end

end

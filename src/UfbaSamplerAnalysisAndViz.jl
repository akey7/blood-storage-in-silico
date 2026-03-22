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
using Chain

export histograms_for_reaction_v2,
    plot_all_histograms_for_reactions,
    diagnose_flux_stats,
    pivot_sampling_df_long,
    map_metabolites_to_sinks,
    net_sink_fluxes,
    net_flux_from_up_and_down,
    calc_median_flux_df,
    combine_and_clean_addititve_final_time,
    prepare_median_flux_vector_matrix,
    prepare_measurements_and_sinks_report_df


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
4. `:flux`: The flux through that reaction at that sample.
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
    map_metabolites_to_sinks(long_df, additive, final_time)

Extracts the sink from the given samples and maps unmeasured metabolite ids to their corresponding sinks.

# Arguments
1. `long_df`: The long sampling DataFrame.
2. `additive`: The additive of interest
3. `final_time`: The final time of interest

# Returns
`Dict{String,String}`

Returns a dictionary mapping strings (metabolite_ids) to their sink ids.
"""
function map_metabolites_to_sinks(long_df, additive, final_time)
    filtered_df = @rsubset(long_df, :additive == additive, :final_time == final_time)
    reaction_ids = sort(unique(filtered_df.reaction_id))
    sink_ids = [
        reaction_id for reaction_id in reaction_ids if contains(reaction_id, "R_REVSK")
    ]
    sink_metabolite_ids =
        sort(unique([join(split(sink_id, "_")[3:end], "_") for sink_id in sink_ids]))
    metabolite_id_sink_map::Dict{String,Dict{Symbol,String}} = Dict()
    for sink_metabolite_id in sink_metabolite_ids
        sink_id = "R_REVSK_$sink_metabolite_id"
        metabolite_id_sink_map[sink_metabolite_id] = sink_id
    end
    return metabolite_id_sink_map
end

"""
    calc_median_flux_df(sampling_df)

Calculate median fluxes for all reactions in the given wide sampling DataFrame.

# Arguments
1. `sampling_df`: Wide sampling DataFrame returned by [`sample_fluxes`](@ref BloodStorageInSilico.UfbaSampler.sample_fluxes)

# Returns
`DataFrame`

Returns a DataFrame with the following columns: `additive`, `final_time`, `reaction_id`, and `median_flux` columns.
"""
function calc_median_flux_df(sampling_df)
    long_df = pivot_sampling_df_long(sampling_df)
    median_df = @chain long_df begin
        @groupby(:additive, :final_time, :reaction_id)
        @combine(:median_flux = median(:flux))
    end
    return median_df
end

"""
    combine_and_clean_addititve_final_time(additive, final_time)

Combine additive and final time specifications into a single lowercase string with `-` and ` ` substituted with `_`.

# Arguments
1. `additive`: The additive string
2. `final_time`: The final time integer

# Returns
`String`

Returns a string formatted in the way specified above.
"""
function combine_and_clean_addititve_final_time(additive, final_time)
    cleaned_additive = @chain additive begin
        lowercase()
        replace("-" => "_", " " => "_")
    end
    combined = "$(cleaned_additive)_$(final_time)"
    return combined
end

"""
    prepare_median_flux_vector_matrix(sampling_df)

Prepare a data matrix of the uFBA results. Each row is a reaction, each column is an additive at a time point, and each element is the median flux for that row and column.

# Arguments
1. `sampling_df`: The wide formatted sampling DataFrame

# Returns
`DataFrame`

Returns a data matrix in the form of a DataFrame as specified above.
"""
function prepare_median_flux_vector_matrix(sampling_df)
    median_df = calc_median_flux_df(sampling_df)
    transformed_df = @chain median_df begin
        @rtransform(
            :additive_final_time =
                combine_and_clean_addititve_final_time(:additive, :final_time)
        )
        @select(:additive_final_time, :reaction_id, :median_flux)
        @orderby(:additive_final_time, :reaction_id)
        unstack(:reaction_id, :additive_final_time, :median_flux)
    end
    return transformed_df
end

"""
    prepare_measurements_and_sinks_report_df(absolute_quant_long_df, fba_model_metabolites_df, ufba_added_sink_ids_df, sampling_df)

Prepares two DataFrames detailing the number of metabolites in each model and how many sinks there are.

The first DataFrame is longer than the second and contains one row per metabolite. On each row, the additive, time, and metabolite are listed, along with boolean columns with whether the metabolite has up/down sinks, or is measured.

The second DataFrame aggregates these per-metabolite rows into per-model rows, outlining overal counts of unique metabolites, how many metabolites are measured, and the number of up and down sinks.

# Arguments
1. `absolute_quant_long_df`: The absolute quant measurements from `output/absolute_quant_long.csv`
2. `fba_model_metabolites_df`: The FBA model metabolites from `output/fba_model_metabolites.csv`
3. `ufba_added_sink_ids_df`: The sinks added to models during the uFBA process from `output/ufba_added_sink_ids.csv`
4. `sampling_df`: The wide sampling DataFrame from all uFBA runs.

# Returns
`Tuple{DataFrame,DataFrame}`

Returns two DataFrames, one for each report. The first element is the per-metabolite report, and the second element is the per-model report.

The first DataFrame contains the following columns
1. `additive`: Additive of the model
2. `final_time`: Final time of the model
3. `fba_metabolite_id`: Metabolite id from the FBA model
4. `is_measured`: `true` if the metabolite was measured
5. `has_up_sink`: `true` if there is an up sink for the metabolite
6. `has_down_sink`: `true` if there is a down sink for the metabolite

The second DataFrame contains the following columns:
1. `additive`: Additive of the model
2. `final_time`: Final time of the model
3. `n_fba_metabolites`: Number of metabolites in the model
4. `n_measured`: The number of metabolites that have absolute quant approximations.
5. `n_at_least_one_sink`: The number of metabolites that have one or both sinks.
6. `n_no_measure_no_sink`: The number of metabolites that have neither an absolute quant approximation nor any sinks.
"""
function prepare_measurements_and_sinks_report_df(
    absolute_quant_long_df,
    fba_model_metabolites_df,
    ufba_optimized_sinks_df,
    sampling_df,
)
    long_sampling_df = pivot_sampling_df_long(sampling_df)
    additives = sort(unique(long_sampling_df.additive))
    final_times = sort(unique(long_sampling_df.final_time))
    fba_metabolite_ids = sort(unique(fba_model_metabolites_df.metabolite_id))
    tasks = product(fba_metabolite_ids, additives, final_times)
    n_tasks = length(tasks)
    prog = Progress(n_tasks, "Matching FBA, measured, and sink metabolite ids")
    rows = map(tasks) do t
        fba_metabolite_id, additive, final_time = t
        measured_df = @rsubset(
            absolute_quant_long_df,
            :Additive == additive,
            :Time == final_time,
            :Metabolite == fba_metabolite_id
        )
        sink_df = @rsubset(
            ufba_optimized_sinks_df,
            :additive == additive,
            :final_time == final_time,
            :metabolite_id == fba_metabolite_id,
        )
        is_measured = nrow(measured_df) > 0
        has_sink = nrow(sink_df) > 0
        next!(prog)
        return (
            additive = additive,
            final_time = final_time,
            fba_metabolite_id = fba_metabolite_id,
            is_measured = is_measured,
            has_sink = has_sink,
        )
    end
    unsorted_report_df = DataFrame(rows)
    report_df = @orderby(unsorted_report_df, :additive, :final_time, :fba_metabolite_id)
    report_by_model_df = @chain report_df begin
        @rtransform(:no_measure_no_sink = !:is_measured && !:has_sink)
        @groupby(:additive, :final_time)
        @combine(
            :n_fba_metabolites = length(unique(:fba_metabolite_id)),
            :n_measured = sum(:is_measured),
            :n_with_sink = sum(:has_sink),
            :n_no_measure_no_sink = sum(:no_measure_no_sink),
        )
        @orderby(:additive, :final_time)
    end
    return report_df, report_by_model_df
end

end

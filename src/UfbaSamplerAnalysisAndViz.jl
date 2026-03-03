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
    prepare_measurements_and_sinks_report_df,
    safely_query_sink_map


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

Extracts the sinks (up and down) from the given samples and maps unmeasured metabolite ids to their corresponding up and down sinks.

# Arguments
1. `long_df`: The long sampling DataFrame.
2. `additive`: The additive of interest
3. `final_time`: The final time of interest

# Returns
`Dict{String,Dict{Symbol,String}}`

Returns a dictionary mapping strings (metabolite_ids) to a second level of dictionaries. The second level of dictionaries contain `:up` and/or `:down` keys which in turn map to reaction ids that are the up or and/or down sinks for the metabolite_id key in the top-level dictionary.
"""
function map_metabolites_to_sinks(long_df, additive, final_time)
    filtered_df = @rsubset(long_df, :additive == additive, :final_time == final_time)
    reaction_ids = sort(unique(filtered_df.reaction_id))
    sink_ids = [
        reaction_id for reaction_id in reaction_ids if contains(reaction_id, "R_UNKNOWN_SK")
    ]
    sink_metabolite_ids =
        sort(unique([join(split(sink_id, "_")[5:end], "_") for sink_id in sink_ids]))
    metabolite_id_sink_map::Dict{String,Dict{Symbol,String}} = Dict()
    for sink_metabolite_id in sink_metabolite_ids
        sink_up_id = "R_UNKNOWN_SK_UP_$sink_metabolite_id"
        sink_down_id = "R_UNKNOWN_SK_DOWN_$sink_metabolite_id"
        metabolite_id_sink_map[sink_metabolite_id] = Dict()
        if sink_up_id in sink_ids
            metabolite_id_sink_map[sink_metabolite_id][:up] = sink_up_id
        end
        if sink_down_id in sink_ids
            metabolite_id_sink_map[sink_metabolite_id][:down] = sink_down_id
        end
    end
    return metabolite_id_sink_map
end

"""
    net_flux_from_up_and_down(up_flux::Union{Float64,Missing}, down_flux::Union{Float64})

Helper function for [`net_sink_fluxes`](@ref BloodStorageInSilico.UfbaSamplerAnalysisAndViz.net_sink_fluxes). Calculates the net flux of two sink fluxes avoiding missing values.

# Arguments
1. `up_flux::Union{Float64,Missing}`: Flux of the up sink. If `missing`, this value is ignored when computing the net flux.
2. `down_flux::Union{Float64}`: Flux of the down sink. If `missing`, this value is ignored when computing the net flux.

# Returns
`Float64`

Returns the net flux of the sinks, calculated by summing the fluxes together and skipping missing values.
"""
function net_flux_from_up_and_down(
    up_flux::Union{Float64,Missing},
    down_flux::Union{Float64},
)
    if !ismissing(up_flux) && !ismissing(down_flux)
        return up_flux + down_flux
    elseif !ismissing(up_flux)
        return up_flux
    else
        return down_flux
    end
end

"""
    net_sink_fluxes(sampling_df)

Calucates the net fluxes between each pair of sinks by summing their values together (when both sinks are present) or selecting only the up or down flux where just one sink is available.

# Arguments
1. `sampling_df`: Wide DataFrame of sampling values.

# Returns
`DataFrame`

Returns a DataFrame with the following columns
1. `additive`: The additive.
2. `final_time`: The final time point of the model
3. `metabolite_id`: The metabolite id matching the sinks.
4. `up_median_flux`: The median flux of the up sink flux distribution.
5. `down_median_flux`: The median flux of the down sink flux distribution.
6. `net_median_flux`: The net flux summed over both sinks.
"""
function net_sink_fluxes(sampling_df)
    @info "Calculating net sink fluxes"
    long_df = pivot_sampling_df_long(sampling_df)
    final_times = sort(unique(long_df.final_time))
    additives = sort(unique(long_df.additive))
    pairs = product(additives, final_times)
    median_fluxes_df = @chain long_df begin
        @groupby(:additive, :final_time, :reaction_id)
        @combine(:median_flux = median(:flux))
    end
    rows = []
    n_calculations = length(pairs)
    prog = Progress(n_calculations, desc = "Calculating net sink fluxes")
    for (additive, final_time) in pairs
        sink_map = map_metabolites_to_sinks(long_df, additive, final_time)
        for (metabolite_id, sinks) in sink_map
            up_id = get(sinks, :up, nothing)
            down_id = get(sinks, :down, nothing)
            up_df =
                !isnothing(up_id) ?
                @rsubset(
                    median_fluxes_df,
                    :additive == additive,
                    :final_time == final_time,
                    :reaction_id == up_id
                ) : nothing
            down_df =
                !isnothing(down_id) ?
                @rsubset(
                    median_fluxes_df,
                    :additive == additive,
                    :final_time == final_time,
                    :reaction_id == down_id
                ) : nothing
            up_median_flux = !isnothing(up_df) ? up_df[1, :median_flux] : missing
            down_median_flux = !isnothing(down_df) ? down_df[1, :median_flux] : missing
            net_median_flux = net_flux_from_up_and_down(up_median_flux, down_median_flux)
            row = (
                metabolite_id = metabolite_id,
                additive = additive,
                final_time = final_time,
                up_median_flux = up_median_flux,
                down_median_flux = down_median_flux,
                net_median_flux = net_median_flux,
            )
            push!(rows, row)
            next!(prog)
        end
    end
    net_flux_df_1 = DataFrame(rows)
    net_flux_df_2 = @chain net_flux_df_1 begin
        @select(
            :additive,
            :final_time,
            :metabolite_id,
            :up_median_flux,
            :down_median_flux,
            :net_median_flux
        )
        @orderby(:additive, :final_time, :metabolite_id)
    end
    return net_flux_df_2
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

function safely_query_sink_map(sink_map, metabolite_id, direction)
    metabolite_sinks = get(sink_map, metabolite_id, nothing)
    if isnothing(metabolite_sinks)
        return missing
    else
        return get(metabolite_sinks, direction, missing)
    end
end

function prepare_measurements_and_sinks_report_df(
    absolute_quant_long_df,
    fba_model_metabolites_df,
    sampling_df,
)
    long_sampling_df = pivot_sampling_df_long(sampling_df)
    additives = sort(unique(long_sampling_df.additive))
    final_times = sort(unique(long_sampling_df.final_time))
    fba_model_metabolite_ids = sort(unique(fba_model_metabolites_df.metabolite_id))
    measured_metabolite_ids = sort(unique(absolute_quant_long_df.Metabolite))
    pairs = product(additives, final_times)
    n_pairs = length(pairs)
    prog = Progress(n_pairs, "Preparing measurements and sinks report")
    report_rows = []
    for (additive, final_time) in pairs
        sink_map = map_metabolites_to_sinks(long_sampling_df, additive, final_time)
        for fba_model_metabolite_id in fba_model_metabolite_ids
            is_measured = fba_model_metabolite_id in measured_metabolite_ids
            has_sinks = fba_model_metabolite_id in keys(sink_map)
            up_sink = safely_query_sink_map(sink_map, fba_model_metabolite_id, :up)
            down_sink = safely_query_sink_map(sink_map, fba_model_metabolite_id, :down)
            report_row = (
                additive = additive,
                final_time = final_time,
                metabolite_id = fba_model_metabolite_id,
                is_measured = is_measured,
                has_sinks = has_sinks,
                up_sink = up_sink,
                down_sink = down_sink,
            )
            push!(report_rows, report_row)
        end
        next!(prog)
    end
    report_df = DataFrame(report_rows)
    report_aggregated_df = @chain report_df begin
        @groupby(:additive, :final_time)
        @combine(:n_measured = sum(:is_measured), :n_have_sinks = sum(:has_sinks))
        @orderby(:additive, :final_time)
    end
    return report_df, report_aggregated_df
end

end

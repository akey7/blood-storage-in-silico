module UfbaSamplerViz

using Base.Iterators
using Random
using CSV
using DataFrames
using DataFramesMeta
using CategoricalArrays
using Chain
using StatsBase
using StatsModels
using Clustering
using MultivariateStats
using Statistics
using PlotlyJS
using KernelDensity
using CairoMakie
using AlgebraOfGraphics
using ColorSchemes
using ProgressMeter

export stacked_flux_kde_3d,
    plot_all_distributions_for_reactions,
    densities_for_reaction,
    histograms_for_reaction_v2,
    load_sampling_results,
    pca_treatment_effects,
    prepare_treatment_effects_dfs,
    k_means_treatment_effects,
    plot_treatment_effects_kmeans_pca,
    prepare_median_fluxes_dfs,
    k_means_median_fluxes

function load_sampling_results()
    sampling_filename = joinpath("output", "ufba_sampling.csv")
    sampling_df = CSV.read(sampling_filename, DataFrame)
    additives = sort(unique(sampling_df.additive))
    final_times = sort(unique(sampling_df.final_time))
    all_possible = product(additives, final_times)
    working_model_rows = []
    for (additive, final_time) in all_possible
        trial_df = @rsubset(sampling_df, :additive == additive, :final_time == final_time)
        if nrow(trial_df) > 0
            row = (additive = additive, final_time = final_time)
            push!(working_model_rows, row)
        end
    end
    working_models_df = @chain working_model_rows begin
        DataFrame()
        @orderby(:additive, :final_time)
    end
    result = (sampling_df = sampling_df, working_models_df = working_models_df)
    return result
end

"""
    pivot_sampling_df_long(sampling_df)

Pivots the sampling DataFrame longer. This function duplicates the function of the same name in `UfbaSamplerVisualization`. It is duplicated so that neither this module nor the visualization module need to import each other, which would mess up the documentation.

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
    stacked_flux_kde_3d(
        df::DataFrame;
        z_spacing::Real = 0.35,
        npoints::Int = 256,
        bandwidth = nothing,
        line_width::Real = 5,
        peak_projection_width::Real = 2,
        plane_opacity::Real = 0.10,
        fig_title::AbstractString = "3D KDE Curves Across Weeks",
    )

Makes a 3D plot of a given flux distribution trajectory through time. Kernel density estimation are plotted up the z axis and a line connecting the mode bins of all the KDEs is plotted behind these histograms to show the trajectory of flux over time. Lines from the peak of each curve to the trajectory line make the plot more readable.

# Arguments
1. `df::DataFrame`: Long format DataFrame of flux sampling results.
2. `z_spacing::Real = 0.35`: How far apart each histogram is up the z axis.
3. `npoints::Int = 256`: Number of points on the KDE curve
4. `bandwidth = nothing`: Bandwidth of KDE, if specified.
5. `line_width::Real = 5`: Line width of the plot.
6. `peak_projection_width::Real = 2`: Width of line between peaks of curves and the mode trajectory line.
7. `plane_opacity::Real = 0.10`: Opacity of reference planes.
8. `fig_title::AbstractString = "3D KDE Curves Across Weeks"`: Figure title if specified.

# Returns
`PlotlyJS.Plot`

A PlotlyJS plot to display or save.
"""
function stacked_flux_kde_3d(
    df::DataFrame;
    z_spacing::Real = 0.35,
    npoints::Int = 256,
    bandwidth = nothing,
    line_width::Real = 5,
    peak_projection_width::Real = 2,
    plane_opacity::Real = 0.10,
    fig_title::AbstractString = "3D KDE Curves Across Weeks",
)
    weeks = sort(unique(df.final_time))
    length(weeks) > 0 || error("No weeks found in flux_df")
    z_spacing > 0 || error("z_spacing must be positive")
    npoints >= 32 || error("npoints should be at least 32 for a smooth KDE curve")
    all_flux = Float64.(df.flux)
    xmin = minimum(all_flux)
    xmax = maximum(all_flux)
    if isapprox(xmin, xmax)
        δ = max(abs(xmin) * 0.05, 1e-6)
        xmin -= δ
        xmax += δ
    end
    week_z = Dict(week => (i - 1) * z_spacing for (i, week) in enumerate(weeks))
    zvals = [week_z[w] for w in weeks]
    palette = [
        "#1f77b4",
        "#d62728",
        "#2ca02c",
        "#ff7f0e",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
    ]
    week_colors =
        Dict(week => palette[mod1(i, length(palette))] for (i, week) in enumerate(weeks))
    kdes = Dict{eltype(weeks),Any}()
    ymax = 0.0

    for week in weeks
        vals = Float64.(df[df.final_time .== week, :flux])

        kd = if bandwidth === nothing
            kde(vals; boundary = (xmin, xmax), npoints = npoints)
        else
            kde(vals; boundary = (xmin, xmax), npoints = npoints, bandwidth = bandwidth)
        end

        kdes[week] = kd
        ymax = max(ymax, maximum(kd.density))
    end
    zmin = minimum(zvals) - 0.12
    zmax = maximum(zvals) + 0.12
    traces = GenericTrace[]
    push!(
        traces,
        PlotlyJS.surface(
            x = [xmin, xmax],
            y = [0.0, ymax],
            z = fill(zmin, 2, 2),
            colorscale = [[0.0, "white"], [1.0, "white"]],
            opacity = plane_opacity,
            showscale = false,
            hoverinfo = "skip",
            showlegend = false,
            name = "XY plane",
        ),
    )
    push!(
        traces,
        PlotlyJS.surface(
            x = [xmin, xmax],
            y = fill(0.0, 2, 2),
            z = [zmin zmin; zmax zmax],
            colorscale = [[0.0, "white"], [1.0, "white"]],
            opacity = plane_opacity,
            showscale = false,
            hoverinfo = "skip",
            showlegend = false,
            name = "XZ plane",
        ),
    )
    peak_x = Float64[]
    peak_y = Float64[]
    peak_z = Float64[]
    for week in weeks
        kd = kdes[week]
        zpos = week_z[week]
        x_curve = Float64.(kd.x)
        y_curve = Float64.(kd.density)
        z_curve = fill(zpos, length(x_curve))
        push!(
            traces,
            scatter3d(
                x = x_curve,
                y = y_curve,
                z = z_curve,
                mode = "lines",
                name = "Week $week",
                line = attr(width = line_width, color = week_colors[week]),
                hovertemplate = "Week $week<br>Flux: %{x:.4f}<br>Density: %{y:.4f}<extra></extra>",
            ),
        )
        peak_idx = argmax(y_curve)
        px = x_curve[peak_idx]
        py = y_curve[peak_idx]
        push!(peak_x, px)
        push!(peak_y, 0.0)
        push!(peak_z, zpos)
        push!(
            traces,
            scatter3d(
                x = [px, px],
                y = [0.0, py],
                z = [zpos, zpos],
                mode = "lines",
                name = "Peak projection",
                showlegend = false,
                line = attr(
                    width = peak_projection_width,
                    color = week_colors[week],
                    dash = "dot",
                ),
                hoverinfo = "skip",
            ),
        )
    end
    push!(
        traces,
        scatter3d(
            x = peak_x,
            y = peak_y,
            z = peak_z,
            mode = "lines+markers",
            name = "KDE peak trajectory",
            line = attr(width = 6, color = "black"),
            marker = attr(size = 5, color = "black"),
            hovertemplate = "Week: %{text}<br>Peak flux: %{x:.4f}<extra></extra>",
            text = string.(weeks),
        ),
    )
    layout = Layout(
        title = fig_title,
        paper_bgcolor = "white",
        plot_bgcolor = "white",
        scene = attr(
            xaxis = attr(
                title = "Flux",
                backgroundcolor = "white",
                gridcolor = "lightgray",
                zerolinecolor = "lightgray",
            ),
            yaxis = attr(
                title = "Density",
                backgroundcolor = "white",
                gridcolor = "lightgray",
                zerolinecolor = "lightgray",
            ),
            zaxis = attr(
                title = "Week",
                tickmode = "array",
                tickvals = zvals,
                ticktext = ["Week $w" for w in weeks],
                backgroundcolor = "white",
                gridcolor = "lightgray",
                zerolinecolor = "lightgray",
            ),
            aspectmode = "manual",
            aspectratio = attr(x = 1.6, y = 1.0, z = 0.7),
            camera = attr(eye = attr(x = 1.7, y = 1.4, z = 1.1)),
        ),
        showlegend = true,
    )
    return PlotlyJS.plot(traces, layout)
end

"""
    plot_all_distributions_for_reactions(sampling_df, rxn_ids_to_strings; bins = 20)

Plots histograms and densities for all reactions in all additives at all time points. This function saves each figure as they are made to `output/uFBA_histograms_v2` or `output/uFBA_densities` as appropriate. Displays a progress meter as the plots are made. Returns nothing since it saves as it goes.

# Arguments
1. `sampling_df`: Wide DataFrame of uFBA sampling results.
2. `rxn_ids_to_strings`: Dictionary mapping reaction ids to human readable strings for plot subtitles.
3. `control_additive = nothing`: If both this and `treatment_additive` are specified, the densities/histograms will just be between these two conditions and the color palette will be fixed.
4. `treatment_additive = nothing`: See `control_additive` above.
5. `bins`: Number of bins to put onto histograms.
"""
function plot_all_distributions_for_reactions(
    sampling_df,
    rxn_ids_to_strings;
    control_additive = nothing,
    treatment_additive = nothing,
    bins = 20,
)
    if nrow(sampling_df) == 0
        @warn "uFBA: Nothing to plot"
    else
        long_sampling_df = pivot_sampling_df_long(sampling_df)
        reaction_ids = unique(long_sampling_df.reaction_id)
        n_reaction_ids = length(reaction_ids)
        prog = Progress(n_reaction_ids, desc = "Writing histograms and densities")
        for reaction_id in reaction_ids
            reaction_string = rxn_ids_to_strings[reaction_id]["rxn_string"]
            subsystem = rxn_ids_to_strings[reaction_id]["subsystem"]
            # reaction_name = rxn_ids_to_strings[reaction_id]["name"]
            fig_hist = histograms_for_reaction_v2(
                long_sampling_df,
                reaction_id,
                reaction_string,
                subsystem;
                control_additive = control_additive,
                treatment_additive = treatment_additive,
                bins = bins,
            )
            fig_density = densities_for_reaction(
                long_sampling_df,
                reaction_id,
                reaction_string,
                subsystem;
                control_additive = control_additive,
                treatment_additive = treatment_additive,
            )
            filename_hist = joinpath(
                "output",
                "uFBA_histograms_v2",
                "$treatment_additive $reaction_id Histograms.png",
            )
            filename_density = joinpath(
                "output",
                "uFBA_densities",
                "$treatment_additive $reaction_id Densities.png",
            )
            save(filename_hist, fig_hist)
            save(filename_density, fig_density)
            next!(prog)
        end
        finish!(prog)
    end
end

"""
    additive_comparison_palette(
        control_additive = nothing,
        treatment_additive = nothing,
    )

Generates the color palette for the histograms and density plots. If **either** parameter is `nothing`, defaults to the color palette for the additives in the first relative quant dataset. If string values are supplied to both parameters, sets consistent color for control and additive treatments.

# Arguments
1. `control_additive = nothing`: If specified, the name of the control additive.
2. `treatment_additive = nothing`: If specified, the name of the treatment additive.

# Returns
`Dict{String,Union{Symbol,String}}`

Returns a dictionary mapping the additive names as keys to their colors as values. The values can be either symbols for Makie or hex codes.
"""
function additive_comparison_palette(
    control_additive = nothing,
    treatment_additive = nothing,
)
    if isnothing(control_additive) || isnothing(treatment_additive)
        additive_palette = [
            "01-Ctrl AS3" => :dodgerblue,
            "02-Adenosine" => :orange,
            "03-Glutamine" => :blueviolet,
            "04-Methionine" => :crimson,
            "07-NAC" => :brown,
            "08-Taurine" => :magenta,
        ]
        return additive_palette
    else
        additive_palette = [control_additive => :dodgerblue, treatment_additive => :orange]
    end
end

"""
    subset_long_sampling_df(
        long_sampling_df,
        reaction_id;
        control_additive = nothing,
        treatment_additive = nothing,
    )
    
Subset the long sampling DataFrame for plotting the given reaction. For the second relative quant dataset, the option of restricting to the given control and treatment additives is available.

Note: **Both** control and treatment additives must be specified to filter down the additives.

# Arguments
1. `long_sampling_df`: Long sampling DataFrame
2. `reaction_id`: Reaction id of interest.
3. `control_additive = nothing`: If specified, this is the control additive to include.
4. `treatment_additive = nothing`: If specified, this is the treatment additive to include.

# Returns
`DataFrame`

Returns a DataFrame subsetted as specified.
"""
function subset_long_sampling_df(
    long_sampling_df,
    reaction_id;
    control_additive = nothing,
    treatment_additive = nothing,
)
    if !isnothing(control_additive) && !isnothing(treatment_additive)
        interesting_additives = [control_additive, treatment_additive]
        plt_df = @chain long_sampling_df begin
            @rsubset(:reaction_id == reaction_id, :additive in interesting_additives)
            @rtransform(:time_span = "Final Time $(:final_time)")
        end
        if nrow(plt_df) == 0
            @warn "Empty plt_df for $reaction_id"
        end
        return plt_df
    else
        plt_df = @chain long_sampling_df begin
            @rsubset(:reaction_id == reaction_id)
            @rtransform(:time_span = "Timespan Ending at $(:final_time)")
        end
        return plt_df
    end
end

"""
    densities_for_reaction(long_sampling_df, reaction_id, reaction_string, subsystem)

Plots KDEs of the flux distributions for the reaction in the various additives.

# Arguments
1. `long_sampling_df`: Sampling DataFrame, pivoted long
2. `reaction_id`: The reaction id for which the samples are being plotted.
3. `reaction_string`: The human-readable reaction string to place as a subtitle on the plot.
4. `subsystem`: Human-readable susbsytem of the reaction
5. `control_additive = nothing`: If both this and `treatment_additive` are specified, the densities/histograms will just be between these two conditions and the color palette will be fixed.
6. `treatment_additive = nothing`: See `control_additive` above.

# Returns
`Figure`

Returns a Makie `Figure` to display or save.
"""
function densities_for_reaction(
    long_sampling_df,
    reaction_id,
    reaction_string,
    subsystem;
    control_additive = nothing,
    treatment_additive = nothing,
)
    plt_df = subset_long_sampling_df(
        long_sampling_df,
        reaction_id;
        control_additive = control_additive,
        treatment_additive = treatment_additive,
    )
    clean_reaction_id = replace(reaction_id, "R_" => "")
    title = "$clean_reaction_id ($subsystem)\n$reaction_string"
    additive_palette = additive_comparison_palette(control_additive, treatment_additive)
    density_layer =
        data(plt_df) *
        mapping(:flux; color = :additive, row = :time_span => nonnumeric) *
        AlgebraOfGraphics.density() *
        visual(alpha = 0.5)
    zero_line_layer =
        data((flux = [0],)) *
        mapping(:flux) *
        visual(VLines; color = :black, linestyle = :dash, linewidth = 3)
    plt = density_layer + zero_line_layer
    return draw(
        plt,
        scales(
            Color = (; palette = additive_palette),
            X = (; label = "Flux (mM/week)"),
            Y = (; label = "Density"),
        );
        facet = (; linkxaxes = :all, linkyaxes = :all),
        figure = (; title = title, size = (700, 700)),
    )
end

"""
    histograms_for_reaction_v2(long_sampling_df, reaction_id, reaction_string; bins = 20)

Plots histograms for a single reaction, with time points as separate panels and additives layered on top of each other in different colors. Draws a thick black dashed vertical line at the 0 point on all rows.

# Arguments
1. `long_sampling_df`: Sampling DataFrame, pivoted long
2. `reaction_id`: The reaction id for which the samples are being plotted.
3. `reaction_string`: The human-readable reaction string to place as a subtitle on the plot.
4. `subsystem`: Human-readable susbsytem of the reaction
5. `control_additive = nothing`: If both this and `treatment_additive` are specified, the densities/histograms will just be between these two conditions and the color palette will be fixed.
6. `treatment_additive = nothing`: See `control_additive` above.
7. `bins`: Number of bins in the histograms.

# Returns
`Figure`

Returns a Makie `Figure` to display or save.
"""
function histograms_for_reaction_v2(
    long_sampling_df,
    reaction_id,
    reaction_string,
    subsystem;
    control_additive = nothing,
    treatment_additive = nothing,
    bins = 20,
)
    plt_df = subset_long_sampling_df(
        long_sampling_df,
        reaction_id;
        control_additive = control_additive,
        treatment_additive = treatment_additive,
    )
    clean_reaction_id = replace(reaction_id, "R_" => "")
    title = "$clean_reaction_id ($subsystem)\n$reaction_string"
    additive_palette = additive_comparison_palette(control_additive, treatment_additive)
    hist_layer =
        data(plt_df) *
        mapping(:flux; color = :additive, row = :time_span => nonnumeric) *
        AlgebraOfGraphics.histogram(bins = bins) *
        visual(alpha = 0.5)
    zero_line_layer =
        data((flux = [0],)) *
        mapping(:flux) *
        visual(VLines; color = :black, linestyle = :dash, linewidth = 3)
    plt = hist_layer + zero_line_layer
    return draw(
        plt,
        scales(
            Color = (; palette = additive_palette),
            X = (; label = "Flux (mM/week)"),
            Y = (; label = "Sample Count"),
        );
        facet = (; linkxaxes = :all, linkyaxes = :all),
        figure = (; title = title, size = (700, 700)),
    )
end

#####################################################################
# EFFECT SIZE PCA/K-MEANS                                           #
#####################################################################

function prepare_treatment_effects_dfs(control_vs_treatments_signif_df)
    final_times = sort(unique(control_vs_treatments_signif_df.final_time))
    effects_dfs = DataFrame[]
    for final_time in final_times
        df = @chain control_vs_treatments_signif_df begin
            @rsubset(:final_time == final_time)
            unstack(
                [:treatment_additive, :final_time],
                :reaction_id,
                :reaction_cohen_effect_z,
            )
            DataFrames.transform(
                Not([:treatment_additive, :final_time]) .=>
                    (x -> coalesce.(x, median(collect(skipmissing(x))))) .=> identity,
            )
        end
        push!(effects_dfs, df)
    end
    result = (final_times = final_times, effects_dfs = effects_dfs)
    return result
end

function k_means_treatment_effects(
    prepared_treatments_result;
    k = 5,
    seed = 123,
    maxiter = 300,
)
    final_times = prepared_treatments_result.final_times
    effects_dfs = prepared_treatments_result.effects_dfs
    feature_cols = names(effects_dfs[1], Not([:treatment_additive, :final_time]))
    cluster_dfs = DataFrame[]
    for (final_time, effects_df) in zip(final_times, effects_dfs)
        X = Matrix{Float64}(effects_df[:, feature_cols])'
        Random.seed!(seed)
        result = kmeans(X, k; maxiter = maxiter, tol = 1.0e-6, display = :none)
        if !result.converged
            @warn "k-means did not converge" k=k seed=seed maxiter=maxiter
        end
        cluster_df = DataFrame(
            treatment_additive = effects_df.treatment_additive,
            final_time = fill(final_time, nrow(effects_df)),
            n_clusters = fill(k, nrow(effects_df)),
            cluster = result.assignments,
        )
        push!(cluster_dfs, cluster_df)
    end
    treatment_k_means_df = @orderby(vcat(cluster_dfs...), :treatment_additive, :final_time)
    return treatment_k_means_df
end

function pca_treatment_effects(prepared_treatments_result; n_pcs = 5)
    final_times = prepared_treatments_result.final_times
    effects_dfs = prepared_treatments_result.effects_dfs
    pc_names = Symbol.("PC", 1:n_pcs)
    pca_dfs = DataFrame[]
    loadings_dfs = DataFrame[]
    for (final_time, df) in zip(final_times, effects_dfs)
        reaction_ids = names(select(df, Not([:treatment_additive, :final_time])))
        X = Matrix(select(df, Not([:treatment_additive, :final_time])))
        Xt = copy(X')
        M = fit(PCA, Xt; maxoutdim = n_pcs, mean = false)
        scores = MultivariateStats.transform(M, Xt)
        pca_df = DataFrame(collect(scores'), pc_names)
        insertcols!(
            pca_df,
            1,
            :treatment_additive => df.treatment_additive,
            :final_time => fill(final_time, nrow(pca_df)),
        )
        push!(pca_dfs, pca_df)
        loadings = projection(M)
        loadings_df = DataFrame(collect(loadings), pc_names)
        insertcols!(
            loadings_df,
            1,
            :reaction_id => reaction_ids,
            :final_time => fill(final_time, nrow(loadings_df)),
        )
        push!(loadings_dfs, loadings_df)
    end
    pca_df = @orderby(vcat(pca_dfs...), :treatment_additive, :final_time)
    loadings_df = @orderby(vcat(loadings_dfs...), :reaction_id, :final_time)
    result = (pca_df = pca_df, loadings_df = loadings_df)
    return result
end

function plot_treatment_effects_kmeans_pca(pca_df, treatment_k_means_df)
    n_clusters = maximum(treatment_k_means_df.cluster)
    all_time_df = @chain pca_df begin
        innerjoin(treatment_k_means_df; on = [:treatment_additive, :final_time])
        @transform(:cluster = categorical(:cluster))
        @orderby(:final_time, :cluster)
        @select(:final_time, :cluster, :PC1, :PC2)
    end
    cluster_colors = get(colorschemes[:okabe_ito], range(0, 1, length = n_clusters))
    final_times = sort(unique(all_time_df.final_time))
    n_plots = length(final_times)
    prog = Progress(n_plots, "Writing effects k-means PCA plots")
    for final_time in final_times
        filename = joinpath(
            "output",
            "viz_effects_kmeans_pca",
            "effects_kmeans_pca_$(final_time).png",
        )
        title = "Effects K-Means PCA $final_time"
        plt_df = @rsubset(all_time_df, :final_time == final_time)
        scatter_plt =
            data(plt_df) *
            mapping(:PC1, :PC2, color = :cluster) *
            visual(Scatter, markersize = 14, alpha = 0.75)
        fig = draw(
            scatter_plt,
            scales(Color = (; palette = cluster_colors)),
            figure = (; size = (500, 500)),
            axis = (; title = title),
        )
        save(filename, fig)
        next!(prog)
    end
end

#####################################################################
# MEDIAN FLUX PCA/K-MEANS                                           #
#####################################################################

function zscore_col(xs)
    μ = mean(skipmissing(xs))
    σ = std(skipmissing(xs))
    if isapprox(σ, 0.0)
        return fill(0.0, length(xs))
    else
        return (xs .- μ) ./ σ
    end
end

function prepare_median_fluxes_dfs(median_fluxes_df)
    final_times = sort(unique(median_fluxes_df.final_time))
    centered_scaled_dfs = DataFrame[]
    for final_time in final_times
        df = @chain median_fluxes_df begin
            @rsubset(:final_time == final_time)
            unstack([:additive, :final_time], :reaction_id, :median_flux)
        end
        feature_cols = names(df, Not([:additive, :final_time]))
        for col in feature_cols
            @transform!(df, $col = zscore_col($col))
        end
        push!(centered_scaled_dfs, df)
    end
    result = (final_times = final_times, centered_scaled_dfs = centered_scaled_dfs)
    return result
end

function k_means_median_fluxes(prepared_medians_result; k = 5, seed = 123, maxiter = 300)
    final_times = prepared_medians_result.final_times
    centered_scaled_dfs = prepared_medians_result.centered_scaled_dfs
    feature_cols = names(centered_scaled_dfs[1], Not([:additive, :final_time]))
    cluster_dfs = DataFrame[]
    for (final_time, median_fluxes_df) in zip(final_times, centered_scaled_dfs)
        X = Matrix{Float64}(median_fluxes_df[:, feature_cols])'
        Random.seed!(seed)
        result = kmeans(X, k; maxiter = maxiter, tol = 1.0e-6, display = :none)
        if !result.converged
            @warn "k-means did not converge" k=k seed=seed maxiter=maxiter
        end
        cluster_df = DataFrame(
            additive = median_fluxes_df.additive,
            final_time = fill(final_time, nrow(median_fluxes_df)),
            n_clusters = fill(k, nrow(median_fluxes_df)),
            cluster = result.assignments,
        )
        push!(cluster_dfs, cluster_df)
    end
    fluxes_k_means_df = @orderby(vcat(cluster_dfs...), :additive, :final_time)
    return fluxes_k_means_df
end

end

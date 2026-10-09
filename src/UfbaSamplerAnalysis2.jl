module UfbaSamplerAnalysis2

using Base.Iterators
using Random
using CSV
using LinearAlgebra
using DataFrames
using DataFramesMeta
using CategoricalArrays
using Chain
using StatsBase
using StatsModels
using Clustering
using MultivariateStats
using Statistics
using KernelDensity
using CairoMakie
using AlgebraOfGraphics
using AlgebraOfGraphics: verbatim
using ColorSchemes
using ProgressMeter

export plot_all_distributions_for_reactions,
    densities_for_reaction,
    histograms_for_reaction_v2,
    load_and_select_sampling_results,
    pca_treatment_effects,
    prepare_treatment_effects_dfs,
    k_means_treatment_effects,
    prepare_median_fluxes_dfs,
    k_means_median_fluxes,
    pca_median_fluxes,
    treatment_distances_from_control,
    plot_treatment_effects_kmeans_pca,
    plot_median_fluxes_kmeans_pca

"""
    load_and_select_sampling_results()

Loads the uFBA sampling results previously generated. Returns the raw sampling results alongside a set of sampling results filtered down to only those additives that have complete models for all possible timepoints.

# Returns
`NamedTuple`

Returns a named tuple with the following fields:
1. `sampling_df`: The wide sampling DataFrame with all additives (including incomplete ones)
2. `working_models_df`: A sampling DataFrame with results from only additives that have working models from all time points.
"""
function load_and_select_sampling_results()
    sampling_filename = joinpath("output", "ufba_sampling.csv")
    sampling_df = CSV.read(sampling_filename, DataFrame)
    final_times = sort(unique(sampling_df.final_time))
    n_final_times = length(final_times)
    complete_additives_df = @chain sampling_df begin
        @groupby(:additive)
        @combine(:n_unique_final_times = length(unique(:final_time)))
        @rsubset(:n_unique_final_times == n_final_times)
        @select(:additive)
    end
    working_models_df = @chain sampling_df begin
        innerjoin(complete_additives_df; on = :additive)
        @orderby(:additive, :final_time)
    end
    result = (sampling_df = sampling_df, working_models_df = working_models_df)
    return result
end

"""
    pivot_sampling_df_long(sampling_df)

Pivots the sampling DataFrame longer. This function duplicates the function of the same name in `UfbaSamplerAnalysis2`. It is duplicated so that neither this module nor the visualization module need to import each other, which would mess up the documentation.

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
    plot_histograms = false,
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
            category = rxn_ids_to_strings[reaction_id]["category"]
            # reaction_name = rxn_ids_to_strings[reaction_id]["name"]
            if plot_histograms
                fig_hist = histograms_for_reaction_v2(
                    long_sampling_df,
                    reaction_id,
                    reaction_string,
                    subsystem;
                    control_additive = control_additive,
                    treatment_additive = treatment_additive,
                    bins = bins,
                )
                save(filename_hist, fig_hist)
            end
            fig_density = densities_for_reaction(
                long_sampling_df,
                reaction_id,
                reaction_string,
                category,
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
            @rtransform(:time_span = "Week $(:final_time)")
        end
        return plt_df
    end
end

"""
    densities_for_reaction(long_sampling_df, reaction_id, reaction_string, category, subsystem)

Plots KDEs of the flux distributions for the reaction in the various additives.

# Arguments
1. `long_sampling_df`: Sampling DataFrame, pivoted long
2. `reaction_id`: The reaction id for which the samples are being plotted.
3. `reaction_string`: The human-readable reaction string to place as a subtitle on the plot.
4. `category`: Human-readable category of the reaction.
5. `subsystem`: Human-readable subsystem of the reaction.
6. `control_additive = nothing`: If both this and `treatment_additive` are specified, the densities/histograms will just be between these two conditions and the color palette will be fixed.
7. `treatment_additive = nothing`: See `control_additive` above.

# Returns
`Figure`

Returns a Makie `Figure` to display or save.
"""
function densities_for_reaction(
    long_sampling_df,
    reaction_id,
    reaction_string,
    category,
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
    clean_reaction_string = replace(reaction_string, r"\s*\([^)]*\)" => "")
    title = "$clean_reaction_id ($category)\n$clean_reaction_string"
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
        axis = (;
            xlabelsize = 18,
            ylabelsize = 18,
            xlabelfont = :bold,
            ylabelfont = :bold,
            xticklabelsize = 15,
            yticklabelsize = 15,
        ),
        legend = (; titlesize = 17, labelsize = 15),
        facet = (; linkxaxes = :all, linkyaxes = :all),
        figure = (; title = title, titlesize = 22, size = (700, 700)),
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

"""
    prepare_treatment_effects_dfs(
        control_vs_treatments_signif_df;
        complete_only = true,
    )

Prepares the DataFrames for PCA and k-means analysis of treatment vs control Cohen's effects.
    
# Arguments
1. `control_vs_treatments_signif_df`: Treatment effects filtered down to significant treatments.
2. `complete_only = true`: If `true`, only gathers results from treatments where all possible timepoints are present.

# Returns
`NamedTuple`

Returns a named tuple with the following elements:
1. `final_times`: Sorted vector of all final times found.
2. `effects_dfs`: Dictionary mapping integer final times to DataFrames, pivoted wide, with columns of Cohen's effects for each reaction.
"""
function prepare_treatment_effects_dfs(
    control_vs_treatments_signif_df;
    complete_only = true,
)
    final_times = sort(unique(control_vs_treatments_signif_df.final_time))
    n_final_times = length(final_times)
    complete_df = @chain control_vs_treatments_signif_df begin
        @groupby(:treatment_additive)
        @combine(:n_unique_final_times = length(unique(:final_time)))
        @rsubset(:n_unique_final_times == n_final_times)
        @select(:treatment_additive)
    end
    complete_vs_treatments_df =
        complete_only ?
        innerjoin(control_vs_treatments_signif_df, complete_df; on = :treatment_additive) :
        control_vs_treatments_signif_df
    effects_dfs = DataFrame[]
    for final_time in final_times
        df = @chain complete_vs_treatments_df begin
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

"""
    k_means_treatment_effects(
        prepared_treatments_result;
        k = 5,
        seed = 123,
        maxiter = 300,
    )

Finds k-means clusters of treatments using treatment effects as features calculated for each time point.

# Arguments
1. `prepared_treatments_result`: Data prepared by [`prepare_treatment_effects_dfs`](@ref BloodStorageInSilico.UfbaSamplerAnalysis2.prepare_treatment_effects_dfs)
2. `k = 5`: The number of clusters to create.
3. `seed = 123`: RNG seed
4. `maxiter = 300`: Maximum iterations of clustering algorithm.

# Returns
`DataFrame`

K-Means cluster assignments for each additive and final time, ordered by additive and final time.
"""
function k_means_treatment_effects(
    prepared_treatments_result;
    k = 5,
    seed = 123,
    maxiter = 300,
)
    final_times = prepared_treatments_result.final_times
    effects_dfs = prepared_treatments_result.effects_dfs
    cluster_dfs = DataFrame[]
    for (final_time, effects_df) in zip(final_times, effects_dfs)
        feature_cols = names(effects_df, Not([:treatment_additive, :final_time]))
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

"""
    pca_treatment_effects(prepared_treatments_result; n_pcs = 5)

Finds principal components of treatments treatment effects as features calculated for each time point.

# Arguments
1. `prepared_treatments_result`: Data prepared by [`prepare_treatment_effects_dfs`](@ref BloodStorageInSilico.UfbaSamplerAnalysis2.prepare_treatment_effects_dfs)
2. `n_pcs = 5`: Number of principal components to calculate.

# Returns
`NamedTuple`

Returns a named tuple with the following elements:
1. `pca_df`: Rows with treatment, final time, and PC values
2. `loadings_df`: Loadings of each reaction effect on each PC.
"""
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
        size_transpose = size(collect(scores'))
        if size_transpose[2] != length(pc_names)
            error(
                """
              PCA output dimension mismatch.

              Expected the number of columns in `size_transpose` to match the number of PC names.

              Observed:
              size(size_transpose, 2) = $(size_transpose[2])
              length(pc_names)        = $(length(pc_names))

              Suggested resolution:
              This is often caused from attempting to run a second phase analysis on an incomplete
              set of models (such as if --nmodels was not set to -1)
          """,
            )
        end
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

#####################################################################
# MEDIAN FLUX PCA/K-MEANS                                           #
#####################################################################

"""
    zscore_col(xs)

Calculates the z-scores of a column of values. Helper function for [`prepare_median_fluxes_dfs`](@ref BloodStorageInSilico.UfbaSamplerAnalysis2.prepare_median_fluxes_dfs)

# Arguments
1. `xs`: Vector of values for which to compute z-scores.

# Returns
`Vector{Float64}`

Returns vector of z-score values.
"""
function zscore_col(xs)
    μ = mean(skipmissing(xs))
    σ = std(skipmissing(xs))
    if isapprox(σ, 0.0)
        return fill(0.0, length(xs))
    else
        return (xs .- μ) ./ σ
    end
end

"""
    prepare_median_fluxes_dfs(median_fluxes_df; complete_only = true)

Prepares centered and scaled median fluxes for each treatment grouped by time.

# Arguments
1. `median_fluxes_df`: DataFrame of median flux values for each additive and time point.
2. `complete_only = true`: If `true`, filters treatments to only those that have results for the complete set of time points.

# Returns
`NamedTuple`

Returns a named tuple with the following elements:
1. `final_times`: The vector of all final times.
2. `centered_scaled_dfs`: Dictionary mapping final times to centered and scaled DataFrames of median fluxes per additive.
"""
function prepare_median_fluxes_dfs(median_fluxes_df; complete_only = true)
    final_times = sort(unique(median_fluxes_df.final_time))
    n_final_times = length(final_times)
    complete_df = @chain median_fluxes_df begin
        @groupby(:additive, :reaction_id)
        @combine(:n_unique_final_times = length(unique(:final_time)))
        @rsubset(:n_unique_final_times == n_final_times)
        @combine(:additive = unique(:additive))
        @select(:additive)
    end
    median_fluxes_df_2 =
        complete_only ? innerjoin(median_fluxes_df, complete_df; on = :additive) :
        median_fluxes_df
    centered_scaled_dfs = DataFrame[]
    for final_time in final_times
        df = @chain median_fluxes_df_2 begin
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

"""
   treatment_distances_from_control(
        prepare_median_fluxes_result;
        control_additive = "AS3",
    )

For each time point and using centered and scaled median fluxes for all reactions as features, computes Euclidean distance from the control for each additive.

# Argument
1. `prepare_median_fluxes_result`: Result from [`prepare_median_fluxes_dfs`](@ref BloodStorageInSilico.UfbaSamplerAnalysis2.prepare_median_fluxes_dfs)
2. `control_additive = "AS3"`: Additive to be treated as the control.

# Returns
`DataFrame`

DataFrame with additives ranked by their distances from control.
"""
function treatment_distances_from_control(
    prepare_median_fluxes_result;
    control_additive = "AS3",
)
    centered_scaled_dfs = prepare_median_fluxes_result.centered_scaled_dfs
    feature_cols = names(centered_scaled_dfs[1], Not([:additive, :final_time]))
    results = DataFrame[]
    for sdf in centered_scaled_dfs
        time_value = sdf[1, :final_time]
        control_df = @rsubset(sdf, :additive == control_additive)
        if nrow(control_df) != 1
            error(
                "Expected exactly one AS3 control row for $(time_col) = $(time_value); " *
                "found $(nrow(control_df)).",
            )
        end
        X = Matrix(select(sdf, feature_cols))
        x_control = vec(Matrix(select(control_df, feature_cols)))
        distances = norm.(eachrow(X .- x_control'))
        out = copy(sdf)
        out[!, :distance_127D] = distances
        out[!, :is_control] = out[!, :additive] .== control_additive
        push!(results, out)
    end
    result_df = vcat(results...)
    ranked_df = @chain result_df begin
        @orderby(:final_time, :distance_127D)
        @groupby(:final_time)
        @transform(:rank_127D_most_control_like = 1:length(:distance_127D))
        @orderby(:final_time, :rank_127D_most_control_like)
        @select(
            :additive,
            :is_control,
            :final_time,
            :distance_127D,
            :rank_127D_most_control_like
        )
    end
    return ranked_df
end

"""
    k_means_median_fluxes(prepared_medians_result; k = 5, seed = 123, maxiter = 300)

Calculates k-means clusters for treatments using median fluxes of reactions as features.

# Arguments
1. `prepare_median_fluxes_result`: Result from [`prepare_median_fluxes_dfs`](@ref BloodStorageInSilico.UfbaSamplerAnalysis2.prepare_median_fluxes_dfs)
2. `k = 5`: The number of clusters to create.
3. `seed = 123`: RNG seed
4. `maxiter = 300`: Maximum iterations of clustering algorithm.

# Returns
`DataFrame`

K-Means cluster assignments for each additive and final time, ordered by additive and final time.
"""
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

"""
    pca_median_fluxes(prepared_medians_result; n_pcs = 5)

Finds principal components of median fluxes for each additive/time point

# Arguments
1. `prepare_median_fluxes_result`: Result from [`prepare_median_fluxes_dfs`](@ref BloodStorageInSilico.UfbaSamplerAnalysis2.prepare_median_fluxes_dfs)
2. `n_pcs = 5`: Number of principal components to calculate.

# Returns
`NamedTuple`

Returns a named tuple with the following elements:
1. `pca_df`: Rows with treatment, final time, and PC values
2. `loadings_df`: Loadings of each reaction effect on each PC.
"""
function pca_median_fluxes(prepared_medians_result; n_pcs = 5)
    final_times = prepared_medians_result.final_times
    centered_scaled_dfs = prepared_medians_result.centered_scaled_dfs
    pc_names = Symbol.("PC", 1:n_pcs)
    pca_dfs = DataFrame[]
    loadings_dfs = DataFrame[]
    for (final_time, df) in zip(final_times, centered_scaled_dfs)
        reaction_ids = names(select(df, Not([:additive, :final_time])))
        X = Matrix(select(df, Not([:additive, :final_time])))
        Xt = copy(X')
        M = fit(PCA, Xt; maxoutdim = n_pcs, mean = false)
        scores = MultivariateStats.transform(M, Xt)
        pca_df = DataFrame(collect(scores'), pc_names)
        insertcols!(
            pca_df,
            1,
            :additive => df.additive,
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
    pca_df = @orderby(vcat(pca_dfs...), :additive, :final_time)
    loadings_df = @orderby(vcat(loadings_dfs...), :reaction_id, :final_time)
    result = (pca_df = pca_df, loadings_df = loadings_df)
    return result
end

#####################################################################
# PLOT PCA                                                          #
#####################################################################

"""
    plot_treatment_effects_kmeans_pca(
        pca_df,
        treatment_k_means_df;
        color_clusters = false,
    )

Plots the PCA of the standardized treatment Cohen's effects with PC2 on vertical axis and PC1 on horizontal axis with a dot for each treatment additive and each dot labeled with the name of the treatment additive. Saves one plot per treatment additive to `output/viz_effects_kmeans_pca` and displays a progress bar as it goes.

# Arguments
1. `pca_df`: PCA DataFrame from [`pca_median_fluxes`](@ref BloodStorageInSilico.UfbaSamplerAnalysis2.pca_median_fluxes)
2. `treatment_k_means_df`: K-Means DataFrame from [`k_means_median_fluxes`](@ref BloodStorageInSilico.UfbaSamplerAnalysis2.k_means_median_fluxes)
3. `color_clusters = false`: If `true`, colors the dots according to the cluster number in the k-means. **Caveat: Just raw k-means clusters will be different per time point.**
"""
function plot_treatment_effects_kmeans_pca(
    pca_df,
    treatment_k_means_df;
    color_clusters = false,
)
    n_clusters = maximum(treatment_k_means_df.cluster)
    all_time_df = @chain pca_df begin
        innerjoin(treatment_k_means_df; on = [:treatment_additive, :final_time])
        @transform(:cluster = categorical(:cluster))
        @orderby(:final_time, :cluster)
        @select(:final_time, :cluster, :PC1, :PC2, :treatment_additive)
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
        title = "Effects K-Means PCA Final Time $final_time"
        plt_df = @rsubset(all_time_df, :final_time == final_time)
        x_min, x_max = extrema(plt_df.PC1)
        x_span = x_max - x_min
        x_span_for_padding = iszero(x_span) ? 1.0 : x_span
        left_pad = 0.05 * x_span_for_padding
        right_pad = 0.30 * x_span_for_padding
        x_limits = (x_min - left_pad, x_max + right_pad)
        y_min, y_max = extrema(plt_df.PC2)
        y_span = y_max - y_min
        y_span_for_padding = iszero(y_span) ? 1.0 : y_span
        y_pad = 0.12 * y_span_for_padding
        y_limits = (y_min - y_pad, y_max + y_pad)
        scatter_plt =
            color_clusters ?
            data(plt_df) * (
                mapping(:PC1, :PC2, color = :cluster) *
                visual(Scatter, markersize = 14, alpha = 0.75) +
                mapping(:PC1, :PC2, text = :treatment_additive => verbatim) *
                visual(Makie.Text, align = (:left, :bottom), offset = (5, 5))
            ) :
            data(plt_df) * (
                mapping(:PC1, :PC2) * visual(Scatter, markersize = 14, alpha = 0.75) +
                mapping(:PC1, :PC2, text = :treatment_additive => verbatim) *
                visual(Makie.Text, align = (:left, :bottom), offset = (5, 5))
            )
        fig =
            color_clusters ?
            draw(
                scatter_plt,
                scales(Color = (; palette = cluster_colors)),
                figure = (; size = (500, 500)),
                axis = (; title = title, limits = (x_limits, y_limits)),
            ) :
            draw(
                scatter_plt,
                figure = (; size = (500, 500)),
                axis = (; title = title, limits = (x_limits, y_limits)),
            )
        save(filename, fig)
        next!(prog)
    end
end

"""
    plot_treatment_effects_kmeans_pca(
        pca_df,
        treatment_k_means_df;
        color_clusters = false,
    )

Plots the PCA of the standardized median fluxes with PC2 on vertical axis and PC1 on horizontal axis with a dot for each treatment additive and each dot labeled with the name of the treatment additive. Saves one plot per treatment additive to `output/viz_fluxes_kmeans_pca` and displays a progress bar as it goes.

# Arguments
1. `pca_df`: PCA DataFrame from [`pca_median_fluxes`](@ref BloodStorageInSilico.UfbaSamplerAnalysis2.pca_median_fluxes)
2. `treatment_k_means_df`: K-Means DataFrame from [`k_means_median_fluxes`](@ref BloodStorageInSilico.UfbaSamplerAnalysis2.k_means_median_fluxes)
3. `color_clusters = false`: If `true`, colors the dots according to the cluster number in the k-means. **Caveat: Just raw k-means clusters will be different per time point.**
"""
function plot_median_fluxes_kmeans_pca(pca_df, fluxes_k_means_df; color_clusters = false)
    n_clusters = maximum(fluxes_k_means_df.cluster)
    all_time_df = @chain pca_df begin
        innerjoin(fluxes_k_means_df; on = [:additive, :final_time])
        @transform(:cluster = categorical(:cluster))
        @orderby(:final_time, :cluster)
        @select(:final_time, :cluster, :PC1, :PC2, :additive)
    end
    cluster_colors = get(colorschemes[:okabe_ito], range(0, 1, length = n_clusters))
    final_times = sort(unique(all_time_df.final_time))
    n_plots = length(final_times)
    prog = Progress(n_plots, "Writing median flux k-means PCA plots")
    for final_time in final_times
        filename = joinpath(
            "output",
            "viz_fluxes_kmeans_pca",
            "median_flux_kmeans_pca_$(final_time).png",
        )
        title = "Median Flux PCA, Week $final_time"
        plt_df = @rsubset(all_time_df, :final_time == final_time)
        x_min, x_max = extrema(plt_df.PC1)
        x_span = x_max - x_min
        x_span_for_padding = iszero(x_span) ? 1.0 : x_span
        left_pad = 0.05 * x_span_for_padding
        right_pad = 0.40 * x_span_for_padding
        x_limits = (x_min - left_pad, x_max + right_pad)
        y_min, y_max = extrema(plt_df.PC2)
        y_span = y_max - y_min
        y_span_for_padding = iszero(y_span) ? 1.0 : y_span
        y_pad = 0.12 * y_span_for_padding
        y_limits = (y_min - y_pad, y_max + y_pad)
        scatter_plt =
            color_clusters ?
            data(plt_df) * (
                mapping(:PC1, :PC2, color = :cluster) *
                visual(Scatter, markersize = 14, alpha = 0.50) +
                mapping(:PC1, :PC2, text = :additive => verbatim) * visual(
                    Makie.Text,
                    align = (:left, :bottom),
                    offset = (5, 5),
                    fontsize = 22,
                )
            ) :
            data(plt_df) * (
                mapping(:PC1, :PC2) * visual(Scatter, markersize = 14, alpha = 0.50) +
                mapping(:PC1, :PC2, text = :additive => verbatim) * visual(
                    Makie.Text,
                    align = (:left, :bottom),
                    offset = (5, 5),
                    fontsize = 22,
                )
            )
        fig =
            color_clusters ?
            draw(
                scatter_plt,
                scales(Color = (; palette = cluster_colors)),
                figure = (; size = (700, 700)),
                axis = (;
                    title = title,
                    titlesize = 28,
                    xlabelsize = 22,
                    ylabelsize = 22,
                    xlabelfont = :bold,
                    ylabelfont = :bold,
                    xticklabelsize = 15,
                    yticklabelsize = 15,
                    limits = (x_limits, y_limits),
                ),
            ) :
            draw(
                scatter_plt,
                figure = (; size = (700, 700)),
                axis = (;
                    title = title,
                    titlesize = 28,
                    xlabelsize = 22,
                    ylabelsize = 22,
                    xlabelfont = :bold,
                    ylabelfont = :bold,
                    xticklabelsize = 15,
                    yticklabelsize = 15,
                    limits = (x_limits, y_limits),
                ),
            )
        save(filename, fig)
        next!(prog)
    end
end

end

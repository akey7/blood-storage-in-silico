module RawRelativeIntensities

using Base.Iterators
using CSV
using DataFrames
using DataFramesMeta
using MultivariateStats
using StatsBase
using StatsModels
using Statistics
using Makie
using GLMakie
using ThreadsX
using ProgressMeter

export load_relative_intensities,
    pca_relative_intensities,
    plot_pca_panels,
    display_pca_scores_3d,
    gather_pca_scores,
    calc_pca_scores_3d_limits,
    pca_loadings_report,
    plot_single_additive_2d_pcas

"""
    load_relative_intensities()

Loads the relative quantification (intensity) and pivots it long.

# Returns
`DataFrame`

Returns a long DataFrame with the following columns: 

1. `:Sample`, the sample id

2. `:Time` the time point of the measurement (in weeks)

3. `:Additive`: Additive the measurement was taken in.

4. `:MixedName`: The name of either a single compound or group of compounds under the same peak.

5. `:Intensity`: The integrated area of the peak.
"""
function load_relative_intensities()
    relative_filename = joinpath("input", "Data Sheet 1.CSV")
    wide_df = CSV.read(relative_filename, DataFrame)
    long_df = stack(
        wide_df,
        Not([:Sample, :Time, :Additive]),
        variable_name = :MixedName,
        value_name = :Intensity,
    )
    return long_df
end

"""
    pca_relative_intensities(long_df, additive)

Perform a robust PCA of the relative intensity data of metabolites within a given additive. Handles NaNs and missing values gracefully. Centers and scales prior to PCA.

# Arguments
1. `long_df`: The long dataframe as returned by [`load_relative_intensities`](@ref BloodStorageInSilico.RawRelativeIntensities.load_relative_intensities)
2. `additive`: The additive for which to perform the PCA

# Returns
`NamedTuple`

1. `model`: PCA model produced, which enables accessing properties of the PCA model downstream.
2. `scores`: Scores of each observation so that principal component scatter plots can be made.
3. `patient_time_labels`: Labels for each observation of patient and time.
4. `kept_columns`: List of columns that were kept for the PCA after data cleaning
5. `wide_df`: Wide DataFrame used to make the `Matrix` for the PCA.
6. `additive`: The additive the PCA was performed for.
"""
function pca_relative_intensities(long_df, additive)
    # @info "Beginning PCA for $additive"
    wide_df = @chain long_df begin
        @rsubset(:Additive == additive)
        @rtransform(:Patient = split(:Sample, "_")[3][1:2])
        @select(:MixedName, :Patient, :Time, :Intensity)
        unstack([:Patient, :Time], :MixedName, :Intensity, combine = first)
        @orderby(:Patient, :Time)
    end
    metabolite_names = names(wide_df)[3:end]
    # display(metabolite_names) 
    patient_time_labels = @select(wide_df, :Patient, :Time)
    # display(patient_time_labels)
    X = Matrix(select(wide_df, Not([:Patient, :Time])))
    colmeans = map(zip(metabolite_names, eachcol(X))) do p
        metabolite_name, c = p
        m = mean(skipmissing(c))
        if isfinite(m)
            return m
        else
            println("$metabolite_name has non-finite mean")
            return missing
        end
    end
    for (metabolite_name, j) in zip(metabolite_names, axes(X, 2))
        if ismissing(colmeans[j])
            continue
        end
        for (patient_time_label, i) in zip(eachrow(patient_time_labels), axes(X, 1))
            if ismissing(X[i, j])
                println("$patient_time_label, $metabolite_name is missing")
                X[i, j] = colmeans[j]
            end
        end
    end
    Xf = Array{Float64}(undef, size(X))
    for j in axes(X, 2), i in axes(X, 1)
        Xf[i, j] = ismissing(X[i, j]) ? NaN : Float64(X[i, j])
    end
    good_cols = trues(size(Xf, 2))
    for j in axes(Xf, 2)
        col = view(Xf, :, j)
        if any(!isfinite, col)
            good_cols[j] = false
            continue
        end
        s = std(col)
        if !isfinite(s) || s == 0.0
            good_cols[j] = false
        end
    end
    # display(good_cols)
    Xf = Xf[:, good_cols]
    if size(Xf, 2) == 0
        error("After filtering, no valid metabolite columns remain for PCA.")
    end
    zt = StatsBase.fit(StatsBase.ZScoreTransform, Xf; dims = 1)
    Xz = StatsBase.transform(zt, Xf)
    # Check for NaN and missing
    for j in axes(X, 2), i in axes(X, 1)
        if isnan(Xf[i, j]) || ismissing(Xf[i, j])
            println("Xf[$i, $j] is NaN or missing")
        end
    end
    Xzt = copy(Xz')
    M = fit(PCA, Xzt; maxoutdim = 6, mean = false)
    # display(M)
    scores = MultivariateStats.transform(M, Xzt)
    kept_columns = findall(good_cols)
    # display(kept_columns)
    # display(M)
    # @info "Finished PCA for $additive"
    return (
        additive = additive,
        model = M,
        scores = scores,
        patient_time_labels = patient_time_labels,
        kept_columns = kept_columns,
        wide_df = wide_df,
    )
end

"""
    plot_all_2d_pcas(long_df)

Plot the 2D PCA multi panel plots.

# Arguments
1. `long_df`: Long DataFrame of relative intensities.
"""
function plot_single_additive_2d_pcas(long_df)
    additives = sort(unique(long_df.Additive))
    n_additives = length(additives)
    prog = Progress(n_additives, desc = "Plotting 2D PCAs")
    for additive in additives
        pca_result = pca_relative_intensities(long_df, additive)
        fig = plot_pca_panels(pca_result, additive)
        filename = joinpath("output", "pca_plots", "PCA $additive.png")
        save(filename, fig)
        next!(prog)
    end
end

"""
    plot_pca_panels(pca_result, super_title)

Using [`plot_pca_scores`](@ref BloodStorageInSilico.RawRelativeIntensities.plot_pca_scores) and [`plot_pca_scree`](@ref BloodStorageInSilico.RawRelativeIntensities.plot_pca_scree), assemble a 2D set of panels for to plot the PCA results.

# Arguments
1. `pca_result`: Result from [`pca_relative_intensities`](@ref BloodStorageInSilico.RawRelativeIntensities.pca_relative_intensities)
2. `super_title`: The super title to put over the top of both panels.

# Returns
`Figure`

Returns a Makie `Figure` object to be shown or saved.
"""
function plot_pca_panels(pca_result, super_title)
    fig = Figure(; size = (1280, 720))
    plot_pca_scores(pca_result, fig)
    plot_pca_scree(pca_result, fig)
    Label(fig[0, :], text = super_title, fontsize = 50)
    return fig
end

"""
    plot_pca_scores(pca_result, fig)

Plot a panel of the first two PCs against each other in a scatter plot.

# Arguments
1. `pca_result`: Result from [`pca_relative_intensities`](@ref BloodStorageInSilico.RawRelativeIntensities.pca_relative_intensities)
2. `fig`: A Makie figure to plot onto.
"""
function plot_pca_scores(pca_result, fig)
    M = pca_result.model
    scores = pca_result.scores
    pc1 = scores[1, :]
    pc2 = scores[2, :]
    time_labels = pca_result.patient_time_labels.Time
    time_color_map = Dict(
        1 => "#006CD1",
        2 => "#E66100",
        3 => "#5D3A9B",
        4 => "#40B0A6",
        5 => "#AFAF01",
        6 => "#222222",
    )
    time_shape_map = Dict(
        1 => :circle,
        2 => :rect,
        3 => :diamond,
        4 => :cross,
        5 => :utriangle,
        6 => :dtriangle,
    )
    var_explained = principalvars(M) ./ tvar(M)
    xlabel = "PC1 $(round(var_explained[1]*100, digits = 2))%"
    ylabel = "PC2 $(round(var_explained[2]*100, digits = 2))%"
    title = "PCA of Timeseries"
    ax_scatter = Axis(fig[1:3, 2:3], xlabel = xlabel, ylabel = ylabel, title = title)
    ax_hist = Axis(fig[4, 2:3])
    # for (x, y, tl) in zip(pc1, pc2, time_labels)
    #     text!(ax, x, y; text = string(tl), offset = (5, -5), align = (:left, :bottom))
    # end
    unique_times = sort(unique(time_labels))
    for t in unique_times
        idxs = findall(==(t), time_labels)
        scatter!(
            ax_scatter,
            pc1[idxs],
            pc2[idxs],
            color = time_color_map[t],
            marker = time_shape_map[t],
            markersize = 20,
            label = string(t),
            alpha = 0.75,
        )
    end
    hist!(ax_hist, pc1; bins = 6)
    axislegend(ax_scatter; position = :rb)
end

"""
    plot_pca_scree(pca_result, fig)

Plots a PCA scree plot panel onto the given figure.

# Arguments
1. `pca_result`: Result from [`pca_relative_intensities`](@ref BloodStorageInSilico.RawRelativeIntensities.pca_relative_intensities)
2. `fig`: Make `Figure` to plot the panel onto.
"""
function plot_pca_scree(pca_result, fig)
    M = pca_result.model
    var_explained = principalvars(M) ./ tvar(M)
    ys = cumsum(var_explained) .* 100
    xs = eachindex(ys)
    yticks = range(0.0, 100.0, 5)
    ytick_labels = string.(round.(yticks))
    xlabel = "Component"
    ylabel = "Percent"
    title = "Cumulative variance explained"
    ax = Axis(
        fig[2:3, 1],
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = (xs, string.(xs)),
        yticks = (yticks, ytick_labels),
        limits = (nothing, nothing, 0.0, 100.0),
    )
    lines!(ax, xs, ys)
    scatter!(ax, xs[2], ys[2], markersize = 20, color = :crimson)
    text!(
        ax,
        xs[2],
        ys[2];
        text = "$(round(ys[2], digits = 2))%",
        offset = (10, -10),
        align = (:left, :bottom),
    )
end

"""
    gather_pca_scores(pca_result)

Using results of the PCA [`pca_relative_intensities`](@ref BloodStorageInSilico.RawRelativeIntensities.pca_relative_intensities), make a DataFrame that can be saved for manual inspection.

# Arguments
1. `pca_result`: The PCA result to gather.

# Returns
`DataFrame`

Returns a DataFrame with the following columns:

1. `patient`: Patient
2. `time`: Time point of observation.
3. `pc1`, `pc2`, `pc3`: Principal components
"""
function gather_pca_scores(pca_result)
    scores = pca_result.scores
    time_labels = pca_result.patient_time_labels.Time
    patient_labels = pca_result.patient_time_labels.Patient
    pc1 = scores[1, :]
    pc2 = scores[2, :]
    pc3 = scores[3, :]
    df = DataFrame(
        patient = patient_labels,
        time = time_labels,
        pc1 = pc1,
        pc2 = pc2,
        pc3 = pc3,
    )
    return df
end

"""
    calc_pca_scores_3d_limits(long_df; margin = 1.1)

Computes 3D axis limits for PCA plots across all additives to set the axis limits of all 3D PCA plots so that plots of different additives can be directly compared.

# Arguments
1. `long_df`: Long dataframe relative intensities.
2. `margin = 1.1`: Multiplier to create margins around the PCA plots.

# Returns
`Tuple{Tuple{Float64,Float64},Tuple{Float64,Float64},Tuple{Float64,Float64}}`

Returns tuple of tuples suitable for passing to GLMakie that define axis limits for each principal component.
"""
function calc_pca_scores_3d_limits(long_df; margin = 1.1)
    additives = sort(unique(long_df.Additive))
    all_pca_results = ThreadsX.map(additives) do additive
        pca_relative_intensities(long_df, additive)
    end
    pc1_min = Inf
    pc2_min = Inf
    pc3_min = Inf
    pc1_max = -Inf
    pc2_max = -Inf
    pc3_max = -Inf
    for pca_result in all_pca_results
        scores = pca_result.scores
        pc1 = scores[1, :]
        pc2 = scores[2, :]
        pc3 = scores[3, :]
        pc1_min = pc1_min > minimum(pc1) ? minimum(pc1) : pc1_min
        pc2_min = pc2_min > minimum(pc2) ? minimum(pc2) : pc2_min
        pc3_min = pc3_min > minimum(pc3) ? minimum(pc3) : pc3_min
        pc1_max = pc1_max < maximum(pc1) ? maximum(pc1) : pc1_max
        pc2_max = pc2_max < maximum(pc2) ? maximum(pc2) : pc2_max
        pc3_max = pc3_max < maximum(pc3) ? maximum(pc3) : pc3_max
    end
    return (
        (pc1_min*margin, pc1_max*margin),
        (pc2_min*margin, pc2_max*margin),
        (pc3_min*margin, pc3_max*margin),
    )
end

"""
    display_pca_scores_3d(pca_result, additive)

Plot **and display** a 3D PCA scatter with GLMakie.

# Arguments
1. `pca_result`: Result from [`pca_relative_intensities`](@ref BloodStorageInSilico.RawRelativeIntensities.pca_relative_intensities).
2. `additive`: Additive to plot the PCA results for.
"""
function display_pca_scores_3d(limits, pca_result, additive)
    @info "Display PCA for $additive"
    M = pca_result.model
    scores = pca_result.scores
    pc1 = scores[1, :]
    pc2 = scores[2, :]
    pc3 = scores[3, :]
    time_labels = pca_result.patient_time_labels.Time
    time_color_map = Dict(
        1 => "#006CD1",
        2 => "#E66100",
        3 => "#5D3A9B",
        4 => "#40B0A6",
        5 => "#AFAF01",
        6 => "#222222",
    )
    time_shape_map = Dict(
        1 => :circle,
        2 => :rect,
        3 => :diamond,
        4 => :cross,
        5 => :utriangle,
        6 => :dtriangle,
    )
    var_explained = principalvars(M) ./ tvar(M)
    xlabel = "PC1 $(round(var_explained[1]*100, digits = 2))%"
    ylabel = "PC2 $(round(var_explained[2]*100, digits = 2))%"
    zlabel = "PC3 $(round(var_explained[3]*100, digits = 2))%"
    title = "$additive PCA"
    fig = Figure(size = (750, 750), figure_padding = 75)
    ax_scatter_3d = Axis3(
        fig[1, 1],
        xlabel = xlabel,
        ylabel = ylabel,
        zlabel = zlabel,
        title = title,
        limits = limits,
    )
    unique_times = sort(unique(time_labels))
    for t in unique_times
        idxs = findall(==(t), time_labels)
        n_points = length(idxs)
        println("$n_points at time $t")
        scatter!(
            ax_scatter_3d,
            pc1[idxs],
            pc2[idxs],
            pc3[idxs],
            color = time_color_map[t],
            marker = time_shape_map[t],
            markersize = 15,
            alpha = 0.75,
            label = string(t),
        )
    end
    axislegend(ax_scatter_3d, "Week"; position = :rb, margin = (-30, -30, -30, -30))
    @info "Finished preparing PCA plot"
    GLMakie.display(fig)
end

"""
    extract_pca_loadings(pca_result, additive)

Extracts the loadings of the metabolite features on each of the PCs. Used by [`pca_loadings_report`](@ref BloodStorageInSilico.RawRelativeIntensities.pca_loadings_report).

# Arguments
1. `pca_result`: Result from The long DataFrame from [`pca_relative_intensities`](@ref BloodStorageInSilico.RawRelativeIntensities.pca_relative_intensities).

# Returns
`DataFrame`

Returns a DataFrame, ordered by the column `pc1_loading`, that has the following columns:

1. `additive`: Additive the PCA was performed for.
2. `metabolite_name`: Names of the metabolites.
3. `pc1_loading`: Loadings on the first PC.
4. `pc2_loading`: Loadings on the second PC.
5. `pc3_loading`: Loadings on the third PC.
6. `pc4_loading`: Loadings on the fourth PC.
7. `pc5_loading`: Loadings on the fifth PC.
8. `pc6_loading`: Loadings on the sixth PC.
"""
function extract_pca_loadings(pca_result)
    additive = pca_result.additive
    M = pca_result.model
    L = loadings(M)
    pc1_loadings = L[:, 1]
    pc2_loadings = L[:, 2]
    pc3_loadings = L[:, 3]
    pc4_loadings = L[:, 4]
    pc5_loadings = L[:, 5]
    pc6_loadings = L[:, 6]
    kept_columns = pca_result.kept_columns
    wide_df = pca_result.wide_df
    metabolite_names = names(select(wide_df, Not(:Time)))[kept_columns]
    loadings_df = DataFrame(
        additive = additive,
        metabolite_name = metabolite_names,
        pc1_loading = pc1_loadings,
        pc2_loading = pc2_loadings,
        pc3_loading = pc3_loadings,
        pc4_loading = pc4_loadings,
        pc5_loading = pc5_loadings,
        pc6_loading = pc6_loadings,
    )
    result_df = @orderby(loadings_df, :pc1_loading)
    return result_df
end

"""
    loadings_report(long_df)

Collect and concatenate all PC loadings in all additives into a DataFrame.

# Arguments
1. `long_df`: Long DataFrame of relative intensities

# Returns
`DataFrame`

Returns a DataFrame, ordered by the column `pc1_loading`, that has the following columns:

1. `additive`: Additive the PCA was performed for.
2. `metabolite_name`: Names of the metabolites.
3. `pc1_loading`: Loadings on the first PC.
4. `pc2_loading`: Loadings on the second PC.
5. `pc3_loading`: Loadings on the third PC.
6. `pc4_loading`: Loadings on the fourth PC.
7. `pc5_loading`: Loadings on the fifth PC.
8. `pc6_loading`: Loadings on the sixth PC.
"""
function pca_loadings_report(long_df)
    additives = sort(unique(long_df.Additive))
    all_pca_results = ThreadsX.map(additives) do additive
        pca_result = pca_relative_intensities(long_df, additive)
        extract_pca_loadings(pca_result)
    end
    result_df = vcat(all_pca_results...)
    return result_df
end

end

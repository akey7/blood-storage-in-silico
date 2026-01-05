module RawRelativeIntensities

using CSV
using DataFrames
using DataFramesMeta
using MultivariateStats
using StatsBase
using StatsModels
using Statistics
using Makie
using GLMakie

export load_relative_intensities,
    pca_relative_intensities, plot_pca_panels, display_pca_scores_3d

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

function pca_relative_intensities(long_df, additive)
    @info "Beginning PCA"
    wide_df = @chain long_df begin
        @rsubset(:Additive == additive)
        @rtransform(:Patient = :Sample[7:8])
        @select(:MixedName, :Patient, :Time, :Intensity)
        unstack([:Patient, :Time], :MixedName, :Intensity, combine = first)
        @orderby(:Patient, :Time)
    end
    metabolite_names = names(wide_df)[3:end]
    # display(metabolite_names) 
    patient_time_labels = @select(wide_df, :Patient, :Time)
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
                println("$patient_time_label, $metabolite is missing")
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
    @info "Finished PCA"
    return (
        model = M,
        scores = scores,
        patient_time_labels = patient_time_labels,
        kept_columns = kept_columns,
        wide_df = wide_df,
    )
end

function plot_pca_panels(pca_result, super_title)
    fig = Figure(; size = (1280, 720))
    plot_pca_scores(pca_result, fig)
    plot_pca_scree(pca_result, fig)
    # plot_pca_loadings(pca_result, fig)
    Label(fig[0, :], text = super_title, fontsize = 50)
    return fig
end

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

function display_pca_scores_3d(pca_result, additive)
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
    ax_scatter_3d =
        Axis3(fig[1, 1], xlabel = xlabel, ylabel = ylabel, zlabel = zlabel, title = title)
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
            markersize = 20,
            alpha = 0.75,
            label = string(t),
        )
    end
    axislegend(ax_scatter_3d, "Week"; position = :rb, margin = (-30, -30, -30, -30))
    @info "Finished preparing PCA plot"
    GLMakie.display(fig)
end

end

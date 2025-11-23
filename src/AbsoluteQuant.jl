module AbsoluteQuant

using Base.Iterators
using CSV
using XLSX
using DataFrames
using DataFramesMeta
using StatsBase
using AlgebraOfGraphics
using CairoMakie
using Makie
using Clustering
using Distances
using ShiftedArrays
using MultivariateStats
using GLM
using StatsModels
using Statistics
using Random
using ThreadsX

export load_absolute_quant,
    load_relative_quant,
    combine_relative_and_absolute_quant,
    cluster_all_additives_all_n_clusters,
    plot_elbows,
    plot_c_means_all_additives,
    plot_all_mM_timeseries,
    diff_mM,
    pca_timeseries,
    regress_concentration_dxdt,
    plot_pca_loadings

function load_absolute_quant()
    absolute_filename = joinpath("input", "Absolute Quant Data Sheet.xlsx")
    cells_day_1_df = DataFrame(XLSX.readtable(absolute_filename, "cells_day_1"))
    absolute_metabolite_ids = DataFrame(XLSX.readtable(absolute_filename, "metabolite_ids"))
    absolute_quant_df = @chain cells_day_1_df begin
        stack(Not(:id), variable_name = :mixed_name, value_name = :mM)
        @rtransform(:sample_set = split(:id, "_")[2])
        innerjoin(absolute_metabolite_ids, on = :mixed_name => :MixedName)
        @rtransform(:prop_mM = :mM * :Proportion)
        select([:sample_set, :id, :Metabolite, :prop_mM])
    end
    absolute_quant_medians_df = @by absolute_quant_df :Metabolite begin
        :median_prop_mM = median(:prop_mM)
    end
    return absolute_quant_df, absolute_quant_medians_df
end

function load_relative_quant()
    relative_filename = joinpath("input", "Data Sheet 1.CSV")
    wide_df = CSV.read(relative_filename, DataFrame)
    proportination_filename = joinpath("input", "Proportionation Sheet 2.csv")
    proportination_df = CSV.read(proportination_filename, DataFrame)
    long_df = stack(
        wide_df,
        Not([:Sample, :Time, :Additive]),
        variable_name = :MixedName,
        value_name = :Intensity,
    )
    control_intensity_df =
        subset(long_df, :Additive => x -> x .== "01-Ctrl AS3", :Time => x -> x .== 1)
    ctrl_time_1_median_df = @by control_intensity_df :MixedName begin
        :CtrlTime1MedianIntensity = median(skipmissing(:Intensity))
    end
    fold_changes_df = @chain long_df begin
        innerjoin(ctrl_time_1_median_df, on = :MixedName)
        @rtransform(:FoldChange = :Intensity / :CtrlTime1MedianIntensity)
        innerjoin(proportination_df, on = :MixedName)
        @select(:Sample, :Time, :Additive, :Metabolite, :FoldChange)
    end
    return fold_changes_df
end

function combine_relative_and_absolute_quant(fold_changes_df, absolute_quant_medians_df)
    long_df = @chain fold_changes_df begin
        innerjoin(absolute_quant_medians_df, on = :Metabolite)
        @rtransform(:absolute_mM = :FoldChange * :median_prop_mM)
        @orderby(:Additive, :Time, :Metabolite)
    end
    wide_df = @chain long_df begin
        @select(:Sample, :Time, :Additive, :Metabolite, :absolute_mM)
        unstack([:Sample, :Time, :Additive], :Metabolite, :absolute_mM, combine = first)
        @orderby(:Additive, :Time)
    end
    return long_df, wide_df
end

function prepare_long_df_for_clustering(long_df, additive)
    long_df_2 = deepcopy(long_df)
    wide_timeseries_df = @chain long_df_2 begin
        @rsubset(:Additive == additive, !isapprox(:absolute_mM, 0.0))
        @select(:Metabolite, :Time, :absolute_mM)
        @orderby(:Metabolite, :Time)
        unstack(:Metabolite, :Time, :absolute_mM, combine = first)
    end
    return wide_timeseries_df
end

function calc_fuzzy_objective(result, X, μ)
    c = result.centers
    W = result.weights
    total = 0.0
    for i in axes(X, 1)
        for j in axes(c, 2)
            total += W[i, j]^μ * sum((X[i, :] ./ -c[:, j]) .^ 2)
        end
    end
    return total
end

function c_means_metabolite_trajectories(
    wide_timeseries_df;
    additive = "01-Ctrl AS3",
    n_clusters = 5,
    μ = 5.0,
)
    X0 = Matrix{Float64}(disallowmissing(wide_timeseries_df[:, Not(:Metabolite)]))
    X = (X0 .- mean(X0, dims = 1)) ./ std(X0, dims = 1)
    nans = count(isnan, X)
    infs = count(isinf, X)
    println("NaN count: $nans, Inf count: $infs")
    @assert nans == 0 "Remove NaNs before clustering."
    @assert infs == 0 "Remove Infs before clustering."
    d, n = size(X')
    println("Shape d x n = $d x $n  (features x observations)")
    constf = sum([iszero(X'[i, :]) for i = 1:d])
    @assert constf == 0 "Drop constant features to avoid zero distances."
    uniq_cols = length(unique(eachcol(X')))
    println("Unique observations: $uniq_cols / $n")
    @assert n_clusters <= uniq_cols "n_clusters must not exceed number of unique observations."
    result = fuzzy_cmeans(
        X',
        n_clusters,
        μ,
        maxiter = 200,
        display = :iter,
        dist_metric = Cityblock(),
    )
    weights_col_names = string.(axes(result.weights, 2))
    memberships_df = DataFrame(result.weights, weights_col_names)
    memberships_df.Metabolite = wide_timeseries_df.Metabolite
    memberships_df[!, :Additive] .= additive
    memberships_df[!, :NClusters] .= n_clusters
    fuzzy_objective = calc_fuzzy_objective(result, X, μ)
    return memberships_df, fuzzy_objective, X
end

function cluster_all_additives_all_n_clusters(long_df; max_clusters = 10)
    additives = unique(long_df.Additive)
    all_memberships_dfs::Dict{Int64,DataFrame} = Dict()
    fuzzy_objectives::Vector{NamedTuple} = []
    for n_clusters = 2:max_clusters
        memberships_dfs::Vector{DataFrame} = []
        for additive in additives
            println(uppercase(additive), " n_clusters ", n_clusters)
            memberships_df, fuzzy_objective, _ = @chain long_df begin
                prepare_long_df_for_clustering(additive)
                c_means_metabolite_trajectories(
                    additive = additive,
                    n_clusters = n_clusters,
                )
            end
            push!(
                fuzzy_objectives,
                (
                    additive = additive,
                    n_clusters = n_clusters,
                    fuzzy_objective = fuzzy_objective,
                ),
            )
            push!(memberships_dfs, memberships_df)
            println(first(memberships_df, 10))
            println("Fuzzy objective: ", fuzzy_objective)
        end
        all_memberships_dfs[n_clusters] = vcat(memberships_dfs...)
    end
    fuzzy_objectives_df = DataFrame(fuzzy_objectives)
    return all_memberships_dfs, fuzzy_objectives_df
end

function plot_elbows(fuzzy_objectives_df)
    fig_filename = joinpath("output", "relative_absolute_c_means", "elbows.png")
    xticks = unique(fuzzy_objectives_df.n_clusters)
    plt =
        data(fuzzy_objectives_df) *
        mapping(
            :n_clusters => "N Clusters",
            :fuzzy_objective => "Fuzzy Objective",
            row = :additive,
        ) *
        visual(Lines)
    figure_options = (; size = (500, 1000), title = "C-Means Objective Elbow Plots")
    axis_options = (; xticks = xticks)
    facet_options = (; linkxaxes = :minimal, linkyaxes = :minimal)
    fig = draw(plt; figure = figure_options, axis = axis_options, facet = facet_options)
    fig_filename = joinpath("output", "relative_absolute_c_means", "elbows.png")
    save(fig_filename, fig)
    println("Wrote $fig_filename")
end

function plot_c_means_for_additive_and_n_clusters(
    long_df,
    all_memberships_dfs,
    additive,
    n_clusters,
)
    df1 = deepcopy(all_memberships_dfs[n_clusters])
    membership_df = @chain df1 begin
        @rsubset(:Additive == additive, :NClusters == n_clusters)
        @orderby(:Metabolite)
    end
    primary_cluster_df = @select membership_df Not([:Metabolite, :Additive, :NClusters])
    primary_cluster_df.primary_cluster =
        [argmax(row) for row in eachrow(primary_cluster_df)]
    membership_df.primary_cluster = primary_cluster_df.primary_cluster
    println(first(membership_df, 10))
    standardization_df = @chain long_df begin
        @rsubset(:Additive == additive)
        @rtransform(:Patient = :Sample[7:8])
        innerjoin(membership_df, on = [:Additive, :Metabolite])
        @select(:primary_cluster, :Patient, :Time, :Metabolite, :absolute_mM)
        unstack(
            [:primary_cluster, :Patient, :Time],
            :Metabolite,
            :absolute_mM,
            combine = first,
        )
        @orderby(:primary_cluster, :Patient, :Time)
    end
    X1 = Matrix(select(standardization_df, Not([:primary_cluster, :Patient, :Time])))
    colmeans = map(c -> mean(skipmissing(c)), eachcol(X1))
    for j = 1:size(X1, 2)
        @inbounds for i = 1:size(X1, 1)
            if ismissing(X1[i, j])
                X1[i, j] = colmeans[j]
            end
        end
    end
    X2 = Float64.(X1)
    zt = StatsBase.fit(StatsBase.ZScoreTransform, X2, dims = 1)
    X3 = StatsBase.transform(zt, X2)
    standardization_df[:, Not([:primary_cluster, :Patient, :Time])] = X3
    plt_df = stack(
        standardization_df,
        Not([:primary_cluster, :Patient, :Time]),
        variable_name = :Metabolite,
        value_name = :standardized_mM,
    )
    time_points = unique(plt_df.Time)
    plt =
        data(plt_df) *
        mapping(
            :Time,
            :standardized_mM => "standardized mmol/L",
            row = :primary_cluster,
            group = :Metabolite,
        ) *
        visual(Lines) *
        visual(alpha = 0.1)
    figure_options = (; size = (500, 1000), title = additive)
    fig = draw(
        plt;
        figure = figure_options,
        axis = (; xticks = time_points),
        facet = (; linkxaxes = :minimal, linkyaxes = :minimal),
    )
    clean_additive = replace(additive, r"[^A-Za-z0-9]" => "_")
    fig_filename = joinpath(
        "output",
        "relative_absolute_c_means",
        "$clean_additive $n_clusters Clusters.png",
    )
    save(fig_filename, fig)
    println("Wrote $fig_filename")
    return membership_df
end

function plot_c_means_all_additives(long_df, all_memberships_dfs, n_clusters)
    additives = unique(long_df.Additive)
    primary_cluster_dfs = []
    for additive in additives
        primary_cluster_df = plot_c_means_for_additive_and_n_clusters(
            long_df,
            all_memberships_dfs,
            additive,
            n_clusters,
        )
        push!(primary_cluster_dfs, primary_cluster_df)
    end
    primary_cluster_df = vcat(primary_cluster_dfs...)
    return primary_cluster_df
end

function plot_all_mM_timeseries(long_df)
    agg_df = @chain long_df begin
        @groupby(:Additive, :Metabolite, :Time)
        @combine(:median_mM = median(skipmissing(:absolute_mM)))
        @orderby(:Additive, :Metabolite, :Time)
    end
    metabolites = unique(agg_df.Metabolite)
    time_points = unique(agg_df.Time)
    for metabolite in metabolites
        clean_metabolite = replace(metabolite, r"[^A-Za-z0-9_]" => "_")
        filename = joinpath("output", "relative_absolute_plots", "$(clean_metabolite).png")
        line_plt_df = @rsubset(agg_df, :Metabolite == metabolite)
        line_plt =
            data(line_plt_df) *
            mapping(:Time, :median_mM => "Median mM", color = :Additive) *
            visual(Lines, linewidth = 2)
        scatter_plt_df = @rsubset(long_df, :Metabolite == metabolite)
        scatter_plt =
            data(scatter_plt_df) *
            mapping(:Time, :absolute_mM, color = :Additive, marker = :Additive) *
            visual(Scatter, markersize = 14, alpha = 0.5)
        plt = line_plt + scatter_plt
        fig = draw(
            plt;
            figure = (; size = (750, 500)),
            axis = (; title = metabolite, xticks = time_points),
        )
        save(filename, fig)
        println("Wrote $filename")
    end
end

function diff_mM(long_df)
    diffed_df = @chain long_df begin
        @rtransform(:Patient = :Sample[7:8])
        @orderby(:Additive, :Metabolite, :Patient, :Time)
        @groupby(:Additive, :Metabolite, :Patient)
        @transform(:diff_mM = :absolute_mM .- ShiftedArrays.lag(:absolute_mM))
        @select(:Additive, :Metabolite, :Patient, :Time, :absolute_mM, :diff_mM)
    end
    return diffed_df
end

# function pca_timeseries(long_df, additive)
#     wide_df = @chain long_df begin
#         @rsubset(:Additive == additive)
#         @rtransform(:Patient = :Sample[7:8])
#         @select(:Metabolite, :Patient, :Time, :absolute_mM)
#         @orderby(:Patient, :Time, :Metabolite)
#         unstack([:Patient, :Time], :Metabolite, :absolute_mM)
#         dropmissing()
#     end
#     patient_labels = wide_df[!, :Patient]
#     time_labels = wide_df[!, :Time]
#     X = Matrix(select(wide_df, Not([:Patient, :Time])))
#     zt = StatsBase.fit(StatsBase.ZScoreTransform, X, dims = 1)
#     Xzt = StatsBase.transform(zt, X)'
#     rows_with_nans = vec(any(isnan, Xzt, dims = 2))
#     display(rows_with_nans)
#     Xzt_no_nans = Xzt[.!(rows_with_nans), :]
#     M = fit(PCA, Xzt_no_nans; pratio = 0.9, mean = 0)
#     display(M)
#     Xzt_transform = MultivariateStats.predict(M, Xzt_no_nans)
#     println("size(X) ", size(X))
#     println("size(Xzt) ", size(Xzt))
#     println("size(Xzt_no_nans) ", size(Xzt_no_nans))
#     println("size(Xzt_transform) ", size(Xzt_transform))
# end

function pca_timeseries(long_df, additive)
    wide_df = @chain long_df begin
        @rsubset(:Additive == additive)
        @rtransform(:Patient = :Sample[7:8])
        @select(:Metabolite, :Patient, :Time, :absolute_mM)
        @groupby(:Time, :Metabolite)
        @combine(:mean_mM = mean(skipmissing(:absolute_mM)))
        unstack(:Time, :Metabolite, :mean_mM, combine = first)
        @orderby(:Time)
    end
    display(first(wide_df, 100))
    time_labels = wide_df.Time
    X = Matrix(select(wide_df, Not(:Time)))
    colmeans = map(eachcol(X)) do c
        m = mean(skipmissing(c))
        return isfinite(m) ? m : missing
    end
    for j in axes(X, 2)
        if ismissing(colmeans[j])
            continue
        end
        @inbounds for i in axes(X, 1)
            if ismissing(X[i, j])
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
    Xf = Xf[:, good_cols]
    if size(Xf, 2) == 0
        error("After filtering, no valid metabolite columns remain for PCA.")
    end
    zt = StatsBase.fit(StatsBase.ZScoreTransform, Xf; dims = 1)
    Xz = StatsBase.transform(zt, Xf)
    Xzt = copy(Xz')
    M = fit(PCA, Xzt; pratio = 0.9, mean = false)
    scores = MultivariateStats.predict(M, Xzt)
    display(M)
    return (
        model = M,
        scores = scores,
        time_labels = time_labels,
        kept_columns = findall(good_cols),
        wide_df = wide_df,
    )
end

function plot_pca_loadings(pca_result)
    kept_columns = pca_result.kept_columns
    wide_df = pca_result.wide_df
    M = pca_result.model
    L = loadings(M)
    pc1 = L[:, 1]
    pc2 = L[:, 2]
    metabolite_names = names(select(wide_df, Not(:Time)))[kept_columns]
    fig = Figure(resolution = (700, 600))
    ax = Axis(fig[1, 1],
        xlabel = "PC1 loading",
        ylabel = "PC2 loading",
        title = "PCA Loadings (Pattern Matrix)",
        aspect = DataAspect()
    )
    scatter!(ax, pc1, pc2, markersize = 12, color = :dodgerblue)
    for (x, y, name) in zip(pc1, pc2, metabolite_names)
        text!(ax, x, y, text = name, offset = (5, 5), align = (:left, :bottom))
    end
    hlines!(ax, [0.0], color = (:gray, 0.4), linewidth = 1)
    vlines!(ax, [0.0], color = (:gray, 0.4), linewidth = 1)
    return fig
end

function additive_metabolite_time_points(long_df, additive, metabolite, tf)
    result_df = @chain long_df begin
        @rsubset(:Additive == additive, :Metabolite == metabolite, :Time >= tf - 1, :Time <= tf)
        @select(:Time, :absolute_mM)
        @orderby(:Time)
    end
    return result_df
end

function regress_concentration_dxdt(long_df, bootstrap_reps)
    Random.seed!(123)
    unique_additives = unique(long_df.Additive)
    unique_metabolites = unique(long_df.Metabolite)
    final_times = [2, 3, 4, 5, 6]
    rows = []
    for (metabolite, additive, final_time) in
        product(unique_metabolites, unique_additives, final_times)
        println("Calculating $additive, $metabolite, $final_time")
        regression_df =
            additive_metabolite_time_points(long_df, additive, metabolite, final_time)
        n = nrow(regression_df)
        slopes = ThreadsX.map(1:bootstrap_reps) do _
            sample_idx = rand(1:n, n)
            boot_df = regression_df[sample_idx, :]
            boot_model = lm(@formula(absolute_mM ~ Time), boot_df)
            boot_coefs = coef(boot_model)
            boot_coefs[2]
        end
        lower_bound = quantile(slopes, 0.025)
        upper_bound = quantile(slopes, 0.975)
        mean_rate = mean(slopes)
        rate_skew = skewness(slopes)
        single_model = lm(@formula(absolute_mM ~ Time), regression_df)
        coefs = coef(single_model)
        single_rate = coefs[2]
        ci = confint(single_model)
        single_lb = ci[2, 1]
        single_ub = ci[2, 2]
        row = (
            additive = additive,
            metabolite = metabolite,
            final_time = final_time,
            mean_rate = mean_rate,
            skew = rate_skew,
            lower_bound = lower_bound,
            upper_bound = upper_bound,
            single_rate = single_rate,
            single_lb = single_lb,
            single_ub = single_ub,
        )
        push!(rows, row)
    end
    result_df = @chain rows begin
        DataFrame()
        @orderby(:additive, :metabolite, :final_time)
    end
    return result_df
end

end

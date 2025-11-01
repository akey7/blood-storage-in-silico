module AbsoluteQuant

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

export load_absolute_quant,
    load_relative_quant,
    combine_relative_and_absolute_quant,
    cluster_all_additives_all_n_clusters,
    plot_elbows,
    plot_c_means_all_additives,
    plot_all_mmol_per_L_timeseries,
    diff_mmol_per_L

function load_absolute_quant()
    absolute_filename = joinpath("input", "Absolute Quant Data Sheet.xlsx")
    cells_day_1_df = DataFrame(XLSX.readtable(absolute_filename, "cells_day_1"))
    absolute_metabolite_ids = DataFrame(XLSX.readtable(absolute_filename, "metabolite_ids"))
    absolute_quant_df = @chain cells_day_1_df begin
        stack(Not(:id), variable_name = :mixed_name, value_name = :mmol_per_L)
        @rtransform(:sample_set = split(:id, "_")[2])
        innerjoin(absolute_metabolite_ids, on = :mixed_name => :MixedName)
        @rtransform(:prop_mmol_per_L = :mmol_per_L * :Proportion)
        select([:sample_set, :id, :Metabolite, :prop_mmol_per_L])
    end
    absolute_quant_medians_df = @by absolute_quant_df :Metabolite begin
        :median_prop_mmol_per_L = median(:prop_mmol_per_L)
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
        @rtransform(:relative_mmol_per_L = :FoldChange * :median_prop_mmol_per_L)
        @orderby(:Additive, :Time, :Metabolite)
    end
    wide_df = @chain long_df begin
        @select(:Sample, :Time, :Additive, :Metabolite, :relative_mmol_per_L)
        unstack(
            [:Sample, :Time, :Additive],
            :Metabolite,
            :relative_mmol_per_L,
            combine = first,
        )
        @orderby(:Additive, :Time)
    end
    return long_df, wide_df
end

function prepare_long_df_for_clustering(long_df, additive)
    long_df_2 = deepcopy(long_df)
    wide_timeseries_df = @chain long_df_2 begin
        @rsubset(:Additive == additive, !isapprox(:relative_mmol_per_L, 0.0))
        @select(:Metabolite, :Time, :relative_mmol_per_L)
        @orderby(:Metabolite, :Time)
        unstack(:Metabolite, :Time, :relative_mmol_per_L, combine = first)
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
        @select(:primary_cluster, :Patient, :Time, :Metabolite, :relative_mmol_per_L)
        unstack(
            [:primary_cluster, :Patient, :Time],
            :Metabolite,
            :relative_mmol_per_L,
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
        value_name = :standardized_mmol_per_L,
    )
    time_points = unique(plt_df.Time)
    plt =
        data(plt_df) *
        mapping(
            :Time,
            :standardized_mmol_per_L => "standardized mmol/L",
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
end

function plot_c_means_all_additives(long_df, all_memberships_dfs, n_clusters)
    additives = unique(long_df.Additive)
    for additive in additives
        plot_c_means_for_additive_and_n_clusters(
            long_df,
            all_memberships_dfs,
            additive,
            n_clusters,
        )
    end
end

function plot_all_mmol_per_L_timeseries(long_df)
    agg_df = @chain long_df begin
        @groupby(:Additive, :Metabolite, :Time)
        @combine(:median_mmol_per_L = median(skipmissing(:relative_mmol_per_L)))
        @orderby(:Additive, :Metabolite, :Time)
    end
    metabolites = unique(agg_df.Metabolite)
    time_points = unique(agg_df.Time)
    for metabolite in metabolites
        clean_metabolite = replace(metabolite, r"[^A-Za-z0-9_]" => "_")
        filename = joinpath("output", "relative_absolute_plots", "$(clean_metabolite).png")
        plt_df = @rsubset(agg_df, :Metabolite == metabolite)
        plt =
            data(plt_df) *
            mapping(:Time, :median_mmol_per_L => "Median mmol/L", color = :Additive) *
            (visual(Lines) + visual(Scatter; markersize = 10))
        fig = draw(
            plt;
            figure = (; size = (750, 500)),
            axis = (; title = metabolite, xticks = time_points),
        )
        save(filename, fig)
        println("Wrote $filename")
    end
end

function diff_mmol_per_L(long_df)
    diffed_df = @chain long_df begin
        @rtransform(:Patient = :Sample[7:8])
        @orderby(:Additive, :Metabolite, :Patient, :Time)
        @groupby(:Additive, :Metabolite, :Patient)
        @transform(
            :diff_mmol_per_L =
                :relative_mmol_per_L .- ShiftedArrays.lag(:relative_mmol_per_L)
        )
        @select(
            :Additive,
            :Metabolite,
            :Patient,
            :Time,
            :relative_mmol_per_L,
            :diff_mmol_per_L
        )
    end
    return diffed_df
end

end

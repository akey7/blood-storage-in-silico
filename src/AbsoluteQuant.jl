module AbsoluteQuant

using Base.Iterators
using CSV
using XLSX
using DataFrames
using DataFramesMeta
using Statistics
using StatsBase
using AlgebraOfGraphics
using CairoMakie
using Makie
using CategoricalArrays
using Clustering
using COBREXA
import JSONFBCModels
using ColorSchemes

export load_absolute_quant,
    load_relative_quant, combine_relative_and_absolute_quant, cluster_all_additives

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
        @subset(:Additive .== additive)
        @select(:Metabolite, :Time, :relative_mmol_per_L)
        @orderby(:Metabolite, :Time)
        unstack(:Metabolite, :Time, :relative_mmol_per_L, combine = first)
    end
    return wide_timeseries_df
end

function calc_fuzzy_objective(result, X, μ = 2.0)
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
    μ = 2.0,
)
    X = Matrix{Float64}(disallowmissing(wide_timeseries_df[:, Not(:Metabolite)]))
    result = fuzzy_cmeans(X', n_clusters, μ, maxiter = 200, display = :iter)
    weights_col_names = string.(axes(result.weights, 2))
    memberships_df = DataFrame(result.weights, weights_col_names)
    memberships_df.Metabolite = wide_timeseries_df.Metabolite
    memberships_df[!, :Additive] .= additive
    memberships_df[!, :NClusters] .= n_clusters
    fuzzy_objective = calc_fuzzy_objective(result, X, μ)
    return memberships_df, fuzzy_objective
end

function cluster_all_additives(long_df; n_clusters = 7)
    additives = unique(long_df.Additive)
    for additive in additives
        println(uppercase(additive))
        memberships_df = @chain long_df begin
            prepare_long_df_for_clustering(additive)
            c_means_metabolite_trajectories(additive = additive, n_clusters = n_clusters)
        end
        println(first(memberships_df, 10))
    end
end

end

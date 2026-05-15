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
using ColorSchemes
using Clustering
using Distances
using ShiftedArrays
using MultivariateStats
using GLM
using StatsModels
using Statistics
using Random
using ThreadsX
using ProgressMeter

export load_absolute_quant,
    load_relative_quant,
    combine_relative_and_absolute_quant_c,
    cluster_all_additives_all_n_clusters,
    plot_elbows,
    plot_c_means_all_additives,
    plot_all_mM_timeseries,
    diff_mM,
    pca_timeseries,
    regress_concentration_vs_time,
    plot_pca_all_additives,
    plot_all_regressions,
    qc,
    load_extracellular_absolute_quant,
    combine_relative_and_absolute_quant_e,
    union_and_pivot_wide,
    load_relative_quant_2

"""
    load_absolute_quant()

Loads the aboslute quantification data from the Excel sheet at "input/Absolute Quant Data Sheet.xlsx". Loads the proportination data from a different sheet in the same workbook and maps individual peaks in the absolute quant data to proportional concentrations for further analysis. Performs initial mapping of absolute quant data in terms of single metabolite ids.

# Returns
`Tuple{DataFrame,DataFrame}`

Returns two dataframes: `absolute_quant_df`, which contains the absolute quant information of each sample, and `absolute_quant_medians_df` which is the median concentration value for each compound aggregated across all samples.

The first DataFrame has the following columns
1. `:sample_set`: The set of samples the row is from
2. `:id`: The sample id.
3. `:Metabolite`: The metabolite id of the compound
4. `:prop_mM`: The proportionated concentration of the row in mM

The second DataFrame has the following columns
1. `:Metabolite`: A metabolite id
2. `:media_prop_mM`: The median concentration for that metabolite id.
"""
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

"""
    load_extracellular_absolute_quant()

Loads the absolute quant data for extracellular meatbolites from the `Sheet1` sheet of `input/Absolute Quant Extracellular Datasheet.xlsx`.

# Returns
`DataFrame`

Returns a DataFrame with the following columns
1. `:metabolite_id`: The extracellular metabolite id
2. `:median_prop_mM`: Median proportinated absolute concentration.
"""
function load_extracellular_absolute_quant()
    filename = joinpath("input", "Absolute Quant Extracellular Datasheet.xlsx")
    df = DataFrame(XLSX.readtable(filename, "Sheet1"))
    select_df = @select(df, :metabolite_id, :median_prop_mM)
    return select_df
end

"""
    load_relative_quant()

Loads relative quant data from the file "input/Data Sheet 1.CSV" from Nemkov et al (2022). Also loads the proporination sheet from "input/Proportionation Sheet 2.csv" which maps combined human-friendly metabolite names to individual metabolite ids with fractions of abundances assigned to each individual id. Calculates fold changes relative to the median of 01-Ctrl AS3, Week 1 measurement for each metabolite.

# Returns
`DataFrame`

Returns a DataFrame (pivoted from wide to long) with the following columns:
1. `:Sample`: Sample id that include, among other things, time point and patient.
2. `:Time`: Time point of measurement in weeks.
3. `:Metabolite`: Metabolite id of measurement
4. `:FoldChange`: Fold change over the median control measurement at week 1 for that metabolite.
"""
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

function load_relative_quant_2()
    relative_filename = joinpath("input", "AS Dev Library Trial 1.csv")
    wide_df = CSV.read(relative_filename, DataFrame)
    # proportination_filename = joinpath("input", "Proportionation Sheet 2.csv")
    # proportination_df = CSV.read(proportination_filename, DataFrame)
    long_df_1 = stack(
        wide_df,
        Not([:Sample, :Day, :Condition]),
        variable_name = :MixedName,
        value_name = :Intensity,
    )
    valid_measurements_df = @chain long_df_1 begin
        @rtransform(:is_valid = !(ismissing(:Intensity) || isapprox(:Intensity, 0.0)))
        @orderby(:Condition, :Day, :MixedName)
        @select(:Condition, :Day, :MixedName, :Intensity, :Sample, :is_valid)
    end
    long_df_2 = @chain long_df_1 begin
        @rtransform(:Time = div(:Day, 7, RoundUp))
        @select(:Additive = :Condition, :Time, :MixedName, :Intensity)
    end
    as3_intensity_df = @rsubset(long_df_2, :Additive == "AS3")
    as3_fold_change_df_1 = @chain as3_intensity_df begin
        @groupby(:Additive, :Time, :MixedName)
        @combine(
            :median_as3_intensity = median(:Intensity),
            :max_as3_intensity = maximum(:Intensity),
            :min_as3_intensity = minimum(:Intensity)
        )
        @rtransform(
            :max_as3_fold_change = :max_as3_intensity / :median_as3_intensity,
            :min_as3_fold_change = :min_as3_intensity / :median_as3_intensity
        )
        @orderby(:Additive, :Time, :MixedName)
        @select(
            :Time,
            :MixedName,
            :median_as3_intensity,
            :min_as3_fold_change,
            :max_as3_fold_change
        )
    end
    as3_fold_change_df_2 = @rsubset as3_fold_change_df_1 begin
        isfinite(:min_as3_fold_change)
        isfinite(:max_as3_fold_change)
    end
    non_as3_conditions_df = @rsubset(long_df_2, :Additive != "AS3", isfinite(:Intensity))
    non_as3_conditions_fold_change_df = @chain non_as3_conditions_df begin
        innerjoin(as3_fold_change_df_2; on = [:Time, :MixedName])
        @rtransform(
            :min_intensity = :Intensity * :min_as3_fold_change,
            :max_intensity = :Intensity * :max_as3_fold_change
        )
        @orderby(:Additive, :Time, :MixedName, :Intensity, :min_intensity, :max_intensity)
        @select(:Additive, :Time, :MixedName, :Intensity, :min_intensity, :max_intensity)
    end
    non_as3_expanded_long_rows = []
    for row in eachrow(non_as3_conditions_fold_change_df)
        row_1 = (
            Additive = row.Additive,
            Time = row.Time,
            MixedName = row.MixedName,
            Intensity = row.min_intensity,
        )
        row_2 = (
            Additive = row.Additive,
            Time = row.Time,
            MixedName = row.MixedName,
            Intensity = row.Intensity,
        )
        row_3 = (
            Additive = row.Additive,
            Time = row.Time,
            MixedName = row.MixedName,
            Intensity = row.max_intensity,
        )
        push!(non_as3_expanded_long_rows, row_1)
        push!(non_as3_expanded_long_rows, row_2)
        push!(non_as3_expanded_long_rows, row_3)
    end
    non_as3_expanded_long_df = DataFrame(non_as3_expanded_long_rows)
    display(first(non_as3_expanded_long_df, 10))
    all_conditions_long_df_1 = vcat(as3_intensity_df, non_as3_expanded_long_df)
    all_conditions_long_df_2 =
        @orderby(all_conditions_long_df_1, :Additive, :Time, :MixedName, :Intensity)
    result = (
        valid_measurements_df = valid_measurements_df,
        as3_fold_change_df = as3_fold_change_df_2,
        conditions_fold_change_df = non_as3_conditions_fold_change_df,
        all_conditions_long_df = all_conditions_long_df_2,
    )
    return result
end

"""
    qc(fold_changes_df, patient_count = 6)

A quality check function. Given the `fold_changes_df` returned by [`load_relative_quant`](@ref BloodStorageInSilico.AbsoluteQuant.load_relative_quant) ensures a consistent count of observations for each time, metabolite, additive combinations by coubnting the number of patients for each.

# Arguments
1. `fold_chnages_df`: `fold_changes_df` returned by [`load_relative_quant`](@ref BloodStorageInSilico.AbsoluteQuant.load_relative_quant)
2. `patient_count = 6`: Number of patients that should be represented at each time, metabolite, additive

# Returns
`Tuple{DataFrame,DataFrame}`

Returns two DataFrames:

1. The first DataFrame contains each time, additive, metabolite that does NOT have the number of patients specified by `patient_count`.
2. The second DataFrame contains the count of fold changes that are approximately 0.0 for each time, additive, metabolite.
"""
function qc(fold_changes_df, patient_count = 6)
    qc_fold_change_counts_df = @chain fold_changes_df begin
        @groupby(:Time, :Additive, :Metabolite)
        combine(nrow => :Count)
        @rsubset(:Count != patient_count)
        sort(:Count)
    end
    qc_fold_change_zeros_df = @chain fold_changes_df begin
        @rsubset(isapprox(:FoldChange, 0.0))
        @groupby(:Time, :Additive, :Metabolite)
        combine(nrow => :Count)
        sort(:Count)
    end
    return qc_fold_change_counts_df, qc_fold_change_zeros_df
end

"""
    combine_relative_and_absolute_quant_c(fold_changes_df, absolute_quant_medians_df)

In this function, the relative quant and absolute quant data for the cytosolic metabolites (metabolites ending with `_c`) are combined to approximate aboslute quantification to put into models.

# Arguments
1. `fold_changes_df`: The relative quant data from [`load_relative_quant`](@ref BloodStorageInSilico.AbsoluteQuant.load_relative_quant)
2. `absolute_quant_medians_df`: The absolute quant data from [`load_absolute_quant`](@ref BloodStorageInSilico.AbsoluteQuant.load_absolute_quant)

# Returns
`DataFrame`

Returns a long DataFrame that contains the following columns and sorted by `:Additive`, `:Time`, and `:Metabolite`:
1. `:Sample`: Sample id
2. `:Time`: Measurement time (in weeks)
3. `:Additive`: Additive
4. `:Metabolite`: Metabolite id
5. `:FoldChange`: Fold change from absolute quant
6. `:absolute_mM`: The approximate mM concentration of that metabolite for that row
"""
function combine_relative_and_absolute_quant_c(fold_changes_df, absolute_quant_medians_df)
    long_df = @chain fold_changes_df begin
        innerjoin(absolute_quant_medians_df, on = :Metabolite)
        @rtransform(:absolute_mM = :FoldChange * :median_prop_mM)
        @orderby(:Additive, :Time, :Metabolite)
        @select(:Sample, :Time, :Additive, :Metabolite, :FoldChange, :absolute_mM)
    end
    return long_df
end

"""
    combine_relative_and_absolute_quant_e(fold_changes_df, absolute_extracellular_quant_df)

In this function, the relative quant and absolute quant data for the extracellular metabolites (metabolites ending with `_e`) are combined to approximate aboslute quantification to put into models.

Since the relative quant fold changes DataFrame is specified by cytosolic metabolite ids, these metabolite ids are transformed to extracellular metabolite ids during the join.

# Returns
`DataFrame`

Returns a long DataFrame that contains the following columns and sorted by `:Additive`, `:Time`, and `:Metabolite`:
1. `:Sample`: Sample id
2. `:Time`: Measurement time (in weeks)
3. `:Additive`: Additive
4. `:Metabolite`: Metabolite id
5. `:FoldChange`: Fold change from absolute quant
6. `:absolute_mM`: The approximate mM concentration of that metabolite for that row
"""
function combine_relative_and_absolute_quant_e(
    fold_changes_df,
    absolute_extracellular_quant_df,
)
    quant_e_df = @chain absolute_extracellular_quant_df begin
        rename(:metabolite_id => :Metabolite)
        @rtransform(:compound = replace(:Metabolite, "_e" => ""))
    end
    fold_changes_compounds_df =
        @rtransform(fold_changes_df, :compound = replace(:Metabolite, "_c" => ""))
    long_df = @chain quant_e_df begin
        innerjoin(fold_changes_compounds_df, on = :compound, makeunique = true)
        @rtransform(:absolute_mM = :FoldChange * :median_prop_mM)
        select(Not(:Metabolite_1, :compound))
        @select(:Sample, :Time, :Additive, :Metabolite, :FoldChange, :absolute_mM)
        @orderby(:Additive, :Time, :Metabolite)
    end
    return long_df
end

"""
    union_and_pivot_wide(absolute_quant_c_long_df, absolute_quant_e_long_df; include_extracellular = true)

Unions the cytosolic and extracellular absolute quant long DataFrames. Returns the union of the long DataFrame and a wide pivoted DataFrame. Optionally exclude extracellular metabolites.

See also: [`combine_relative_and_absolute_quant_c`](@ref BloodStorageInSilico.AbsoluteQuant.union_and_pivot_wide) and [`combine_relative_and_absolute_quant_e`](@ref BloodStorageInSilico.AbsoluteQuant.union_and_pivot_wide)

# Arguments
1. `absolute_quant_c_long_df`: Cytosolic absolute quant long DataFrame
2. `absolute_quant_e_long_df`: Extracellular absolute quant long DataFrame
3. `include_extracellular = true`: If `true` (the default), extracellular metabolites are included in the result. If `false`, extracellular metabolites are excluded from the results.

# Returns
`Tuple{DataFrame,DataFrame}`

First element of the tuple is a long DataFrame with the following columns:
1. `:Sample`: Sample id
2. `:Time`: Measurement time (in weeks)
3. `:Additive`: Additive
4. `:Metabolite`: Metabolite id
5. `:FoldChange`: Fold change from absolute quant
6. `:absolute_mM`: The approximate mM concentration of that metabolite for that row

Second element of the tuple is a wide DataFrame that contains the following columns and is sorted by `:Additive`, `:Time`, and `:Sample`:
1. `:Sample`: Sample id of that row
2. `:Time`: Time in weeks of that observation
3. `:Additive`: Additive
4. A subsequent column for each metabolite
"""
function union_and_pivot_wide(
    absolute_quant_c_long_df,
    absolute_quant_e_long_df;
    include_extracellular = true,
)
    union_long_df =
        include_extracellular ? vcat(absolute_quant_c_long_df, absolute_quant_e_long_df) :
        absolute_quant_c_long_df
    union_wide_df = @chain union_long_df begin
        @select(:Sample, :Time, :Additive, :Metabolite, :absolute_mM)
        unstack([:Sample, :Time, :Additive], :Metabolite, :absolute_mM, combine = first)
        @orderby(:Additive, :Time, :Sample)
    end
    return union_long_df, union_wide_df
end

"""
    prepare_long_df_for_clustering(long_df, additive)

Prepare the long DataFrame from [`union_and_pivot_wide`](@ref BloodStorageInSilico.AbsoluteQuant.union_and_pivot_wide) for timeseries analysis by c-means clustering. The original long DataFrame is filtered down to a single additive. If a measurement is duplicated, the first conflicting measurement will be used and zero values are excluded. The result is a wide DataFrame of timeseries, with a column for each timepoint, that can be used for clustering.

# Arguments
1. `long_df`: The long DataFrame to pivot.
2. `additive`: The additive to make time series for.

# Returns
`DataFrame`

Returns a wide DataFrame with the following columns
1. `:Metabolite`: The metabolite for which the row is a time series.
2. Subsequent columns: A column for each timepoint measured for that metabolite.
"""
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

"""
    calc_fuzzy_objective(result, X, μ)

This function is used to make elbow plots to determine the optimal number of clusters for the clustering results. It is a reproduction of the objective function of c-means clustering.

[This is the objective function](https://juliastats.org/Clustering.jl/stable/fuzzycmeans.html#fuzzy_cmeans_def) calculated here.

# Arguments
1. `result`: C-means clustering resuly with centers and weights to evaluate.
2. `X`: The matrix of features that were clustered.
3. `μ`: The fuzziness factor that was used for clustering

# Returns
`Float64`

Returns the value of the c-means objective.
"""
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

"""
    c_means_metabolite_trajectories(wide_timeseries_df; additive = "01-Ctrl AS3", n_clusters = 5, μ = 5.0)

Calculates the c-means clusters of the metabolite time series in the given wide DataFrame. Cleans the DataFrame before clustering to remove problematic values, such as constant features, NaNs and Infs. Sensible clustering conditions are enforced with `@assert`. Prints diagnostic logging messages durinng operation. Uses CityBlock distances because proved to create more robust results.

# Arguments
1. `wide_timeseries_df`: DataFrame created by [`prepare_long_df_for_clustering`](@ref BloodStorageInSilico.AbsoluteQuant.prepare_long_df_for_clustering) for timeseries analysis.
2. `additive`: Additive to use when setting the `:Additive` column of the final output. NOTE: This does not affect `wide_timeseries_df`, which assumed to be filtered before being passed to this function. Rather, this argument just affects the OUTPUT DataFrame so it can be combined with other clustering runs later. Defaults to "01-Ctrl AS3".
3. `n_clusters`: The number of clusters to split the metabolites into. Defaults to 5
4. `μ`: The fuzziness factor for the clustering. Defaults to 5.0

# Returns
`Tuple{DataFrame,Float64,Matrix}`

First, returns a DataFrame with the clustering results that has the following columns:
1. `:Metabolite`: The metabolite
2. `:Additive`: The additive the source data is from
3. `:NClusters`: The number of clusters the data was split into for this run
4. Other columns: A column with the membership weight in each cluster for that metabolite.

Second, returns the fuzzy objective value for making an elbow plot as calculated by [`calc_fuzzy_objective`](@ref BloodStorageInSilico.AbsoluteQuant.calc_fuzzy_objective).

Third, returns the Matrix used for the clustering.
"""
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

"""
    cluster_all_additives_all_n_clusters(long_df; max_clusters = 10)

To make a complete clustering analysis of this dataset, clustering must be performed for each additive, different numbers of clusters must be attempted, and the results need to be aggregated to make figures. This function iterates through all additives and numbers of clusters to aggregate all of these runs into one place for further analysis.

# Arguments
1. `long_df`: The long DataFrame from [`union_and_pivot_wide`](@ref BloodStorageInSilico.AbsoluteQuant.union_and_pivot_wide) that is the source of the data to be clustered.
2. `max_clusters = 10`: Defaults to 10. For example, if left at 10, clustering into 2, 3, 4, 5, 6, 7, 8, 9, and 10 clusters will be attempted for selection of the optimal number of clusters.

# Returns
`Tuple{Vector{DataFrame},DataFrame}`

The first element of the tuple is a Vector of DataFrames. The contents are the membership weights DataFrames resulting from calling [`c_means_metabolite_trajectories`](@ref BloodStorageInSilico.AbsoluteQuant.c_means_metabolite_trajectories) and taking the first DataFrame of the resulting Tuple.

The second element of the Tuple is a DataFrame that is the minimized objective value for each number of clusters. The columns of this DataFrame are:
1. `:additive`: The additive
2. `:n_clusters`: The number of clusters the data were split into.
3. `:fuzzy_objective`: The minimized objective value for that number of clusters, suitable for making an elbow plot.
"""
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

"""
    plot_elbows(fuzzy_objectives_df)

Creates and saves the elbow plots for each additive to determine the optimal number of clusters based on fuzzy objective values created by [`cluster_all_additives_all_n_clusters`](@ref BloodStorageInSilico.AbsoluteQuant.cluster_all_additives_all_n_clusters).
    
The final output is saved to `output/relative_absolute_c_means/elbows.png`

# Arguments
1. `fuzzy_objectives_df`: DataFrame of objective values
"""
function plot_elbows(fuzzy_objectives_df)
    xticks = unique(fuzzy_objectives_df.n_clusters)
    plt =
        data(fuzzy_objectives_df) *
        mapping(
            :n_clusters => "N Clusters",
            :fuzzy_objective => "Fuzzy Objective",
            row = :additive,
        ) *
        visual(Lines)
    figure_options = (; size = (300, 700), title = "Objective Elbows")
    axis_options = (; xticks = xticks)
    facet_options = (; linkxaxes = :all, linkyaxes = :minimal)
    fig = draw(plt; figure = figure_options, axis = axis_options, facet = facet_options)
    fig_filename = joinpath("output", "relative_absolute_c_means", "elbows.png")
    save(fig_filename, fig)
    println("Wrote $fig_filename")
end

"""
    plot_c_means_for_additive_and_n_clusters(long_df, all_memberships_dfs, additive, n_clusters)

Clusters metabolite timeline trajectories in the given additive into the given number of clusters, standardizes the concentrations, saves a plot to the `output/relative_absolute_c_means` folder, and returns the memberships DataFrame that made the plot. The requested number of clusters and additive must be in the data passed.

# Arguments
1. `long_df`: The long DataFrame from [`union_and_pivot_wide`](@ref BloodStorageInSilico.AbsoluteQuant.union_and_pivot_wide) that is the source of the data to be clustered.
2. `all_memberships_dfs`: Vector of DataFrames for clustering into various numbers of clusters as returned by [`cluster_all_additives_all_n_clusters`](@ref BloodStorageInSilico.AbsoluteQuant.cluster_all_additives_all_n_clusters)
3. `additive`: Additive for the clustering and plotting.
4. `n_clusters`: Number of clusters to plot the trajectories into.

# Returns
`DataFrame`

Returns the cluster membership DataFrame that was plotted. This is useful for manual inspection if needed.
"""
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
    stacked_df = stack(
        standardization_df,
        Not([:primary_cluster, :Patient, :Time]),
        variable_name = :Metabolite,
        value_name = :standardized_mM,
    )
    plt_df = @rtransform(stacked_df, :cluster_label = "Cluster $(:primary_cluster)")
    time_points = unique(plt_df.Time)
    plt =
        data(plt_df) *
        mapping(
            :Time => "Time (weeks)",
            :standardized_mM => "standardized mM",
            row = :cluster_label,
            group = :Metabolite,
            color = :cluster_label,
        ) *
        visual(Lines) *
        visual(alpha = 0.3)
    figure_options = (; size = (300, 700), title = additive)
    base_palettes = Dict(
        "01-Ctrl AS3" => ColorSchemes.devon,
        "02-Adenosine" => ColorSchemes.buda,
        "03-Glutamine" => ColorSchemes.berlin,
        "04-Methionine" => ColorSchemes.batlow,
        "07-NAC" => ColorSchemes.acton,
        "08-Taurine" => ColorSchemes.bamako,
    )
    cluster_palette = get(base_palettes[additive], range(0, 0.6, length = n_clusters))
    fig = draw(
        plt,
        scales(Color = (; legend = false, palette = cluster_palette));
        figure = figure_options,
        axis = (; xticks = time_points, limits = (nothing, nothing, -10.0, 10.0)),
        facet = (; linkxaxes = :all, linkyaxes = :all),
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

"""
    plot_c_means_all_additives(long_df, all_memberships_dfs, n_clusters)

Iterates through all additives and calls [`plot_c_means_for_additive_and_n_clusters`](@ref BloodStorageInSilico.AbsoluteQuant.plot_c_means_for_additive_and_n_clusters) to make a plot for each additive with the given number of clusters.

# Arguments
1. `long_df`: The long DataFrame from [`union_and_pivot_wide`](@ref BloodStorageInSilico.AbsoluteQuant.union_and_pivot_wide) that is the source of the data to be clustered.
2. `all_memberships_dfs`: Vector of DataFrames for clustering into various numbers of clusters as returned by [`cluster_all_additives_all_n_clusters`](@ref BloodStorageInSilico.AbsoluteQuant.cluster_all_additives_all_n_clusters)
3. `n_clusters`: Number of clusters to plot the trajectories into.

# Returns
`DataFrame`

Consolidated DataFrame of all primary cluster assignments for all metabolites in all additives.
"""
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

"""
    plot_all_mM_timeseries(long_df)

Plots the absolute quant approximations for all metabolites in all additives. Saves each plot to the `output/relative_absolute_plots` folder as it goes. Dislpays a nifty progress bar as it writes plots.

# Arguments
1. `long_df`: The long DataFrame from [`union_and_pivot_wide`](@ref BloodStorageInSilico.AbsoluteQuant.union_and_pivot_wide). The source of the data that will be plotted.
"""
function plot_all_mM_timeseries(long_df)
    agg_df = @chain long_df begin
        @groupby(:Additive, :Metabolite, :Time)
        @combine(:median_mM = median(skipmissing(:absolute_mM)))
        @orderby(:Additive, :Metabolite, :Time)
    end
    metabolites = unique(agg_df.Metabolite)
    time_points = unique(agg_df.Time)
    additives = sort(unique(agg_df.Additive))
    cluster_palette = Dict(
        "01-Ctrl AS3" => ColorSchemes.devon10[5],
        "02-Adenosine" => ColorSchemes.buda10[5],
        "03-Glutamine" => ColorSchemes.berlin10[5],
        "04-Methionine" => ColorSchemes.batlow10[5],
        "07-NAC" => ColorSchemes.acton10[5],
        "08-Taurine" => ColorSchemes.bamako10[5],
    )
    pal_vec = [cluster_palette[a] for a in additives]
    n_plots = length(metabolites)
    prog = Progress(n_plots, "Writing mM timeseries plots")
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
            mapping(
                :Time => "Time (week)",
                :absolute_mM => "Concentration (mM)",
                color = :Additive,
                marker = :Additive,
            ) *
            visual(Scatter, markersize = 14, alpha = 0.5)
        plt = line_plt + scatter_plt
        fig = draw(
            plt,
            scales(Color = (; palette = pal_vec));
            figure = (; size = (750, 500)),
            axis = (; title = metabolite, xticks = time_points),
        )
        save(filename, fig)
        next!(prog)
    end
end

"""
    additive_metabolite_time_points(long_df, additive, metabolite, tf)

Used by [`regress_concentration_vs_time`](@ref BloodStorageInSilico.AbsoluteQuant.regress_concentration_vs_time) to get a time course for a particular metabolite in a specified additive.

# Arguments
1. `long_df`: The long DataFrame from [`union_and_pivot_wide`](@ref BloodStorageInSilico.AbsoluteQuant.union_and_pivot_wide).
2. `additive`: Additive to select.
3. `metabolite`: Metabolite id to select
4. `tf`: Final time point (2, 3, 4, 5, 6) to select

# Returns
`DataFrame`

Each time point with the approximated mM concentration, ordered by time.
"""
function additive_metabolite_time_points(long_df, additive, metabolite, tf)
    result_df = @chain long_df begin
        @rsubset(:Additive == additive, :Metabolite == metabolite, :Time >= tf - 1, :Time <= tf)
        @select(:Time, :absolute_mM)
        @orderby(:Time)
    end
    return result_df
end

"""
    regress_concentration_vs_time(long_df)

Regresses the concentration vs time to find the rate of metabolite concentration change (95% confidence interval upper and lower bounds) for all the metabolites and additives in `long_df`. Uses ThreadsX to split this task into multiple threads if multiple threads are available.

# Arguments
1. `long_df`: The long DataFrame from [`union_and_pivot_wide`](@ref BloodStorageInSilico.AbsoluteQuant.union_and_pivot_wide).

# Returns
`DataFrame`

Returns a DataFrame with the concentration rate regression results. The DataFrame is sorted by additive, final_time, and metabolite. The following columns are available:

1. `additive`: The additive the data for the regression is from.
2. `metabolite`: Metabolite id of the row.
3. `final_time`: Final time point of the regression. The timespan of the regression is final_time - 1 to final_time.
4. `intercept`: Intercept of the regression
5. `rate`: Slope of the regression.
6. `lb`: Lower bound of the 95% confidence interval of the slope.
7. `ub`: Upper bound of the 95% confidence interval of the slope.
8. `lb_ub_different_sign`: `true` if the lower bound and upper bound are different signs.
"""
function regress_concentration_vs_time(long_df)
    Random.seed!(123)
    unique_additives = unique(long_df.Additive)
    unique_metabolites = unique(long_df.Metabolite)
    final_times = [2, 3, 4, 5, 6]
    tasks = product(unique_metabolites, unique_additives, final_times)
    rows = ThreadsX.map(tasks) do t
        metabolite, additive, final_time = t
        println("Calculating $additive, $metabolite, $final_time")
        regression_df =
            additive_metabolite_time_points(long_df, additive, metabolite, final_time)
        single_model = lm(@formula(absolute_mM ~ Time), regression_df)
        coefs = coef(single_model)
        intercept = coefs[1]
        rate = coefs[2]
        ci = confint(single_model)
        lb = ci[2, 1]
        ub = ci[2, 2]
        (
            additive = additive,
            metabolite = metabolite,
            final_time = final_time,
            intercept = intercept,
            rate = rate,
            lb = lb,
            ub = ub,
        )
    end
    regressions_df = @chain rows begin
        DataFrame()
        @rtransform(:lb_ub_different_sign = sign(:lb) != sign(:ub))
        @orderby(:additive, :final_time, :metabolite)
        @select(
            :additive,
            :final_time,
            :metabolite,
            :intercept,
            :rate,
            :lb,
            :ub,
            :lb_ub_different_sign
        )
    end
    return regressions_df
end

"""
    plot_all_regressions(long_df)

Uses [`plot_regression`](@ref BloodStorageInSilico.AbsoluteQuant.plot_regression) for all metabolites in all additives to plot regression results. Saves plots to `output/regression_plots`. Displays a nifty progress bar as it writes plots.

# Arguments
1. `long_df`: The long DataFrame from [`union_and_pivot_wide`](@ref BloodStorageInSilico.AbsoluteQuant.union_and_pivot_wide).
"""
function plot_all_regressions(long_df)
    additives = unique(long_df.Additive)
    metabolites = unique(long_df.Metabolite)
    pairs = product(additives, metabolites)
    n_pairs = length(pairs)
    prog = Progress(n_pairs, "Writing regression plots")
    for (additive, metabolite) in pairs
        filename = joinpath("output", "regression_plots", "$additive $metabolite.png")
        fig = plot_regression(long_df, additive, metabolite)
        save(filename, fig)
        next!(prog)
    end
end

"""
    plot_regression(long_df, additive, metabolite)

# Arguments
1. `long_df`: The long DataFrame from [`union_and_pivot_wide`](@ref BloodStorageInSilico.AbsoluteQuant.union_and_pivot_wide).
2. `additive`: Additive
3. `metabolite`: Metabolite id

# Returns
`Figure`

Returns a Makie figure that can be displayed or saved.
"""
function plot_regression(long_df, additive, metabolite)
    super_title = "$additive $metabolite"
    fig = Figure(; size = (360, 720))
    Label(fig[0, :], text = super_title, fontsize = 25)
    final_times_to_figure_map =
        Dict(2 => fig[1, 1], 3 => fig[2, 1], 4 => fig[3, 1], 5 => fig[4, 1], 6 => fig[5, 1])
    for (final_time, fig_ref) in final_times_to_figure_map
        plot_data = scatter_plot_df(long_df, additive, metabolite, final_time)
        if final_time < 6
            plot_conc_vs_time_from_plot_data!(plot_data, fig_ref, false)
        else
            plot_conc_vs_time_from_plot_data!(plot_data, fig_ref, true)
        end
    end
    return fig
end

"""
    scatter_plot_df(long_df, additive, metabolite, final_time)

Returns data for the scatter layers for [`plot_regression`](@ref BloodStorageInSilico.AbsoluteQuant.plot_regression).

# Arguments
1. `long_df`: The long DataFrame from [`union_and_pivot_wide`](@ref BloodStorageInSilico.AbsoluteQuant.union_and_pivot_wide).
2. `additive`: Additive
3. `metabolite`: Metabolite id

# Returns
`NamedTuple`

NamedTuple with the following fields:

1. `scatter_df`: The DataFrame for the scatter points
2. `ylims`: The y limits for the scatter plot.
"""
function scatter_plot_df(long_df, additive, metabolite, final_time)
    scatter_df = @chain long_df begin
        @rsubset(
            :Additive == additive,
            :Metabolite == metabolite,
            :Time <= final_time,
            :Time >= final_time - 1
        )
        @select(:Time, :absolute_mM)
    end
    ylims_df = @chain long_df begin
        @rsubset(:Additive == additive, :Metabolite == metabolite)
        @combine(:ymin = minimum(:absolute_mM), :ymax = maximum(:absolute_mM))
    end
    ylims = (ylims_df[1, :ymin], ylims_df[1, :ymax])
    return (scatter_df = scatter_df, ylims = ylims)
end

"""
    plot_conc_vs_time_from_plot_data!(plot_data, fig_ref, time_label)

Used by [`plot_regression`](@ref BloodStorageInSilico.AbsoluteQuant.plot_regression). Draw the scatter and line layers onto a regression plot using the given data onto the specified section of the figure, optionally with time axis label. This function mutates the `Axis` object, so its return value isn't necessary.

# Arguments
1. `plot_data`: DataFrame with `:Time` and `:absolute_mM` columns to draw onto the plot
2. `fig_ref`: A Vector with the indices of the section of the figure to draw the plot onto (see calling function). Such as [1, 1], [2, 1], etc.
3. `time_label`: `true` will write the time axis label. `false` will suppress the time axis label.

# Returns
`GridLayout`

The `GridLayout` object upon which the drawing was made.
"""
function plot_conc_vs_time_from_plot_data!(plot_data, fig_ref, time_label)
    ax =
        time_label ? Axis(fig_ref, ylabel = "mM", xlabel = "Time (week)") :
        Axis(fig_ref, ylabel = "mM")
    plt =
        data(plot_data.scatter_df) *
        mapping(:Time, :absolute_mM) *
        (visual(Scatter) + linear(level = 0.95))
    draw!(ax, plt)
end

end

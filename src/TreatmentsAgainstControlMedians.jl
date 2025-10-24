module TreatmentsAgainstControlMedians

using CSV
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

export load_and_clean_2,
    c_means_metabolite_trajectories,
    plot_c_means_for_all_additives,
    plot_fuzzy_objectives_elbow,
    cluster_enrichment_analysis,
    load_gem_and_subsystems,
    plot_cluster_analysis

function load_and_clean_2()
    filename = joinpath("input", "Data Sheet 1.CSV")
    df1 = CSV.read(filename, DataFrame)
    df2 = stack(
        df1,
        Not([:Sample, :Time, :Additive]),
        variable_name = :MixedName,
        value_name = :Intensity,
    )
    df3 = subset(df2, :Additive => x -> x .== "01-Ctrl AS3")
    df4 = @combine(
        groupby(df3, [:MixedName, :Time]),
        :ControlMedianIntensity = median(skipmissing(:Intensity))
    )
    df5 = innerjoin(df2, df4, on = [:MixedName, :Time])
    df6 = transform(
        df5,
        [:Intensity, :ControlMedianIntensity] =>
            ByRow((x, y) -> x / y) => :ControlMedianNormalizedIntensity,
    )
    proportination_filename = joinpath("input", "Proportionation Sheet 2.csv")
    proportination_df = CSV.read(proportination_filename, DataFrame)
    df8 = innerjoin(df6, proportination_df, on = :MixedName)
    df9 = transform(
        df8,
        [:ControlMedianNormalizedIntensity, :Proportion] =>
            ByRow((x, y) -> x * y) => :SplitIntensity,
    )
    df10 = select(df9, [:Sample, :Time, :Additive, :Metabolite, :SplitIntensity])
    return df10
end

function load_gem_and_subsystems()
    gem_filename = joinpath("input", "RBC-GEM.json")
    model = load_model(gem_filename)
    reactions_rows = map(keys(model.reactions)) do r
        (
            RxnId = model.reactions[r]["id"],
            RxnName = model.reactions[r]["name"],
            Subsystem = model.reactions[r]["subsystem"],
        )
    end
    reactions_df1 = DataFrame(reactions_rows)
    metabolites_rows = []
    for r in keys(model.reactions)
        for m in keys(model.reactions[r]["metabolites"])
            row = (RxnId = model.reactions[r]["id"], Metabolite = m)
            push!(metabolites_rows, row)
        end
    end
    metabolites_df = DataFrame(metabolites_rows)
    subsystems_filename = joinpath("input", "Subsystem Category Map.csv")
    subsystems_df = CSV.read(subsystems_filename, DataFrame)
    reactions_df2 = leftjoin(reactions_df1, subsystems_df, on = :Subsystem => :name)
    return reactions_df2, metabolites_df
end

function prepare_everything_df_for_clustering(everything_df, additive)
    df0 = deepcopy(everything_df)
    df1 = subset(df0, :Additive => x -> x .== additive)
    df2 = select(df1, [:Metabolite, :Time, :SplitIntensity])
    df3 = unstack(df2, :Time, :SplitIntensity, combine = mean)
    df4 =
        filter(row -> all(!isnan, skipmissing([row[col] for col in names(df3)[2:7]])), df3)
    wide_timeseries_df = sort(df4, :Metabolite)
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

function c_means_metabolite_trajectories_in_additive(
    wide_timeseries_df,
    additive;
    n_clusters = 5,
    μ = 2.0,
)
    X = Matrix{Float64}(disallowmissing(wide_timeseries_df[:, Not(:Metabolite)]))
    println("Feature matrix: ", size(X, 1), " Metabolites, ", size(X, 2), " Time Points")
    result = fuzzy_cmeans(X', n_clusters, μ, maxiter = 200, display = :iter)
    weights_col_names = string.(axes(result.weights, 2))
    memberships_df = DataFrame(result.weights, weights_col_names)
    memberships_df.Metabolite = wide_timeseries_df.Metabolite
    memberships_df[!, :Additive] .= additive
    memberships_df[!, :NClusters] .= n_clusters
    fuzzy_objective = calc_fuzzy_objective(result, X, μ)
    return memberships_df, fuzzy_objective
end

function c_means_metabolite_trajectories(everything_df, max_clusters)
    # additives = unique(everything_df.Additive)
    additives_for_iterator =
        ["02-Adenosine", "01-Ctrl AS3", "03-Glutamine", "07-NAC", "08-Taurine"]
    c_means_long_dfs = []
    wide_timeseries_dfs = []
    fuzzy_objectives = []
    additives_rows = []
    n_clusters_rows = []
    for additive in additives_for_iterator
        for n_clusters in collect(2:max_clusters)
            wide_timeseries_df =
                prepare_everything_df_for_clustering(everything_df, additive)
            println("=" ^ 60)
            println(uppercase(additive), " ", n_clusters, " clusters ")
            c_means_df, fuzzy_objective = c_means_metabolite_trajectories_in_additive(
                wide_timeseries_df,
                additive,
                n_clusters = n_clusters,
            )
            wide_timeseries_df[!, :Additive] .= additive
            wide_timeseries_df[!, :NClusters] .= n_clusters
            c_means_long_df = stack(
                c_means_df,
                Not([:Additive, :Metabolite, :NClusters]),
                variable_name = :Cluster,
                value_name = :Weight,
            )
            println(first(c_means_long_df, 10))
            push!(fuzzy_objectives, fuzzy_objective)
            push!(additives_rows, additive)
            push!(n_clusters_rows, n_clusters)
            push!(wide_timeseries_dfs, wide_timeseries_df)
            push!(c_means_long_dfs, c_means_long_df)
        end
    end
    fuzzy_objectives_df = DataFrame(
        Additive = additives_rows,
        NClusters = n_clusters_rows,
        FuzzyObjective = fuzzy_objectives,
    )
    all_wide_timeseries_df = vcat(wide_timeseries_dfs...)
    all_c_means_df = vcat(c_means_long_dfs...)
    return all_c_means_df, all_wide_timeseries_df, fuzzy_objectives_df
end

function cluster_counts_for_additive(c_means_df, additive)
    println(additive)
    df1 = subset(c_means_df, :Additive => x -> x .== additive)
    df2 = DataFrames.combine(groupby(df1, :PrimaryCluster), nrow => :Count)
    return df2
end

function plot_c_means_for_additive(additive, c_means_df, wide_timeseries_df)
    println("Plotting c-means plot for $additive")
    c_means_df_copy = deepcopy(c_means_df)
    df0 = unstack(c_means_df_copy, :Cluster, :Weight)
    df_clusters = select(df0, Not([:Metabolite, :Additive, :NClusters]))
    df0.PrimaryCluster = [argmax(row) for row in eachrow(df_clusters)]
    df1 = select(df0, [:Metabolite, :Additive, :PrimaryCluster])
    df2 = innerjoin(df1, wide_timeseries_df, on = [:Additive, :Metabolite])
    df3 = stack(
        df2,
        Not([:Additive, :Metabolite, :PrimaryCluster, :NClusters]),
        variable_name = :Time,
        value_name = :MeanNormalizedIntensity,
    )
    df4 = sort(df3, [:Additive, :Metabolite, :PrimaryCluster, :Time])
    df5 = subset(df4, :Additive => x -> x .== additive)
    df5.Time = parse.(Int, df5.Time)
    plt_df = dropmissing(df5, :MeanNormalizedIntensity)
    time_points = unique(plt_df.Time)
    cluster_counts_df = cluster_counts_for_additive(df0, additive)
    println(first(cluster_counts_df, 5))
    cluster_counts_subtitle = join(
        [
            "Cluster $c, n=$n" for (c, n) in
            zip(cluster_counts_df[!, :PrimaryCluster], cluster_counts_df[!, :Count])
        ],
        ";\n",
    )
    plt =
        data(plt_df) *
        mapping(
            :Time,
            :MeanNormalizedIntensity,
            row = :PrimaryCluster,
            group = :Metabolite,
        ) *
        visual(Lines) *
        visual(alpha = 0.1)
    figure_options =
        (; size = (500, 1000), title = additive, subtitle = cluster_counts_subtitle)
    fig = draw(
        plt;
        figure = figure_options,
        axis = (; xticks = time_points),
        facet = (; linkxaxes = :minimal, linkyaxes = :minimal),
    )
    clean_additive = replace(additive, r"[^A-Za-z0-9]" => "_")
    fig_filename = joinpath("output", "c_means_plots", "$clean_additive.png")
    save(fig_filename, fig)
    println("Wrote $fig_filename")
end

function plot_c_means_for_all_additives(n_clusters, all_c_means_df, all_wide_timeseries_df)
    c_means_df = subset(all_c_means_df, :NClusters => x -> x .== n_clusters)
    additives = ["02-Adenosine", "01-Ctrl AS3", "03-Glutamine", "07-NAC", "08-Taurine"]
    for additive in additives
        println(">" ^ 60)
        println(uppercase(additive))
        wide_timeseries_df =
            subset(all_wide_timeseries_df, :Additive => x -> x .== additive)
        plot_c_means_for_additive(additive, c_means_df, wide_timeseries_df)
    end
end

function plot_fuzzy_objectives_elbow(fuzzy_objectives_df)
    xticks = unique(fuzzy_objectives_df.NClusters)
    plt =
        data(fuzzy_objectives_df) *
        mapping(:NClusters, :FuzzyObjective, row = :Additive) *
        visual(Lines)
    figure_options = (; size = (500, 1000), title = "C-Means Objective Elbow Plots")
    axis_options = (; xticks = xticks)
    facet_options = (; linkxaxes = :minimal, linkyaxes = :minimal)
    fig = draw(plt; figure = figure_options, axis = axis_options, facet = facet_options)
    fig_filename = joinpath("output", "c_means_plots", "elbows.png")
    save(fig_filename, fig)
    println("Wrote $fig_filename")
end

function cluster_enrichment_analysis(
    n_clusters,
    all_c_means_df,
    gem_reactions_df,
    gem_metabolites_df,
)
    all_c_means_df_copy = deepcopy(all_c_means_df)
    df0 = subset(all_c_means_df_copy, :NClusters => x -> x .== n_clusters)
    df05 = unstack(df0, :Cluster, :Weight)
    df_clusters = select(df05, Not([:Metabolite, :Additive, :NClusters]))
    df05.PrimaryCluster = [argmax(row) for row in eachrow(df_clusters)]
    df1 = select(df05, [:Metabolite, :Additive, :PrimaryCluster])
    df2 = innerjoin(df1, gem_metabolites_df, on = :Metabolite)
    df3 = innerjoin(df2, gem_reactions_df, on = :RxnId)
    gdf4 = @groupby(df3, [:Additive, :PrimaryCluster, :category])
    df5 = DataFrames.combine(gdf4, nrow => :Count)
    enrichment_df =
        sort(df5, [:Additive, :PrimaryCluster, :Count], rev = [false, false, true])
    metabolites_subsystems_df = sort(df3, [:Additive, :PrimaryCluster, :Metabolite])
    top3_df =
        DataFrames.combine(groupby(enrichment_df, [:Additive, :PrimaryCluster])) do sdf
            no_transport_df = subset(sdf, :category => x -> x .!= "Transport reactions")
            sort(no_transport_df, :Count, rev = true)[1:min(3, nrow(no_transport_df)), :]
        end
    return enrichment_df, metabolites_subsystems_df, top3_df
end

function plot_cluster_analysis(top3_df, additive)
    plt_df = subset(top3_df, :Additive => x -> x .== additive)
    println(plt_df)
    plt =
        data(plt_df) * mapping(:PrimaryCluster, :Count, color = :category) * visual(BarPlot)
    figure_options =
        (; size = (1000, 1000), title = "Top 3 Categories of Reactions in Each Cluster")
    fig = draw(plt; figure = figure_options)
    fig_filename = joinpath("output", "c_means_plots", "top3.png")
    save(fig_filename, fig)
    println("Wrote $fig_filename")
end

end

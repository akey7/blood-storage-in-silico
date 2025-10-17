module TreatmentsAgainstControlMedians

using CSV
using DataFrames
using DataFramesMeta
using Statistics
using AlgebraOfGraphics
using CairoMakie
using Makie
using GLM
using CategoricalArrays
using Base.Iterators
using ThreadsX
using MultipleTesting
using Clustering

export load_and_clean_2,
    plot_loess_for_all_metabolites,
    test_mixed_models,
    find_significant_metabolites_additives,
    c_means_metabolite_trajectories_in_additive

function load_and_clean_2()
    filename = joinpath("input", "Data Sheet 1.CSV")
    df1 = CSV.read(filename, DataFrame)
    df2 = stack(
        df1,
        Not([:Sample, :Time, :Additive]),
        variable_name = :Metabolite,
        value_name = :Intensity,
    )
    df3 = subset(df2, :Additive => x -> x .== "01-Ctrl AS3")
    df4 = @combine(
        groupby(df3, [:Metabolite, :Time]),
        :ControlMedianIntensity = median(skipmissing(:Intensity))
    )
    df5 = innerjoin(df2, df4, on = [:Metabolite, :Time])
    df6 = transform(
        df5,
        [:Intensity, :ControlMedianIntensity] =>
            ByRow((x, y) -> x / y) => :ControlMedianNormalizedIntensity,
    )
    df7 = select(
        df6,
        [:Sample, :Time, :Additive, :Metabolite, :ControlMedianNormalizedIntensity],
    )
    return df7
end

function plot_loess_for_metabolite(everything_df, metabolite)
    df = subset(everything_df, :Metabolite => x -> x .== metabolite)
    time_points = unique(df.Time)
    plt =
        data(df) *
        mapping(
            :Time => "Time",
            :ControlMedianNormalizedIntensity => "Control Median Normalized Intensity",
            color = :Additive => "Additive",
        ) *
        (visual(Scatter; markersize = 10, alpha = 0.3) + linear())
    fig = draw(
        plt;
        figure = (; size = (750, 500)),
        axis = (; title = metabolite, xticks = time_points),
    )
    return fig
end

function plot_loess_for_all_metabolites(df)
    metabolites = unique(df.Metabolite)
    for metabolite in metabolites
        fig = plot_loess_for_metabolite(df, metabolite)
        clean_metabolite = replace(metabolite, r"[^A-Za-z0-9]" => "_")
        filename = joinpath(
            "output",
            "control_median_normalized_intensity_plots",
            "$(clean_metabolite).png",
        )
        save(filename, fig)
        println("Wrote $filename")
    end
end

function find_significant_metabolites_additives(everything_df)
    control = "01-Ctrl AS3"
    fdr_threshold = 0.05
    frm = @formula(ControlMedianNormalizedIntensity ~ Time + AdditiveC)
    metabolites = unique(everything_df.Metabolite)
    additives = [
        additive for
        additive in unique(everything_df.Additive) if !contains(additive, control)
    ]
    df1 = deepcopy(everything_df)
    df1.AdditiveC = categorical(df1.Additive)
    pairs = vec(collect(product(additives, metabolites)))
    rows = ThreadsX.map(pairs) do pair
        additive, metabolite = pair
        println(additive, " ", metabolite)
        df2 = subset(
            df1,
            :Additive => x -> x .== additive .|| x .== control,
            :Metabolite => x -> x .== metabolite,
        )
        df3 = select(df2, [:ControlMedianNormalizedIntensity, :Time, :AdditiveC])
        model = lm(frm, df3)
        ct = coeftable(model)
        names = coefnames(model)
        pvals_vec = ct.cols[4]
        pvals = Dict(names .=> pvals_vec)
        p_value =
            isnan(pvals["AdditiveC: $additive"]) ? 1.0 : pvals["AdditiveC: $additive"]
        return (additive = additive, metabolite = metabolite, p_value = p_value)
    end
    result_df = DataFrame(rows)
    result_df.adj_p_value = adjust(result_df.p_value, BenjaminiHochberg())
    result_df.significant = result_df.adj_p_value .< fdr_threshold
    final_df = sort(result_df, :adj_p_value)
    return final_df
end

function c_means_metabolite_trajectories_in_additive(everything_df, additive)
    df1 = subset(everything_df, :Additive => x -> x .== additive)
    df2 = select(df1, [:Metabolite, :Time, :ControlMedianNormalizedIntensity])
    df3 = unstack(df2, :Time, :ControlMedianNormalizedIntensity, combine = mean)
    df4 =
        filter(row -> all(!isnan, skipmissing([row[col] for col in names(df3)[2:7]])), df3)
    df5 = sort(df4, :Metabolite)
    X = Matrix{Float64}(disallowmissing(df5[:, Not(:Metabolite)]))
    println(
        "Feature matrix: ",
        size(X, 1),
        " Metabolites, ",
        size(X, 2),
        " Time Points",
    )
    R = fuzzy_cmeans(X', 3, 2.0, maxiter = 200, display = :iter)
    M = R.centers
    memberships = R.weights
    println(memberships)
end

end

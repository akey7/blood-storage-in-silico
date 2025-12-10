module RawRelativeIntensities

using CSV
using DataFrames
using DataFramesMeta
using MultivariateStats
using StatsBase
using StatsModels
using Statistics

export load_relative_intensities, pca_relative_intensities

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
    wide_df = @chain long_df begin
        @rsubset(:Additive == additive)
        @rtransform(:Patient = :Sample[7:8])
        @select(:MixedName, :Patient, :Time, :Intensity)
        unstack([:Patient, :Time], :MixedName, :Intensity, combine = first)
        @orderby(:Patient, :Time)
    end
    patient_time_labels = @select(wide_df, :Patient, :Time)
    X = Matrix(select(wide_df, Not([:Patient, :Time])))
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
    M = fit(PCA, Xzt; maxoutdim = 6, mean = false)
    # display(M)
    scores = MultivariateStats.transform(M, Xzt)
    return (
        model = M,
        scores = scores,
        patient_time_labels = patient_time_labels,
        kept_columns = findall(good_cols),
        wide_df = wide_df,
    )
end

end

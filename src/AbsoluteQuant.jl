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

export load_absolute_quant, load_relative_quant

function load_absolute_quant()
    absolute_filename = joinpath("input", "Absolute Quant Data Sheet.xlsx")
    cells_day_1_df = DataFrame(XLSX.readtable(absolute_filename, "cells_day_1"))
    absolute_metabolite_ids = DataFrame(XLSX.readtable(absolute_filename, "metabolite_ids"))
    absolute_quant_df = @chain cells_day_1_df begin
        stack(Not(:id), variable_name = :mixed_name, value_name = :mmol_per_L)
        transform(:id => ByRow(x -> split(x, "_")[2]) => :sample_set)
        innerjoin(absolute_metabolite_ids, on = :mixed_name => :MixedName)
        transform([:mmol_per_L, :Proportion] => ByRow((x, y) -> x * y) => :prop_mmol_per_L)
        select([:sample_set, :id, :Metabolite, :prop_mmol_per_L])
    end
    normalization_df = @by absolute_quant_df [:sample_set, :Metabolite] begin
        :median_prop_mmol_per_L = median(:prop_mmol_per_L)
    end
    return absolute_quant_df, normalization_df
end

function load_relative_quant()
    relative_filename = joinpath("input", "Data Sheet 1.CSV")
    wide_df = CSV.read(relative_filename, DataFrame)
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
    relative_fold_change_df = @chain long_df begin
        innerjoin(ctrl_time_1_median_df, on = [:MixedName])
        @transform(@byrow :RelativeFoldChange = :Intensity / :CtrlTime1MedianIntensity)
    end
    return relative_fold_change_df
end

end

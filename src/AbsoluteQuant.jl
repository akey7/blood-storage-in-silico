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
        @transform(@byrow :sample_set = split(:id, "_")[2])
        innerjoin(absolute_metabolite_ids, on = :mixed_name => :MixedName)
        @transform(@byrow :prop_mmol_per_L = :mmol_per_L * :Proportion)
        select([:sample_set, :id, :Metabolite, :prop_mmol_per_L])
    end
    normalization_df = @by absolute_quant_df :Metabolite begin
        :median_prop_mmol_per_L = median(:prop_mmol_per_L)
    end
    return absolute_quant_df, normalization_df
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
        @transform(@byrow :FoldChange = :Intensity / :CtrlTime1MedianIntensity)
        innerjoin(proportination_df, on = :MixedName)
        @select(:Sample, :Time, :Additive, :Metabolite, :FoldChange)
    end
    return fold_changes_df
end

function combine_relative_and_absolute_quant(relative_fold_changes_df)
end

end

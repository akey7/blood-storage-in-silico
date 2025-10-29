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

export load_and_clean_3

function load_and_clean_3()
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
    normalization_df = @chain absolute_quant_df begin
        @groupby(:sample_set, :Metabolite)
        @combine(:median_prop_mmol_per_L = median(:prop_mmol_per_L))
    end
    return absolute_quant_df, normalization_df
end

end

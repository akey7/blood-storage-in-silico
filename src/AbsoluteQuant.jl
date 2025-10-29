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
    cells_long_df1 = stack(
        cells_day_1_df,
        Not(:id),
        variable_name = :metabolite_name,
        value_name = :mmol_per_L,
    )
    cells_long_df2 =
        transform(cells_long_df1, :id => ByRow(x -> split(x, "_")[2]) => :sample_set)
    println(first(cells_long_df2, 10))
    println(first(absolute_metabolite_ids, 10))
end

end

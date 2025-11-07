using CSV
using DataFrames

include("src/MasspyInterface.jl")
using .MasspyInterface

println(">" ^ 10, " STEADY-STATE SIMULATION ", "<" ^ 10)
ss_conc_filename = joinpath("input", "masspy_interface", "conc_sol_ss.csv")
ss_conc_df = CSV.read(ss_conc_filename, DataFrame)
ss_flux_filename = joinpath("input", "masspy_interface", "flux_sol_ss.csv")
ss_flux_df = CSV.read(ss_flux_filename, DataFrame)
plot_solutions("Steady State", ss_conc_df, ss_flux_df)

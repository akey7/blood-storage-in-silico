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

println(">" ^ 10, " PERTURBED SIMULATION ", "<" ^ 10)
perturbed_conc_filename = joinpath("input", "masspy_interface", "conc_sol_perturbed.csv")
perturbed_conc_df = CSV.read(ss_conc_filename, DataFrame)
perturbed_flux_filename = joinpath("input", "masspy_interface", "flux_sol_perturbed.csv")
perturbed_flux_df = CSV.read(ss_flux_filename, DataFrame)
plot_solutions("Perturbed", perturbed_conc_df, perturbed_flux_df)

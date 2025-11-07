module MasspyInterface

using DataFrames
using DataFramesMeta
using CairoMakie
using AlgebraOfGraphics

export plot_solutions

function plot_solutions(base_title, conc_df, flux_df)
    plt_conc_df = stack(conc_df, Not(:t), variable_name = :metabolite, value_name = :conc)
    plt_flux_df = stack(flux_df, Not(:t), variable_name = :flux, value_name = :rate)
    figure_options = (; size = (750, 500))
    plt_conc = data(plt_conc_df) * mapping(:t, :conc, color = :metabolite) * visual(Lines)
    fig_conc = draw(
        plt_conc;
        figure = figure_options,
        axis = (; title = "$base_title Concentrations"),
    )
    fig_conc_filename = joinpath("output", "masspy_model", "$base_title concentrations.png")
    save(fig_conc_filename, fig_conc)
    println("Wrote $fig_conc_filename")
    plt_flux = data(plt_flux_df) * mapping(:t, :rate, color = :flux) * visual(Lines)
    fig_flux =
        draw(plt_flux; figure = figure_options, axis = (; title = "$base_title Fluxes"))
    fig_flux_filename = joinpath("output", "masspy_model", "$base_title fluxes.png")
    save(fig_flux_filename, fig_flux)
    println("Wrote $fig_flux_filename")
end

end

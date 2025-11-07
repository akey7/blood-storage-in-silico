module MasspyInterface

using DataFrames
using DataFramesMeta
using CairoMakie
using AlgebraOfGraphics

export plot_solutions

function plot_solutions(base_title, conc_df, flux_df)
    plt_conc_df = @chain conc_df begin
        stack(Not(:t), variable_name = :metabolite, value_name = :conc)
        @rsubset(!isapprox(:t, 0.0))
    end
    plt_flux_df = @chain flux_df begin
        stack(Not(:t), variable_name = :flux, value_name = :rate)
        @rsubset(!isapprox(:t, 0.0))
    end
    figure_options = (; size = (750, 500))
    plt_conc = data(plt_conc_df) * mapping(:t, :conc, color = :metabolite) * visual(Lines)
    fig_conc = draw(
        plt_conc;
        figure = figure_options,
        axis = (;
            xscale = log10,
            yscale = log10,
            xlabel = "log10(Time)",
            ylabel = "log10(Concentration)",
            title = "$base_title Concentrations",
        ),
    )
    fig_conc_filename = joinpath("output", "masspy_interface", "$base_title concentrations.png")
    save(fig_conc_filename, fig_conc)
    println("Wrote $fig_conc_filename")
    plt_flux = data(plt_flux_df) * mapping(:t, :rate, color = :flux) * visual(Lines)
    fig_flux = draw(
        plt_flux;
        figure = figure_options,
        axis = (;
            xscale = log10,
            xlabel = "log10(Time)",
            ylabel = "Flux",
            title = "$base_title Fluxes",
        ),
    )
    fig_flux_filename = joinpath("output", "masspy_interface", "$base_title fluxes.png")
    save(fig_flux_filename, fig_flux)
    println("Wrote $fig_flux_filename")
end

end

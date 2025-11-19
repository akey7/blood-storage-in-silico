module DynamicModel

using Catalyst
using Catalyst: species, parameters, reactions, reactionrates
using ModelingToolkit
using DifferentialEquations
using CairoMakie
using GraphMakie
using NetworkLayout
using Latexify

export run_glycolysis, plot_glycolysis

# k1(; k2, q10, t2 = 27.0, t1 = 4.0) = k2 / q10^((t2-t1)/10)

k1(; k2, q10, t2 = 27.0, t1 = 4.0) = 1.0 * k2

function run_glycolysis()
    # Q10 values from Yurkovich et al, 2017 Table 1 and Figure 3
    #
    # Reaction directionality from Yurkovich et al, 2017 Figure 3
    #
    # k2 values from PERC values in Ch. 10 of Systems Biology:
    # Simulation of Dynamic Network States by Palsson.
    # https://masspy.readthedocs.io/en/latest/education/sb2/chapters/sb2_chapter10.html

    @info "Creating reaction network..."
    rn = @reaction_network glycolysis begin
        @require_declaration

        @parameters begin
            k_hex1_f
            k_pgi_f
            k_pgi_r
            k_pfk_f
            k_fba_f
            k_tpi_f
            k_tpi_r
            k_gapd_f
            k_gapd_r
            k_pgk_f
            k_pgk_r
            k_pgm_f
            k_pgm_r
            k_eno_f
            k_eno_r
            k_pyk_f
            k_ldh_f
            k_ldh_r
            k_SK_glc__D_c_f
            k_SK_lac__L_f
            k_SK_amp_c
            k_adk_f
            k_atpm_f
            k_DM_amp_c_f
            k_DM_nadh_c_f
            k_SK_pyr_c_f
            k_SK_pyr_c_r
            k_SK_h_c_f
            k_SK_h_c_r
            k_SK_h2o_c_f
            k_SK_h2o_c_r
        end

        @species begin
            glc__D_c(t)
            g6p_c(t)
            f6p_c(t)
            fdp_c(t)
            dhap_c(t)
            g3p_c(t)
            _13dpg_c(t)
            _3pg_c(t)
            _2pg_c(t)
            pep_c(t)
            pyr_c(t)
            lac__L_c(t)
            nad_c(t)
            nadh_c(t)
            amp_c(t)
            adp_c(t)
            atp_c(t)
            pi_c(t)
            h_c(t)
            h2o_c(t)
            a_tot(t)
        end

        @observables begin
            a_tot ~ amp_c + adp_c + atp_c
        end

        # Boundary reactions
        k_DM_amp_c_f, amp_c --> 0
        (k_SK_pyr_c_f, k_SK_pyr_c_r), pyr_c <--> 0
        k_SK_lac__L_f, lac__L_c --> 0
        k_SK_glc__D_c_f, 0 --> glc__D_c
        k_SK_amp_c, 0 --> amp_c
        (k_SK_h_c_f, k_SK_h_c_r), h_c <--> 0
        (k_SK_h2o_c_f, k_SK_h2o_c_r), h2o_c <--> 0

        k_DM_nadh_c_f, nadh_c --> h_c + nad_c
        k_hex1_f, atp_c + glc__D_c --> adp_c + g6p_c + h_c
        (k_pgi_f, k_pgi_r), g6p_c <--> f6p_c
        k_pfk_f, atp_c + f6p_c --> adp_c + fdp_c + h_c
        k_fba_f, fdp_c --> dhap_c + g3p_c
        (k_tpi_f, k_tpi_r), dhap_c <--> g3p_c
        (k_gapd_f, k_gapd_r), g3p_c + nad_c + pi_c <--> _13dpg_c + h_c + nadh_c
        (k_pgk_f, k_pgk_r), _13dpg_c + adp_c <--> _3pg_c + atp_c
        (k_pgm_f, k_pgm_r), _3pg_c <--> _2pg_c
        (k_eno_f, k_eno_r), _2pg_c <--> h2o_c + pep_c
        k_pyk_f, adp_c + h_c + pep_c --> atp_c + pyr_c
        (k_ldh_f, k_ldh_r), h_c + nadh_c + pyr_c <--> lac__L_c + nad_c
        k_adk_f, 2*adp_c --> amp_c + atp_c
        k_atpm_f, atp_c + h2o_c --> adp_c + h_c + pi_c
    end

    println("Reaction network name: ", nameof(rn))
    println("Reaction network parameters: ", parameters(rn))
    println("Reaction network species: ", species(rn))

    ps = [
        :k_SK_glc__D_c_f => 1.12,
        :k_hex1_f => k1(k2 = 0.7, q10 = 2.60),
        :k_pgi_f => k1(k2 = 3644.444, q10 = 2.72),
        :k_pgi_r => k1(k2 = 3644.444, q10 = 2.72),
        :k_pfk_f => k1(k2 = 35.369, q10 = 2.65),
        :k_fba_f => k1(k2 = 2834.568, q10 = 2.65),
        :k_tpi_f => k1(k2 = 34.356, q10 = 2.65),
        :k_tpi_r => k1(k2 = 34.356, q10 = 2.65),
        :k_gapd_f => k1(k2 = 3376.749, q10 = 2.63),
        :k_gapd_r => k1(k2 = 3376.749, q10 = 2.63),
        :k_pgk_f => k1(k2 = 1273531.270, q10 = 2.53),
        :k_pgk_r => k1(k2 = 1273531.270, q10 = 2.53),
        :k_pgm_f => k1(k2 = 4868.589, q10 = 2.57),
        :k_pgm_r => k1(k2 = 4868.589, q10 = 2.57),
        :k_eno_f => k1(k2 = 1763.741, q10 = 2.57),
        :k_eno_r => k1(k2 = 1763.741, q10 = 2.57),
        :k_pyk_f => k1(k2 = 454.386, q10 = 2.59),
        :k_ldh_f => k1(k2 = 1112.574, q10 = 2.61),
        :k_ldh_r => k1(k2 = 1112.574, q10 = 2.61),
        :k_atpm_f => 1.400,
        :k_SK_lac__L_f => 10.0,
        :k_adk_f => 1.0e6,
        :k_SK_amp_c => 0.014,
        :k_DM_amp_c_f => 0.161,
        :k_DM_nadh_c_f => 7.442,
        :k_SK_pyr_c_f => 744.186,
        :k_SK_pyr_c_r => 744.186,
        :k_SK_h_c_f => 1.0e6,
        :k_SK_h_c_r => 1.0e6,
        :k_SK_h2o_c_f => 1.0e6,
        :k_SK_h2o_c_r => 1.0e6,
    ]

    u0 = [
        :glc__D_c => 1.0,
        :g6p_c => 0.049,
        :f6p_c => 0.02,
        :fdp_c => 0.015,
        :dhap_c => 0.16,
        :g3p_c => 0.007,
        :_13dpg_c => 0.0,
        :_3pg_c => 0.077,
        :_2pg_c => 0.011,
        :pep_c => 0.017,
        :pyr_c => 0.06,
        :lac__L_c => 1.36,
        :nad_c => 0.059,
        :nadh_c => 0.03,
        :amp_c => 0.087,
        :adp_c => 0.29,
        :atp_c => 1.6,
        :pi_c => 2.5,
        :h_c => 0.0,
        :h2o_c => 1.0,
    ]

    tspan = (0.0, 1.0)
    @info "Formulating glycolysis ODEProblem..."
    prob = ODEProblem(rn, u0, tspan, ps)
    @info "Solving ODEs..."
    sol = solve(prob, Rodas5P(); reltol = 1.0e-8, abstol = 1.0e-10)

    @info "Plotting main metabolites..."
    species2 = [
        :glc__D_c,
        :g6p_c,
        :f6p_c,
        :fdp_c,
        :dhap_c,
        :g3p_c,
        :_13dpg_c,
        :_3pg_c,
        :_2pg_c,
        :pep_c,
        :pyr_c,
        :lac__L_c,
    ]
    title2 = "Glycolysis Main Metabolites"
    plot_metabolites(sol, species2, title2)

    @info "Plotting cofactors..."
    title3 = "Glycolysis NAD and NADH"
    species3 = [:nad_c, :nadh_c]
    plot_metabolites(sol, species3, title3)

    @info "Plotting ATP/ADP..."
    species4 = [:pi_c, :amp_c, :adp_c, :atp_c, :a_tot]
    plot_metabolites(sol, species4, "Glycolysis AMP ADP ATP")

    @info "Plotting reaction graph..."
    plot_reaction_network_graph(rn)

    @info "Calling Latexify..."
    copy_to_clipboard(true)
    latexify(rn; form = :ode)
end

function plot_metabolites(sol, species, title)
    size = (900, 600)
    fig = Figure(; size = size)
    ax = Axis(fig[1, 1]; xlabel = "Time", ylabel = "Concentration", title = title)
    for sp in species
        lines!(ax, sol.t, sol[sp, :]; label = string(sp), linewidth = 2)
    end
    axislegend(ax; position = :rb, framevisible = false)
    fig_filename = joinpath("output", "kinetic_model", "$(title).png")
    save(fig_filename, fig)
end

function plot_reaction_network_graph(rn)
    g = plot_network(rn)
    g_filename = joinpath("output", "kinetic_model", "Glycolysis Network.png")
    save(g_filename, g)
end

end

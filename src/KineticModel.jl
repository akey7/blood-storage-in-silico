module KineticModel

using Catalyst
using Catalyst: species, parameters, reactions, reactionrates
using ModelingToolkit
using DifferentialEquations
using CairoMakie

export glycolysis

function glycolysis()
    @parameters k_hex1_f, k_hex1_r, k_pgi_f, k_pgi_r, k_pfk_f, k_pfk_r
    @parameters k_fba_f, k_fba_r, k_tpi_f, k_tpi_r, k_gapd_f, k_gapd_r
    @parameters k_pgk_f, k_pgk_r, k_pgm_f, k_pgm_r, k_eno_f, k_eno_r
    @parameters k_pyk_f, k_pyk_r, k_ldh_f, k_ldh_r
    @variables t
    @species glc__D_c(t) g6p_c(t) f6p_c(t) fdp_c(t) dhap_c(t)
    @species g3p_c(t) _13dpg_c(t) _3pg_c(t) _2pg_c(t) pep_c(t)
    @species pep_c(t) pyr_c(t) lac__L_c(t) nad_c(t) nadh_c(t)
    @species amp_c(t) adp_c(t) atp_c(t) pi_c(t) h_c(t)
    @species h2o_c(t)

    glycolysis = @reaction_network begin
        (k_hex1_f, k_hex1_r), atp_c + glc__D_c <--> adp_c + g6p_c + h_c
        (k_pgi_f, k_pgi_r), g6p_c <--> f6p_c
        (k_pfk_f, k_pfk_r), atp_c + f6p_c <--> adp_c + fdp_c + h_c
        (k_fba_f, k_fba_r), fdp_c <--> dhap_c + g3p_c
        (k_tpi_f, k_tpi_r), dhap_c <--> g3p_c
        (k_gapd_f, k_gapd_r), g3p_c + nad_c + pi_c <--> _13dpg_c + h_c + nadh_c
        (k_pgk_f, k_pgk_r), _13dpg_c + adp_c <--> _3pg_c + atp_c
        (k_pgm_f, k_pgm_r), _3pg_c <--> _2pg_c
        (k_eno_f, k_eno_r), _2pg_c <--> h2o_c + pep_c
        (k_pyk_f, k_pyk_r), adp_c + h_c + pep_c <--> atp_c + pyr_c
        (k_ldh_f, k_ldh_r), h_c + nadh_c + pyr_c <--> lac__L_c + nad_c
    end

    p = [
        k_hex1_f => 0.7,
        k_hex1_r => 0.0,
        k_pgi_f => 3644.444,
        k_pgi_r => 0.0,
        k_pfk_f => 35.369,
        k_pfk_r => 0.0,
        k_fba_f => 2834.568,
        k_fba_r => 0.0,
        k_tpi_f => 34.356,
        k_tpi_r => 0.0,
        k_gapd_f => 3376.749,
        k_gapd_r => 0.0,
        k_pgk_f => 1273531.270,
        k_pgk_r => 0.0,
        k_pgm_f => 4868.589,
        k_pgm_r => 0.0,
        k_eno_f => 1763.741,
        k_eno_r => 0.0,
        k_pyk_f => 454.386,
        k_pyk_r => 0.0,
        k_ldh_f => 1112.574,
        k_ldh_r => 0.0,
    ]

    u0 = [
        glc__D_c => 1.0,
        g6p_c => 0.0486,
        f6p_c => 0.0198,
        fdp_c => 0.0146,
        dhap_c => 0.16,
        g3p_c => 0.00728,
        _13dpg_c => 0.000243,
        _3pg_c => 0.0773,
        _2pg_c => 0.0113,
        pep_c => 0.017,
        pyr_c => 0.060301,
        lac__L_c => 1.36,
        nad_c => 0.0589,
        nadh_c => 0.0301,
        amp_c => 0.0867281,
        adp_c => 0.29,
        atp_c => 1.6,
        pi_c => 2.5,
        h_c => 8.99757e-05,
        h2o_c => 1.0,
    ]

    tspan = (0.0, 1.0)
    @info "Formulating glycolysis ODEProblem..."
    prob = ODEProblem(glycolysis, u0, tspan, p)
    @info "Solving ODEs..."
    sol = solve(prob, Rodas5(); reltol = 1.0e-8, abstol = 1.0e-10)

    @info "Plotting solution..."
    fig1 = Figure(resolution = (900, 600))
    ax1 = Axis(
        fig1[1, 1];
        xlabel = "time (s)",
        ylabel = "concentration (mM)",
        title = "Glucose In, Lactate Out",
    )

    species1 = [glc__D_c, lac__L_c]
    labels = ["glc__D_c", "lac__L_c"]
    for (sp, label) in zip(species1, labels)
        lines!(ax1, sol.t, sol[sp, :]; label = label, linewidth = 2)
    end
    axislegend(ax1; position = :rb, framevisible = false)
    fig1_filename = joinpath("output", "glucose_in_lactate_out.png")
    save(fig1_filename, fig1)
end

end

module KineticModel 

using Catalyst
using Catalyst: species, parameters, reactions, reactionrates
using ModelingToolkit
using DifferentialEquations

export glycolysis

function glycolysis()
    @parameters k_hex1_f, k_hex1_r, k_pgi_f, k_pgi_r, k_pfk_f, k_pfk_r
    @parameters k_fba_f, k_fba_r, k_tpi_f, k_tpi_r, k_gapd_f, k_gapd_r
    @parameters k_pgk_f, k_pgk_r, k_pgm_f, k_pgm_r, k_eno_f, k_eno_r
    @parameters k_eno_f, k_eno_r, k_pyk_f, k_pyk_r, k_ldh_f, k_ldh_r
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
        (k_ldh_f, k_ldh_r), h_c + nadh_c + pyr_c <=> lac__L_c + nad_c
    end
end

end

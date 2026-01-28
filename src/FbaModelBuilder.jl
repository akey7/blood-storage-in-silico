module FbaModelBuilder

using COBREXA
import SBMLFBCModels as S
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel: Model, Reaction, Metabolite, Gene, Coupling

export load_base_rbc_gem, create_fba_model

"""
    load_base_rbc_gem()

Load the base RBC-GEM from which the model for uFBA sampling will be made

# Returns
`A.CanonicalModel.Model`

New model with the entire RBC-GEM.
"""
function load_base_rbc_gem()
    println("> Loading RBC-GEM")
    rbc_gem_path = joinpath("input", "RBC-GEM.xml")
    rbc_gem = load_model(S.SBMLFBCModel, rbc_gem_path, A.CanonicalModel.Model)

    # println("Metabolites: $(length(rbc_gem.metabolites))")
    # println("Reactions: $(length(rbc_gem.reactions))")

    return rbc_gem
end

"""
    create_fba_model(base_gem::Union{A.CanonicalModel.Model,Nothing}; add_exchanges::Bool = true)

Creates the three pathway (glycolysis, pentose phosphate, purine salvage) model
for the uFBA study.

# Arguments
1. `base_gem::Union{A.CanonicalModel.Model,Nothing}`: The base gem loaded by `load_base_rbc_gem`. If left as `nothing`, this function will call `load_base_rbc_gem` directly.

2. `add_exchanges::Bool = true`: If `true`, this function will add exchanges to the model. `false` will skip adding exchanges.

# Returns
`A.CanonicalModel.Model`

Returns the newly constructed three pathway model.
"""
function create_fba_model(
    base_gem::Union{A.CanonicalModel.Model,Nothing};
    add_exchanges::Bool = true,
)
    if add_exchanges
        @info "Building 3P model and adding exchanges"
    else
        @info "Building 3P model without exchanges"
    end

    rbc_gem = isnothing(base_gem) ? load_base_rbc_gem() : deepcopy(base_gem)

    println("> Glycolysis")

    glycolysis_reaction_ids = [
        "R_HEX1",
        "R_PGI",
        "R_PFK",
        "R_FBA",
        "R_TPI",
        "R_GAPD",
        "R_PGK",
        "R_PGM",
        "R_ENO",
        "R_PYK",
        "R_LDH_L",
    ]

    # println(glycolysis_reaction_ids)

    println("> RL Shunt")

    rl_shunt_reaction_ids = ["R_DPGM", "R_DPGase"]
    # println(rl_shunt_reaction_ids)

    println("> Pentose phosphate pathway")

    ppp_reaction_ids =
        ["R_G6PDH2", "R_PGL", "R_GND", "R_RPI", "R_RPE", "R_TKT1", "R_TALA", "R_TKT2"]

    # println(ppp_reaction_ids)

    println("> Purine metabolism")

    purine_metabolism_reaction_ids = [
        "R_PRPPS",
        "R_PPM",
        "R_HXPRT",
        "R_ADPT",
        "R_PUNP5",
        "R_NTD11",
        "R_AMPDA",
        "R_NTD7",
        "R_ADA",
        "R_ADNK1",
        "R_ADK1",
        "R_PPA",
    ]

    # println(purine_metabolism_reaction_ids)

    println("> Transporters")

    transporter_reactions_ids = [
        "R_GLC_Dt",
        "R_PYRt2",
        "R_L_LACt2",
        "R_HYXNt",
        "R_INSt",
        "R_ADEt",
        "R_ADNt",
        "R_CO2t",
        "R_NH4c",
        "R_NH4e",
        "R_NH3t",
        "R_PIt",
        "R_Ht",
        "R_H2Ot",
    ]

    # println(transporter_reactions_ids)

    if add_exchanges
        println("> Exchanges")

        exchange_reactions_ids = [
            "R_EX_glc__D_e",
            "R_EX_pyr_e",
            "R_EX_lac__L_e",
            "R_EX_hxan_e",
            "R_EX_ins_e",
            "R_EX_ade_e",
            "R_EX_adn_e",
            "R_EX_co2_e",
            "R_EX_pi_e",
            "R_EX_nh4_e",
            "R_EX_nh3_e",
            "R_EX_h_e",
            "R_EX_h2o_e",
        ]

        # println(exchange_reactions_ids)
    else
        println("> Skipping exchanges")
    end

    println("> Collecting reactions and discovering metabolites")

    if add_exchanges
        all_reaction_ids = [
            glycolysis_reaction_ids
            rl_shunt_reaction_ids
            ppp_reaction_ids
            purine_metabolism_reaction_ids
            transporter_reactions_ids
            exchange_reactions_ids
        ]
    else
        all_reaction_ids = [
            glycolysis_reaction_ids
            rl_shunt_reaction_ids
            ppp_reaction_ids
            purine_metabolism_reaction_ids
            transporter_reactions_ids
        ]
    end

    discovered_metabolite_ids::Vector{String} = []
    for reaction_id ∈ all_reaction_ids
        for metabolite_id ∈ keys(rbc_gem.reactions[reaction_id].stoichiometry)
            push!(discovered_metabolite_ids, metabolite_id)
        end
    end

    println("Discovered $(length(discovered_metabolite_ids)) metabolites.")

    model = Model()

    println(model)

    for discovered_metabolite_id ∈ discovered_metabolite_ids
        model.metabolites[discovered_metabolite_id] =
            deepcopy(rbc_gem.metabolites[discovered_metabolite_id])
        # println("$discovered_metabolite_id")
    end

    println("> Adding reactions and exchanges to model")

    for reaction_id ∈ all_reaction_ids
        model.reactions[reaction_id] = deepcopy(rbc_gem.reactions[reaction_id])
        lower_bound = model.reactions[reaction_id].lower_bound
        upper_bound = model.reactions[reaction_id].upper_bound
        println(reaction_id, " (", lower_bound, ", ", upper_bound, ")")
    end

    println("> Add ATP load")

    model.reactions["R_LOAD_ATP"] = Reaction(
        name = "LOAD_ATP",
        stoichiometry = Dict(
            "M_atp_c" => -1.0,
            "M_h2o_c" => -1.0,
            "M_adp_c" => 1.0,
            "M_pi_c" => 1.0,
            "M_h_c" => 1.0,
        ),
        objective_coefficient = 1.0,
        lower_bound = 0.0,
        upper_bound = 1.0,
    )

    println(model.reactions["R_LOAD_ATP"])

    println("> Adding NADH load")

    # Load due to methemoglobin reduction via CytB5
    model.reactions["R_LOAD_NADH"] = Reaction(
        name = "LOAD_NADH",
        stoichiometry = Dict("M_nadh_c" => -1.0, "M_h_c" => 1.0, "M_nad_c" => 1.0),
        objective_coefficient = 1.0,
        lower_bound = 0.0,
        upper_bound = 1.0,
    )

    println(model.reactions["R_LOAD_NADH"])

    println("> Adding NADPH load")

    # Load due to glutathione reduction from GSSG to GSH
    model.reactions["R_LOAD_NADPH"] = Reaction(
        name = "LOAD_NADPH",
        stoichiometry = Dict("M_nadph_c" => -1.0, "M_h_c" => 1.0, "M_nadp_c" => 1.0),
        objective_coefficient = 1.0,
        lower_bound = 0.0,
        upper_bound = 1.0,
    )

    println(model.reactions["R_LOAD_NADPH"])

    return model
end

end

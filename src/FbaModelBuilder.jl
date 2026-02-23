module FbaModelBuilder

using COBREXA
using DataFrames
using DataFramesMeta
using Accessors
import SBMLFBCModels as S
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel: Model, Reaction, Metabolite, Gene, Coupling

export load_base_rbc_gem, create_fba_model, default_exchanges, as3_exchanges

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

function default_exchanges()
    exchange_reaction_ids = [
        # Original exchanges
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

        # Exchanges added for AS-3
        "R_EX_pi_e",
        "R_EX_cit_e",
        "R_EX_na1_e",
        "R_EX_cl_e",
    ]
    return exchange_reaction_ids
end

function as3_exchanges()
    exchange_reaction_ids = [
        "R_EX_pi_e",
        "R_EX_cit_e",
        "R_EX_na1_e",
        "R_EX_cl_e",
        "R_EX_ade_e",
        "R_EX_glc__D_e",
    ]
    return exchange_reaction_ids
end

"""
   find_flux_bounds_overrides(flux_bounds_overrides_df::Union{Nothing,DataFrame}, reaction_id::String) 

Find flux bounds to override RBC-GEM bounds if such an override exists. If an override exists, return a `Tuple{Float64,Float64}` with reaction bounds override. If the flux overrides DataFrame is `nothing` or the reaction id does not exist in the provided DataFrame, this function returns `nothing`.

# Arguments
1. `flux_bounds_overrides_df::Union{Nothing,DataFrame}`: DataFrame that contains the flux bounds overrides. If `nothing`, then this function will simply return `nothing`.
2. `reaction_id::String`: Reaction id to search for an override.

# Returns
`Union{Nothing,Tuple{Float64,Float64}}`

Returns either `nothing` or flux bounds as described above.
"""
function find_flux_bounds_overrides(
    flux_bounds_overrides_df::Union{Nothing,DataFrame},
    reaction_id::String,
)
    if isnothing(flux_bounds_overrides_df)
        return nothing
    else
        override_df = @rsubset(flux_bounds_overrides_df, :reaction_id == reaction_id)
        if nrow(override_df) < 1
            return nothing
        else
            lower_bound = override_df[1, :lb]
            upper_bound = override_df[1, :ub]
            return (lower_bound, upper_bound)
        end
    end
end

"""
    create_fba_model(base_gem::Union{A.CanonicalModel.Model,Nothing}; exchanges::Union{Nothing,Vector{String}} = nothing, flux_bounds_overrides_df::Union{Nothing,DataFrame} = nothing)

Creates the three pathway (glycolysis, pentose phosphate, purine salvage) model
for the uFBA study.

# Arguments
1. `base_gem::Union{A.CanonicalModel.Model,Nothing}`: The base gem loaded by `load_base_rbc_gem`. If left as `nothing`, this function will call `load_base_rbc_gem` directly.
2. `exchanges::Union{Nothing,Vector{String}} = nothing`: If `nothing`, this function will not add exchanges to the model. If specified, the listed exchanges are added.
3. `flux_bounds_overrides_df::Union{Nothing,DataFrame} = nothing`: If specified, this DataFrame contains flux bounds for reactions that will override the RBC-GEM's flux bounds.

# Returns
`A.CanonicalModel.Model`

Returns the newly constructed three pathway model.
"""
function create_fba_model(
    base_gem::Union{A.CanonicalModel.Model,Nothing};
    exchanges::Union{Nothing,Vector{String}} = nothing,
    flux_bounds_overrides_df::Union{Nothing,DataFrame} = nothing,
)
    if !isnothing(exchanges)
        @info "Building FBA model and adding exchanges"
    else
        @info "Building FBA model without exchanges"
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
        "R_PEPCK",
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
        "R_NTDIMP",
        "R_AMPDA",
        "R_NTDAMP",
        "R_ADA",
        "R_ADNK1",
        "R_ADK1",
        "R_PPA",
        "R_PUNP3",
        # "R_XAO2",  # Zero flux
        # "R_XAO",  # Zero flux
    ]

    # println(purine_metabolism_reaction_ids)

    println("> Methionine Salvage and Metabolism")
    met_salvage_reaction_ids =
        ["R_UNK3", "R_AHC", "R_MDRPD", "R_METAT", "R_MTRI", "R_ARDFE2"]

    println("> Citric Acid Cycle")
    # citric_reaction_ids = ["R_ACITL", "R_FUM", "R_MDH"]  # All citric reactions have zero flux
    citric_reaction_ids = []

    println("> Arginine and Proline Metabolism")
    arg_pro_reaction_ids = ["R_ADMDC", "R_MTAP"]

    println("> Nucleotide Metabolism")
    nucleotide_reaction_ids = [
        # "R_ADNCYC",  # Broken reaction
        "R_GMPR",
        "R_GMPS2",
        "R_GUACYC",
        "R_IMPD",
        "R_NDPK1",
        "R_NDPK2",
        "R_NTDGMP",
        "R_PDEG",
        "R_UMPK1",
        "R_GK1",
    ]

    println("> Glutamate Metabolism")
    glutamate_reaction_ids = ["R_ALATA_L", "R_GLNS", "R_GLUN"]

    println("> Glutathione Metabolism")
    # glutathione_reaction_ids =
    #     ["R_AMPTASECG", "R_GLUCYS", "R_GTHP", "R_GTHS", "R_GTHOy", "R_GGLUCTC"]  # GGLUCTC has zero flux
    glutathione_reaction_ids = ["R_AMPTASECG", "R_GLUCYS", "R_GTHP", "R_GTHS", "R_GTHOy"]

    println("> Urea cycle/amino group metabolism")
    urea_reaction_ids = ["R_ARGN", "R_ORNDC", "R_SPMS", "R_SPRMS"]

    println("> Glycine, Serine, and Threonine Metabolism")
    # glycine_serine_threonine_reaction_ids = ["R_GHMT2"]  # GHMT2 zero flux
    glycine_serine_threonine_reaction_ids = []

    println("> Folate Metabolism")
    folate_reaction_ids = ["R_FTHFL", "R_MTHFC", "R_MTHFD"]

    println("> Fructose and Mannose Metabolism")
    fructose_mannose_reaction_ids = ["R_HEX4", "R_HEX7", "R_MAN6PI", "R_SBTD_D2", "R_SBTRa"]

    println("> Pyrimidine Catabolism")
    pyrimdine_reaction_ids = ["R_NTDUMP"]

    println("> Sodium-Potassium Pump")
    na_k_pump_reaction_ids = ["R_NaKt", "R_NAt"]

    println("> Other reactions")
    other_reaction_ids = ["R_GUAPRT"]

    println("> Transporters")
    transporter_reactions_ids = [
        # Original transporters
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

        # Expanded transporters from Bordbar 2016
        "R_AKGtec",
        "R_ARGtec",
        "R_CAATPS1",
        "R_CAMPtec",
        "R_CGMPtec",
        "R_FRUt1r",
        # "R_FUMtr",  # Zero flux
        "R_GSNt",
        "R_HCYSte",
        # "R_MALt",  # Zero flux
        "R_MANt1r",
        "R_METtec",
        "R_PTRCtex2",
        "R_SPMDtex2",
        "R_SPRMt2",
        # "R_URATEt",  # Zero flux
        "R_UREAt",
        "R_URIt",
        # "R_XANt",  # Zero flux
        "R_GTHOXABCte",
        "R_GLY_Cl_2Nat",

        # Additional transporters for AS-3 not listed above
        # "R_CITt",  # Zero flux
        "R_Clt",
    ]

    # println(transporter_reactions_ids)

    if !isnothing(exchanges)
        println("> Adding exchanges")
    else
        println("> Skipping exchanges")
    end

    exchange_reactions_ids = isnothing(exchanges) ? [] : exchanges

    println("> Collecting reactions and discovering metabolites")

    all_reaction_ids = [
        glycolysis_reaction_ids
        rl_shunt_reaction_ids
        ppp_reaction_ids
        purine_metabolism_reaction_ids
        met_salvage_reaction_ids
        citric_reaction_ids
        arg_pro_reaction_ids
        nucleotide_reaction_ids
        glutamate_reaction_ids
        glutathione_reaction_ids
        urea_reaction_ids
        folate_reaction_ids
        glycine_serine_threonine_reaction_ids
        fructose_mannose_reaction_ids
        pyrimdine_reaction_ids
        other_reaction_ids
        transporter_reactions_ids
        exchange_reactions_ids
        na_k_pump_reaction_ids
    ]

    discovered_metabolite_ids::Vector{String} = []
    for reaction_id ∈ all_reaction_ids
        for metabolite_id ∈ keys(rbc_gem.reactions[reaction_id].stoichiometry)
            push!(discovered_metabolite_ids, metabolite_id)
        end
    end

    println("Discovered $(length(discovered_metabolite_ids)) metabolites.")

    model = Model()

    for discovered_metabolite_id ∈ discovered_metabolite_ids
        copied_metabolite = deepcopy(rbc_gem.metabolites[discovered_metabolite_id])
        model.metabolites[discovered_metabolite_id] = copied_metabolite
    end

    println("> Adding reactions and exchanges to model")

    for reaction_id ∈ all_reaction_ids
        bounds_override = find_flux_bounds_overrides(flux_bounds_overrides_df, reaction_id)
        if !isnothing(bounds_override)
            rxn = deepcopy(rbc_gem.reactions[reaction_id])
            lower_bound, upper_bound = bounds_override
            rxn = setproperties(rxn; lower_bound = lower_bound, upper_bound = upper_bound)
            model.reactions[reaction_id] = rxn
        else
            model.reactions[reaction_id] = deepcopy(rbc_gem.reactions[reaction_id])
        end
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

    # println(model.reactions["R_LOAD_ATP"])

    println("> Adding NADH load")

    # Load due to methemoglobin reduction via CytB5
    model.reactions["R_LOAD_NADH"] = Reaction(
        name = "LOAD_NADH",
        stoichiometry = Dict("M_nadh_c" => -1.0, "M_h_c" => 1.0, "M_nad_c" => 1.0),
        objective_coefficient = 1.0,
        lower_bound = 0.0,
        upper_bound = 1.0,
    )

    # println(model.reactions["R_LOAD_NADH"])

    println("> Adding NADPH load")

    # Load due to glutathione reduction from GSSG to GSH
    model.reactions["R_LOAD_NADPH"] = Reaction(
        name = "LOAD_NADPH",
        stoichiometry = Dict("M_nadph_c" => -1.0, "M_h_c" => 1.0, "M_nadp_c" => 1.0),
        objective_coefficient = 1.0,
        lower_bound = 0.0,
        upper_bound = 1.0,
    )

    # println(model.reactions["R_LOAD_NADPH"])

    println("> Setting NaKt load")
    model.reactions["R_NaKt"].objective_coefficient = 1.0

    return model
end

end

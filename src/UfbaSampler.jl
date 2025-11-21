module UfbaSampler

using Distributed
@everywhere using Pkg
@everywhere Pkg.activate(".")
@info "Distributed.jl nprocs: $(nprocs())"
@everywhere using COBREXA, HiGHS, JuMP, MathOptInterface

import ConstraintTrees as C
import SBMLFBCModels as S
import AbstractFBCModels as A
import AbstractFBCModels: stoichiometry
import AbstractFBCModels.CanonicalModel: Model, Reaction, Metabolite, Gene, Coupling
using CSV
using DataFrames
using DataFramesMeta

export create_3p_model, sample_fluxes, ufba, load_metabolite_bounds, query_metabolite_bounds

function create_3p_model()
    println("############################################################")
    println("# LOAD RBC-GEM                                             #")
    println("############################################################")

    rbc_gem_path = joinpath("input", "RBC-GEM.xml")
    rbc_gem = load_model(S.SBMLFBCModel, rbc_gem_path, A.CanonicalModel.Model)

    println("Metabolites: $(length(rbc_gem.metabolites))")
    println("Reactions: $(length(rbc_gem.reactions))")

    println("############################################################")
    println("# GLYCOLYSIS REACTIONS                                     #")
    println("############################################################")

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

    println(glycolysis_reaction_ids)

    println("############################################################")
    println("# RL SHUNT                                                 #")
    println("############################################################")

    rl_shunt_reaction_ids = ["R_DPGM", "R_DPGase"]
    println(rl_shunt_reaction_ids)

    println("############################################################")
    println("# PENTOSE PHOSPHATE PATHWAY                                #")
    println("############################################################")

    ppp_reaction_ids =
        ["R_G6PDH2", "R_PGL", "R_GND", "R_RPI", "R_RPE", "R_TKT1", "R_TALA", "R_TKT2"]

    println(ppp_reaction_ids)

    println("############################################################")
    println("# PURINE METABOLISM                                        #")
    println("############################################################")

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

    println(purine_metabolism_reaction_ids)

    println("############################################################")
    println("# TRANSPORTERS                                             #")
    println("############################################################")

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

    println(transporter_reactions_ids)

    println("############################################################")
    println("# EXCHANGES                                                #")
    println("############################################################")

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

    println(exchange_reactions_ids)

    println("############################################################")
    println("# COLLECT REACTIONS IDS                                    #")
    println("############################################################")

    all_reaction_ids = [
        glycolysis_reaction_ids
        rl_shunt_reaction_ids
        ppp_reaction_ids
        purine_metabolism_reaction_ids
        transporter_reactions_ids
        exchange_reactions_ids
    ]

    println("############################################################")
    println("# DISCOVER METABOLITES                                     #")
    println("############################################################")

    discovered_metabolite_ids::Vector{String} = []
    for reaction_id ∈ all_reaction_ids
        for metabolite_id ∈ keys(rbc_gem.reactions[reaction_id].stoichiometry)
            push!(discovered_metabolite_ids, metabolite_id)
        end
    end

    println("Discovered $(length(discovered_metabolite_ids)) metabolites.")

    println("############################################################")
    println("# CREATE THREE PATHWAY MODEL                               #")
    println("############################################################")

    model = Model()

    println(model)

    println("############################################################")
    println("# ADD DISCOVERED METABOLITES                               #")
    println("############################################################")

    for discovered_metabolite_id ∈ discovered_metabolite_ids
        model.metabolites[discovered_metabolite_id] =
            deepcopy(rbc_gem.metabolites[discovered_metabolite_id])
        println("$discovered_metabolite_id")
    end

    println("############################################################")
    println("# ADD REACTIONS AND EXCHANGES                              #")
    println("############################################################")

    for reaction_id ∈ all_reaction_ids
        model.reactions[reaction_id] = deepcopy(rbc_gem.reactions[reaction_id])
        lower_bound = model.reactions[reaction_id].lower_bound
        upper_bound = model.reactions[reaction_id].upper_bound
        println(reaction_id, " (", lower_bound, ", ", upper_bound, ")")
    end

    println("############################################################")
    println("# ADD ATP LOAD                                             #")
    println("############################################################")

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

    println("############################################################")
    println("# ADD NADH LOAD                                            #")
    println("############################################################")

    # Load due to methemoglobin reduction via CytB5
    model.reactions["R_LOAD_NADH"] = Reaction(
        name = "LOAD_NADH",
        stoichiometry = Dict("M_nadh_c" => -1.0, "M_h_c" => 1.0, "M_nad_c" => 1.0),
        objective_coefficient = 1.0,
        lower_bound = 0.0,
        upper_bound = 1.0,
    )

    println(model.reactions["R_LOAD_NADH"])

    println("############################################################")
    println("# ADD NADPH LOAD                                           #")
    println("############################################################")

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

function load_metabolite_bounds()
    metabolite_bounds_filename = joinpath("output", "concentration_rates.csv")
    metabolite_bounds_df = CSV.read(metabolite_bounds_filename, DataFrame)
    return metabolite_bounds_df
end

function query_metabolite_bounds(metabolite_bounds_df, additive, metabolite, final_time)
    query_df = @rsubset(
        metabolite_bounds_df,
        :additive == additive,
        :metabolite == metabolite,
        :final_time == final_time
    )
    if nrow(query_df) > 0
        lower_bound = query_df[1, :lower_bound]
        upper_bound = query_df[1, :upper_bound]
        return (lower_bound, upper_bound)
    else
        return nothing
    end
end

function ufba(
    model::A.AbstractFBCModel,
    metabolite_bounds_df::DataFrame,
    additive::String,
    final_time::Int64
)
    println("\n############################################################")
    println("# uFBA Sampling $additive, final time: $final_time")
    println("############################################################")

    # Place bounds for Sv = b_lb, Sv = b_ub
    ct = flux_balance_constraints(model)
    unfound_metabolites = []
    found_metabolites = []
    for k ∈ keys(ct.flux_stoichiometry)
        short_metabolite_id = string(k)[3:end]
        bounds = query_metabolite_bounds(metabolite_bounds_df, additive, short_metabolite_id, final_time)
        if isnothing(bounds)
            push!(unfound_metabolites, short_metabolite_id)
            ct.flux_stoichiometry[k].bound = C.Between(-1000.0, 1000.0)
        else
            push!(found_metabolites, short_metabolite_id)
            lb, ub = bounds
            ct.flux_stoichiometry[k].bound = C.Between(lb, ub)
        end
    end
    println("\n>>>>>>>>> FOUND METABOLITES <<<<<<<<<")
    display(found_metabolites)
    println("\n>>>>>>>>> UNFOUND METABOLITES <<<<<<<<<")
    display(unfound_metabolites)
    println("\n>>>>>>>>> CONSTRAINT TREE <<<<<<<<<")
    C.pretty(ct)
    println("\n>>>>>>>>> SIMPLE OPTIMIZATION ATTEMPT <<<<<<<<<")
    objective_flux = optimized_values(
        ct;
        objective = ct.objective.value,
        output = ct.objective,
        optimizer = HiGHS.Optimizer,
        settings = [],
    )
    if isnothing(objective_flux)
        println("Simple optimization failed")
        return nothing
    else
        println("Simple optimization succeeded!")
        println("\n>>>>>>>>> FLUX SAMPLING <<<<<<<<<")
        samples_df = sample_fluxes(model)
        return samples_df
    end
end

"""
    sample_fluxes(model; n_chains::Int64, tolerance::Float64)

Sample the allowable flux space of the `model`. Use the `julia -p X...` -p command line option to set the number of workers for this operation.

# Arguments
1. `model`: Model to be sampled.
2. `n_chains::Int64`: The number of chains to calculate, with each chain producing ~126 samples. Defaults to 10 chains.
3. `tolerance::Float64`: The tolerance bounds on the objective.

# Returns
`DataFrame`
1. Returns a `DataFrame` with each reaction as a column and each row a flux sample.
"""
function sample_fluxes(model; n_chains::Int64 = 10, tolerance::Float64 = 0.99)
    s = flux_sample(
        model,
        optimizer = HiGHS.Optimizer,
        objective_bound = relative_tolerance_bound(tolerance),
        n_chains = n_chains,
        workers = workers(),
        collect_iterations = [10],
    )
    s_dict = Dict()
    for reaction_id ∈ keys(s)
        s_dict[reaction_id] = s[reaction_id]
    end
    DataFrame(s_dict)
end

end

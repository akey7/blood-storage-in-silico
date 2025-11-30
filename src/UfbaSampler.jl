module UfbaSampler

using Distributed
@everywhere using Pkg
@everywhere Pkg.activate(".")
@info "Distributed.jl nprocs: $(nprocs())"
@everywhere using COBREXA, HiGHS, JuMP, MathOptInterface

using Base.Iterators
import ConstraintTrees as C
import SBMLFBCModels as S
import AbstractFBCModels as A
import AbstractFBCModels: stoichiometry
import AbstractFBCModels.CanonicalModel: Model, Reaction, Metabolite, Gene, Coupling
using CSV
using DataFrames
using DataFramesMeta
using CairoMakie

export create_3p_model,
    sample_fluxes,
    ufba_all_additives_all_times,
    load_metabolite_bounds,
    query_metabolite_bounds,
    histograms_for_reaction_in_additive,
    plot_all_histograms,
    fba,
    add_sinks_for_unmatched_metabolites!,
    find_metabolite_matches,
    is_metabolite_in_exchange,
    add_case_1_constraints!

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
        lower_bound = query_df[1, :lb]
        upper_bound = query_df[1, :ub]
        return (lower_bound, upper_bound)
    else
        return nothing
    end
end

function fba(model::A.AbstractFBCModel; n_chains::Int64 = 10)
    println("\n############################################################")
    println("# STANDARD FBA SAMPLING                                    #")
    println("############################################################")

    println("\n>>>>>>>>> SIMPLE OPTIMIZATION ATTEMPT <<<<<<<<<")
    solution = flux_balance_analysis(model; optimizer = HiGHS.Optimizer)
    if isnothing(solution)
        println("Simple optimization failed")
        return nothing, nothing
    else
        println("Simple optimization succeeded!")
        display(solution.fluxes)
        println("\n>>>>>>>>> FLUX SAMPLING <<<<<<<<<")
        samples_df = sample_fluxes(model; n_chains = n_chains)
        return solution, samples_df
    end
end

function is_metabolite_in_exchange(model::A.AbstractFBCModel, metabolite::AbstractString)
    exchange_substring = "EX_$(metabolite[1:end-2])"
    for rxn in keys(model.reactions)
        if contains(rxn, exchange_substring)
            return true
        end
    end
    return false
end

function find_metabolite_matches(
    model::A.AbstractFBCModel,
    metabolite_bounds_df::DataFrame,
    additive::AbstractString,
    final_time::Int64,
)
    println("\n############################################################")
    println("# MATCHING METABOLITES FROM $additive, t_f = $final_time")
    println("############################################################")

    ct = flux_balance_constraints(model)
    status_rows = []
    found_count = 0
    not_found_count = 0
    in_exchange_count = 0
    for k ∈ keys(ct.flux_stoichiometry)
        short_metabolite_id = string(k)[3:end]
        bounds = query_metabolite_bounds(
            metabolite_bounds_df,
            additive,
            short_metabolite_id,
            final_time,
        )
        if isnothing(bounds)
            status_row = (
                additive = additive,
                metabolite = short_metabolite_id,
                status = "not found",
            )
            push!(status_rows, status_row)
            not_found_count += 1
            # ct.flux_stoichiometry[k].bound = C.Between(-1000.0, 1000.0)
        elseif is_metabolite_in_exchange(model, short_metabolite_id)
            status_row = (
                additive = additive,
                metabolite = short_metabolite_id,
                status = "in exchange",
            )
            in_exchange_count += 1
        else
            status_row =
                (additive = additive, metabolite = short_metabolite_id, status = "found")
            push!(status_rows, status_row)
            found_count += 1
            # lb, ub = bounds
            # ct.flux_stoichiometry[k].bound = C.Between(lb, ub)
        end
    end
    metabolite_status_df = DataFrame(status_rows)
    println(
        "Found $found_count, in exchange $in_exchange_count, not found $not_found_count",
    )
    return metabolite_status_df
end

function add_sinks_for_unmatched_metabolites!(
    model::A.AbstractFBCModel,
    metabolite_status_df::DataFrame,
    additive::AbstractString,
)
    println("\n############################################################")
    println("# ADD SINKS FOR UNMATCHED METABOLITES                      #")
    println("############################################################")

    not_found_df = @chain metabolite_status_df begin
        @rsubset(:status == "not found", :additive == additive)
        @select(:metabolite)
    end
    for metabolite in sort(unique(not_found_df.metabolite))
        sink_up_name = "R_SK_UP_$metabolite"
        sink_up = Reaction(
            name = sink_up_name,
            stoichiometry = Dict("M_$(metabolite)" => -1.0),
            lower_bound = -1000.0,
            upper_bound = 0.0,
        )
        model.reactions[sink_up_name] = sink_up
        display(sink_up)
        sink_down_name = "R_SK_DOWN_$metabolite"
        sink_down = Reaction(
            name = sink_down_name,
            stoichiometry = Dict("M_$(metabolite)" => 1.0),
            lower_bound = 0.0,
            upper_bound = 1000.0,
        )
        model.reactions[sink_down_name] = sink_down
        display(sink_down)
    end
end

function add_case_1_constraints!(model::A.AbstractFBCModel)
    rxn_ids = keys(model.reactions)
    relaxation_sinks =
        [rxn for rxn in rxn_ids if contains(rxn, "R_SK_UP") || contains(rxn, "R_SK_DOWN")]
    ct = flux_balance_constraints(model)
    ct *= :case_1^C.ConstraintTree()
    display(ct)
end

function ufba_additive_at_final_time(
    model::A.AbstractFBCModel,
    metabolite_bounds_df::DataFrame,
    additive::AbstractString,
    final_time::Int64;
    n_chains::Int64 = 10,
)
    println("\n############################################################")
    println("# uFBA Sampling $additive, final time: $final_time")
    println("############################################################")

    # Place bounds for Sv = b_lb, Sv = b_ub
    ct = flux_balance_constraints(model)
    status_rows = []
    for k ∈ keys(ct.flux_stoichiometry)
        short_metabolite_id = string(k)[3:end]
        bounds = query_metabolite_bounds(
            metabolite_bounds_df,
            additive,
            short_metabolite_id,
            final_time,
        )
        if isnothing(bounds)
            status_row = (
                additive = additive,
                metabolite = short_metabolite_id,
                status = "not found",
            )
            push!(status_rows, status_row)
            ct.flux_stoichiometry[k].bound = C.Between(-1000.0, 1000.0)
        else
            status_row =
                (additive = additive, metabolite = short_metabolite_id, status = "found")
            push!(status_rows, status_row)
            lb, ub = bounds
            ct.flux_stoichiometry[k].bound = C.Between(lb, ub)
        end
    end
    metabolite_status_df = DataFrame(status_rows)
    # println("\n>>>>>>>>> CONSTRAINT TREE <<<<<<<<<")
    # C.pretty(ct)
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
        return nothing, metabolite_status_df
    else
        println("Simple optimization succeeded!")
        println("\n>>>>>>>>> FLUX SAMPLING <<<<<<<<<")
        samples_df = sample_fluxes(model; n_chains = n_chains)
        return samples_df, metabolite_status_df
    end
end

function ufba_all_additives_all_times(
    model::A.AbstractFBCModel,
    metabolite_bounds_df::DataFrame;
    n_chains::Int64 = 10,
)
    println("\n############################################################")
    println("# uFBA: QUEUEING ADDITIVES AND FINAL TIMES                 #")
    println("############################################################")

    additives = unique(metabolite_bounds_df.additive)
    final_times = unique(metabolite_bounds_df.final_time)
    pairs = product(additives, final_times)
    status_rows = []
    pair_results = []
    metabolite_status_dfs = []
    println("Number of pairs: ", length(pairs))
    for (additive, final_time) in pairs
        pair_result, metabolite_status_df = ufba_additive_at_final_time(
            model,
            metabolite_bounds_df,
            additive,
            final_time;
            n_chains = n_chains,
        )
        push!(metabolite_status_dfs, metabolite_status_df)
        if isnothing(pair_result)
            status = (additive = additive, final_time = final_time, status = "fail")
            push!(status_rows, status)
        else
            status = (additive = additive, final_time = final_time, status = "ok")
            push!(status_rows, status)
            pair_result[!, :additive] .= additive
            pair_result[!, :final_time] .= final_time
            push!(pair_results, pair_result)
        end
    end
    sampling_df = vcat(pair_results...)
    status_df = DataFrame(status_rows)
    all_metabolite_status_df = vcat(metabolite_status_dfs...)
    return sampling_df, status_df, all_metabolite_status_df
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
    println("N Chains: $n_chains")
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

function histograms_for_reaction_in_additive(long_sampling_df, additive, reaction_id)
    plt_df = @chain long_sampling_df begin
        @rsubset(:additive == additive, :reaction_id == reaction_id)
        select(:final_time, :flux)
    end
    title = "$additive $reaction_id"
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Flux (mM/week)", ylabel = "Density", title = title)
    final_times = [2, 3, 4, 5, 6]
    colors = [:dodgerblue, :orange, :blueviolet, :crimson, :deeppink]
    for (final_time, color) in zip(final_times, colors)
        hist_df = @rsubset(plt_df, :final_time == final_time)
        hist!(
            ax,
            hist_df.flux;
            bins = 50,
            color = (color, 0.33),
            label = string(final_time),
        )
    end
    axislegend(ax)
    return fig
end

function plot_all_histograms(sampling_df)
    println("\n############################################################")
    println("# uFBA: PLOTTING HISTOGRAMS                                #")
    println("############################################################")

    long_sampling_df = stack(
        sampling_df,
        Not([:additive, :final_time]),
        variable_name = :reaction_id,
        value_name = :flux,
    )
    additives = unique(long_sampling_df.additive)
    reaction_ids = unique(long_sampling_df.reaction_id)
    pairs = product(additives, reaction_ids)
    n_pairs = length(pairs)
    for (i, (additive, reaction_id)) in enumerate(pairs)
        fig = histograms_for_reaction_in_additive(long_sampling_df, additive, reaction_id)
        filename = joinpath("output", "uFBA_histograms", "$additive $(reaction_id).png")
        save(filename, fig)
        println("Wrote $i of $n_pairs: $filename")
    end
end

end

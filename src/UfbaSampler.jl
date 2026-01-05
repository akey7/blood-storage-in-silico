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
using ThreadsX
using OrderedCollections

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
    case_3_constraint_tree!,
    list_objectives_in_model,
    optimize_case_3,
    display_jump_results,
    analyze_case_3,
    make_ufba_models_for_additives_and_times,
    execute_all_ufba_jobs,
    load_base_rbc_gem,
    map_reaction_ids_to_reaction_strings,
    extract_case3_sinks

function load_base_rbc_gem()
    println("> Loading RBC-GEM")
    rbc_gem_path = joinpath("input", "RBC-GEM.xml")
    rbc_gem = load_model(S.SBMLFBCModel, rbc_gem_path, A.CanonicalModel.Model)

    # println("Metabolites: $(length(rbc_gem.metabolites))")
    # println("Reactions: $(length(rbc_gem.reactions))")

    return rbc_gem
end

function create_3p_model(
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

function map_reaction_ids_to_reaction_strings(model::A.AbstractFBCModel)
    result = OrderedDict()
    for rxn_id in sort(string.(keys(model.reactions)))
        stoi = model.reactions[rxn_id].stoichiometry
        rxn = model.reactions[rxn_id]
        lhs = replace(
            join(
                [!isapprox(v, -1.0) ? "$(abs(v)) $k" : k for (k, v) in stoi if v < 0],
                " + ",
            ),
            "M_" => "",
        )
        rhs = replace(
            join(
                [!isapprox(v, 1.0) ? "$(abs(v)) $k" : k for (k, v) in stoi if v > 0],
                " + ",
            ),
            "M_" => "",
        )
        if rxn.lower_bound < 0.0 && isapprox(rxn.upper_bound, 0.0)
            result[rxn_id] = "$lhs <-- $rhs ($(rxn.lower_bound), $(rxn.upper_bound))"
        elseif isapprox(rxn.lower_bound, 0.0) && rxn.upper_bound > 0.0
            result[rxn_id] = "$lhs --> $rhs ($(rxn.lower_bound), $(rxn.upper_bound))"
        else
            result[rxn_id] = "$lhs <-> $rhs ($(rxn.lower_bound), $(rxn.upper_bound))"
        end
    end
    return result
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
    @info "Standard FBA sampling, N chains $n_chains"
    println("> Simple optimization attempt")
    solution = flux_balance_analysis(model; optimizer = HiGHS.Optimizer)
    if isnothing(solution)
        println("Simple optimization failed")
        return nothing, nothing
    else
        println("Simple optimization succeeded!")
        display(solution.fluxes)
        println("> Flux sampling")
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
    @info "Matching metabolites, additive: $additive, final_time: $final_time"
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
    prune_zero_sinks::Union{Vector{String},Nothing},
)
    if isnothing(prune_zero_sinks)
        @info "Add sinks for unmatched metabolites, not pruning any sinks"
    else
        @info "Add sinks for unmatched metabolites, pruning $(length(prune_zero_sinks))"
    end

    not_found_df = @chain metabolite_status_df begin
        @rsubset(:status == "not found", :additive == additive)
        @select(:metabolite)
    end
    for metabolite in sort(unique(not_found_df.metabolite))
        sink_up_name = "R_UNKNOWN_SK_UP_$metabolite"
        if !isnothing(prune_zero_sinks) && sink_up_name in prune_zero_sinks
            println("Skipping zero flux sink $sink_up_name")
        else
            sink_up = Reaction(
                name = sink_up_name,
                stoichiometry = Dict("M_$(metabolite)" => -1.0),
                lower_bound = -1000.0,
                upper_bound = 0.0,
            )
            model.reactions[sink_up_name] = sink_up
            # display(sink_up)
        end
        sink_down_name = "R_UNKNOWN_SK_DOWN_$metabolite"
        if !isnothing(prune_zero_sinks) && sink_down_name in prune_zero_sinks
            println("Skipping zero flux sink $sink_down_name")
        else
            sink_down = Reaction(
                name = sink_down_name,
                stoichiometry = Dict("M_$(metabolite)" => -1.0),
                lower_bound = 0.0,
                upper_bound = 1000.0,
            )
            model.reactions[sink_down_name] = sink_down
            # display(sink_down)
        end
    end
end

function list_objectives_in_model(model::A.AbstractFBCModel)
    obj_coeffs = A.AbstractFBCModels.objective(model)
    rxn_ids = keys(model.reactions)
    for (id, coeff) in zip(rxn_ids, obj_coeffs)
        if !isapprox(coeff, 0.0)
            println("Reaction $id has objective coefficient: $coeff")
        end
    end
end

function case_3_constraint_tree!(
    model::A.AbstractFBCModel,
    metabolite_status_df::DataFrame,
    additive::AbstractString,
)
    ct = flux_balance_constraints(model)
    flux_ids = collect(keys(ct.fluxes))
    unfound_metabolites =
        @rsubset(metabolite_status_df, :status == "not found", :additive == additive)
    unfound_metabolite_ids = Symbol.(unique(unfound_metabolites.metabolite))

    abs_flux_vars = C.variables(keys = flux_ids, bounds = C.Between(0.0, Inf))
    ct_flux_pos = C.zip(abs_flux_vars, ct.fluxes) do abs_flux_var, flux
        C.Constraint(abs_flux_var.value + flux.value, C.Between(0.0, Inf))
    end
    ct_flux_neg = C.zip(abs_flux_vars, ct.fluxes) do abs_flux_var, flux
        C.Constraint(abs_flux_var.value - flux.value, C.Between(0.0, Inf))
    end

    abs_metabolite_vars =
        C.variables(keys = unfound_metabolite_ids, bounds = C.Between(0.0, Inf))
    ct_metabolite_pos =
        C.zip(abs_metabolite_vars, ct.flux_stoichiometry) do abs_metabolite_var, stoi
            C.Constraint(abs_metabolite_var.value + stoi.value, C.Between(0.0, Inf))
        end
    ct_metabolite_neg =
        C.zip(abs_metabolite_vars, ct.flux_stoichiometry) do abs_metabolite_var, stoi
            C.Constraint(abs_metabolite_var.value - stoi.value, C.Between(0.0, Inf))
        end

    ct =
        ct +
        C.ConstraintTree(:abs_flux_vars => abs_flux_vars) +
        C.ConstraintTree(:flux_pos => ct_flux_pos) +
        C.ConstraintTree(:flux_neg => ct_flux_neg) +
        C.ConstraintTree(:abs_metabolite_vars => abs_metabolite_vars) +
        C.ConstraintTree(:metabolite_pos => ct_metabolite_pos) +
        C.ConstraintTree(:metabolite_neg => ct_metabolite_neg)

    ct.objective = C.Constraint(
        C.sum(a.value for (rxn_id, a) in abs_flux_vars; init = 0.0) +
        C.sum(a.value for (metabolite_id, a) in abs_metabolite_vars; init = 0.0),
    )

    return ct
end

function optimize_case_3(ct::C.ConstraintTree, objective::C.LinearValue)
    @info "Optimizing case 3"

    num_vars = C.variable_count(ct)
    model = JuMP.Model(HiGHS.Optimizer)
    JuMP.@variable(model, x[1:num_vars])
    JuMP.@objective(model, JuMP.MIN_SENSE, C.substitute(objective, x))

    C.traverse(ct) do c
        b = c.bound
        if b isa C.EqualTo
            JuMP.@constraint(model, C.substitute(c.value, x) == b.equal_to)
        elseif b isa C.Between
            val = C.substitute(c.value, x)
            isinf(b.lower) || JuMP.@constraint(model, val >= b.lower)
            isinf(b.upper) || JuMP.@constraint(model, val <= b.upper)
        end
    end

    JuMP.optimize!(model)
    if is_solved_and_feasible(model)
        println("Case 3 optimization success!")
        result_ct = deepcopy(ct)
        var_values = JuMP.value.(model[:x])
        solution_tree = C.substitute_values(result_ct, var_values)
        return solution_tree
    else
        println("Case 3 optimization failed")
        return nothing
    end
end

function analyze_case_3(case_3_optimize_result::C.Tree{Float64})
    zero_case3_sinks = [
        k for (k, v) in case_3_optimize_result.fluxes if
        isapprox(v, 0.0) && contains(string(k), "R_UNKNOWN_SK")
    ]
    nonzero_case3_sinks = [
        k for (k, v) in case_3_optimize_result.fluxes if
        !isapprox(v, 0.0) && contains(string(k), "R_UNKNOWN_SK")
    ]
    sink_status_rows = []
    for zero_case3_sink in zero_case3_sinks
        row = (sink = zero_case3_sink, is_non_zero = false)
        push!(sink_status_rows, row)
    end
    for nonzero_case3_sink in nonzero_case3_sinks
        row = (sink = nonzero_case3_sink, is_non_zero = true)
        push!(sink_status_rows, row)
    end
    sink_status_df = DataFrame(sink_status_rows)
    return zero_case3_sinks, nonzero_case3_sinks, sink_status_df
end

function execute_ufba_job(job, n_chains = 10)
    additive = job.additive
    final_time = job.final_time
    pruned_model = job.pruned_model
    metabolite_bounds_df = job.metabolite_bounds_df
    @info "execute_ufba_job: additive: $additive, final_time: $final_time"
    ct = flux_balance_constraints(pruned_model)
    for k in keys(ct.flux_stoichiometry)
        short_metabolite_id = string(k)[3:end]
        bounds = query_metabolite_bounds(
            metabolite_bounds_df,
            additive,
            short_metabolite_id,
            final_time,
        )
        if isnothing(bounds)
            ct.flux_stoichiometry[k].bound = C.EqualTo(0.0)
        else
            lb, ub = bounds
            ct.flux_stoichiometry[k].bound = C.Between(lb, ub)
        end
    end
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
        println("> Flux sampling")
        samples_df = sample_fluxes(pruned_model; n_chains = n_chains)
        samples_df[!, :additive] .= additive
        samples_df[!, :final_time] .= final_time
        return samples_df
    end
end

function execute_all_ufba_jobs(jobs, n_chains = 10)
    all_sampling_dfs_1 = map(jobs) do job
        execute_ufba_job(job, n_chains)
    end
    all_sampling_dfs_2 = [df for df in all_sampling_dfs_1 if !isnothing(df)]
    status_rows = vcat(
        eachrow([
            (
                additive = job.additive,
                final_time = job.final_time,
                status = isnothing(sdf) ? "fail" : "ok",
            ) for (job, sdf) in zip(jobs, all_sampling_dfs_1)
        ])...,
    )
    status_df = DataFrame(status_rows)
    return vcat(all_sampling_dfs_2...), status_df
end

function make_ufba_models_for_additives_and_times(
    metabolite_bounds_df::DataFrame,
    n_models::Int64,
)
    base_rbc_gem = load_base_rbc_gem()
    final_times = sort(unique(metabolite_bounds_df.final_time))
    additives = sort(unique(metabolite_bounds_df.additive))
    pairs =
        n_models == -1 ? collect(product(additives, final_times)) :
        collect(product(additives, final_times))[1:n_models]
    n_pairs = length(pairs)
    result = map(enumerate(pairs)) do p
        (i, (additive, final_time)) = p
        @info "make_ufba_models_for_additives_and_times: $i of $n_pairs"
        full_model = create_3p_model(base_rbc_gem; add_exchanges = false)
        metabolite_status_df = find_metabolite_matches(
            full_model,
            metabolite_bounds_df,
            additive,
            final_time,
        )
        add_sinks_for_unmatched_metabolites!(
            full_model,
            metabolite_status_df,
            additive,
            nothing,
        )
        ct = case_3_constraint_tree!(full_model, metabolite_status_df, additive)
        case_3_optimize_result_ct = optimize_case_3(ct, ct.objective.value)
        if isnothing(case_3_optimize_result_ct)
            @error "Failed to optimize case 3 for additive: $additive, final_time: $final_time"
        end
        zero_case3_sinks, nonzero_case3_sinks, sink_status_df =
            analyze_case_3(case_3_optimize_result_ct)
        sink_status_df[!, :additive] .= additive
        sink_status_df[!, :final_time] .= final_time
        pruned_model = create_3p_model(base_rbc_gem; add_exchanges = false)
        add_sinks_for_unmatched_metabolites!(
            pruned_model,
            metabolite_status_df,
            additive,
            string.(zero_case3_sinks),
        )
        (
            additive = additive,
            final_time = final_time,
            full_model = deepcopy(full_model),
            pruned_model = deepcopy(pruned_model),
            sink_status_df = sink_status_df,
            metabolite_bounds_df = deepcopy(metabolite_bounds_df),
            zero_case3_sinks = zero_case3_sinks,
            nonzero_case3_sinks = nonzero_case3_sinks,
        )
    end
    return result
end

function extract_case3_sinks(ufba_jobs)
    rows = []
    for ufba_job in ufba_jobs
        for zero_case3_sink in ufba_job.zero_case3_sinks
            row = (
                additive = ufba_job.additive,
                final_time = ufba_job.final_time,
                sink = zero_case3_sink,
                status = "zero",
            )
            push!(rows, row)
        end
        for nonzero_case3_sink in ufba_job.nonzero_case3_sinks
            row = (
                additive = ufba_job.additive,
                final_time = ufba_job.final_time,
                sink = nonzero_case3_sink,
                status = "nonzero",
            )
            push!(rows, row)
        end
    end
    return DataFrame(rows)
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

end

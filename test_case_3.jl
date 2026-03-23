using CSV
using DataFrames
using DataFramesMeta
using COBREXA
using HiGHS
using Distributed

include("src/UfbaSampler.jl")
using .UfbaSampler
include("src/FbaModelBuilder.jl")
using .FbaModelBuilder
include("src/PruningOptimizations.jl")
using .PruningOptimizations
include("src/MetaboliteBounds.jl")
using .MetaboliteBounds

function init_workers!(; project::AbstractString = Base.active_project())
    for p in workers()
        Distributed.remotecall_eval(
            Main,
            p,
            quote
                import Pkg
                Pkg.activate($project)
                using COBREXA, HiGHS, JuMP, MathOptInterface
                # include("src/UfbaSampler.jl")
                # using .UfbaSampler
            end,
        )
    end
    return nothing
end

init_workers!()
workers_config = workers()

@info "Loading base RBC-GEM"
base_rbc_gem = load_base_rbc_gem()

@info "Loading metabolite bounds"
metabolite_bounds_df = load_metabolite_bounds()

@info "Loading flux bounds overrides"
flux_bounds_overrides_df = load_flux_bounds_overrides()

@info "Loading metabolite measurement opt outs"
metabolites_to_ignore = load_metabolite_measurement_opt_outs()

@info "Create FBA model and map metabolites onto that model"
fba_model, _ = create_fba_model(
    base_rbc_gem;
    exchanges = default_exchanges(),
    flux_bounds_overrides_df = flux_bounds_overrides_df,
)

@info "Reference additive and time point"
additive = "01-Ctrl AS3"
final_time = 2
println("additive: $additive, final_time: $final_time")

@info "Adding sinks to model"
first_model, _ = create_fba_model(
    base_rbc_gem;
    exchanges = default_exchanges(),
    flux_bounds_overrides_df = flux_bounds_overrides_df,
)
metabolite_status_df =
    find_metabolite_matches(first_model, metabolite_bounds_df, additive, final_time)
first_sink_specifications = (
    metabolite_status_df = metabolite_status_df,
    additive = additive,
    prune_zero_sinks = nothing,
    sink_opt_outs = nothing,
)
first_added_sink_ids =
    add_sinks_for_unmatched_metabolites!(first_model, first_sink_specifications)

@info "Add metabolite bounds to ConstraintTree"
second_ct = flux_balance_constraints(first_model)
metabolites_to_ignore = nothing
measured_unmeasured = add_metabolite_bounds_to_constraint_tree!(
    second_ct,
    metabolite_bounds_df,
    additive,
    final_time;
    metabolites_to_ignore = metabolites_to_ignore,
)
unmeasured_metabolite_ids = measured_unmeasured.unmeasured_metabolites

@info "Case 3: Sinks list"
sink_ids = find_sinks_on_ct(second_ct)
display(first(sink_ids, 10))

@info "Case 3: Optimize"
optimize_case_3_ok_fail, optimize_case_3_result =
    optimize_case_3(second_ct, unmeasured_metabolite_ids)
if optimize_case_3_ok_fail == :fail
    display(optimize_case_3_result)
    error("Case 3 optimization failed.")
end

@info "Case 3: Prune zero sinks"
third_fba_model, _ = create_fba_model(
    base_rbc_gem;
    exchanges = default_exchanges(),
    flux_bounds_overrides_df = flux_bounds_overrides_df,
)
case_3_analysis = analyze_pruning_optimization(optimize_case_3_result)
prune_zero_sinks = case_3_analysis.prune
third_sink_specifications = (
    metabolite_status_df = metabolite_status_df,
    additive = additive,
    prune_zero_sinks = prune_zero_sinks,
    sink_opt_outs = nothing,
)
third_added_sink_ids =
    add_sinks_for_unmatched_metabolites!(third_fba_model, third_sink_specifications)
third_ct = flux_balance_constraints(third_fba_model)
add_metabolite_bounds_to_constraint_tree!(
    third_ct,
    metabolite_bounds_df,
    additive,
    final_time;
    metabolites_to_ignore = metabolites_to_ignore,
)
n_pruned_sinks = length(case_3_analysis.prune)
n_kept_sinks = length(case_3_analysis.keep)
println("Pruned $n_pruned_sinks, kept $n_kept_sinks")

@info "Zeroth test case: Sampling, no sinks, no metabolite bounds"
zeroth_ct = flux_balance_constraints(fba_model)
zeroth_samples, _ = sample_fluxes(zeroth_ct, workers_config; n_chains = 5)
zeroth_n_zero_fluxes, _ = count_n_all_zero_fluxes(zeroth_samples)
println("zeroth_n_zero_fluxes: $zeroth_n_zero_fluxes")

@info "First test case: Sampling, all sinks (no pruning), no metabolite bounds"
first_ct = flux_balance_constraints(first_model)
first_samples, _ = sample_fluxes(first_ct, workers_config; n_chains = 5)
first_n_zero_fluxes, _ = count_n_all_zero_fluxes(first_samples)
println("first_n_zero_fluxes: $first_n_zero_fluxes")

@info "Second test case: all sinks (no pruning), all metabolite bounds"
second_samples, _ = sample_fluxes(second_ct, workers_config; n_chains = 5)
second_n_zero_fluxes, _ = count_n_all_zero_fluxes(second_samples)
println("second_n_zero_fluxes: $second_n_zero_fluxes")

@info "Third test case: Pruned sinks, all metabolite bounds"
third_samples, _ = sample_fluxes(third_ct, workers_config; n_chains = 5)
third_n_zero_fluxes, _ = count_n_all_zero_fluxes(third_samples)
println("third_n_zero_fluxes: $third_n_zero_fluxes")

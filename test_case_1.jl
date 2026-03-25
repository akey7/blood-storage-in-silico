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
relax_quantile = 0.1
println("additive: $additive, final_time: $final_time, relax_quantile: $relax_quantile")

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
# print_sinks_in_model(first_model)

@info "Add metabolite bounds to ConstraintTree"
second_ct = flux_balance_constraints(first_model)
add_metabolite_bounds_to_constraint_tree!(
    second_ct,
    metabolite_bounds_df,
    additive,
    final_time;
    metabolites_to_ignore = metabolites_to_ignore,
)
# print_metabolite_bounds_on_constraint_tree(second_ct)

@info "Optimize for Case 1 and create pruned model"
# optimize_case_1_result = optimize_case_1(
#     case1_ct;
#     force_first_sink_on = true,
#     force_first_sink_lb = 0.1,
#     print_objective_value = true,
# )
optimize_case_1_ok_fail, optimize_case_1_result = optimize_case_1(second_ct)
if optimize_case_1_ok_fail == :fail
    display(optimize_case_1_result)
    error("Case 1 optimization failed. Conflicting constraints are listed above. Stopping.")
end
case_1_analysis = analyze_pruning_optimization(optimize_case_1_result)

@info "Prune zero sinks according to Case 1"
third_fba_model, _ = create_fba_model(
    base_rbc_gem;
    exchanges = default_exchanges(),
    flux_bounds_overrides_df = flux_bounds_overrides_df,
)
prune_zero_sinks = case_1_analysis.prune
third_sink_specifications = (
    metabolite_status_df = metabolite_status_df,
    additive = additive,
    prune_zero_sinks = prune_zero_sinks,
    sink_opt_outs = nothing,
)
third_added_sink_ids =
    add_sinks_for_unmatched_metabolites!(third_fba_model, third_sink_specifications)
# println("Added the following sinks")
# display(third_added_sink_ids)
third_ct = flux_balance_constraints(third_fba_model)
measured_unmeasured = add_metabolite_bounds_to_constraint_tree!(
    third_ct,
    metabolite_bounds_df,
    additive,
    final_time;
    metabolites_to_ignore = metabolites_to_ignore,
    relax_quantile = relax_quantile,
)
unmeasured_metabolite_ids = measured_unmeasured.unmeasured_metabolites

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

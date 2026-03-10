include("src/UfbaSampler.jl")
using .UfbaSampler
include("src/FbaModelBuilder.jl")
using .FbaModelBuilder
include("src/PruningOptimizations.jl")
using .PruningOptimizations

@info "Metabolite bounds"
metabolites_bounds_df = load_metabolite_bounds()
display(first(metabolites_bounds_df, 10))

@info "Flux bounds overrides"
flux_bounds_overrides_df = load_flux_bounds_overrides()

@info "Create FBA model and map metabolites onto that model"
fba_model, fba_model_metabolites_df = create_fba_model(
    load_base_rbc_gem();
    exchanges = default_exchanges(),
    flux_bounds_overrides_df = flux_bounds_overrides_df,
)
mapping_additive = "01-Ctrl AS3"
metabolite_status_df =
    find_metabolite_matches(fba_model, metabolites_bounds_df, mapping_additive, 2)
display(first(metabolite_status_df, 10))
first_sink_specifications = (
    metabolite_status_df = metabolite_status_df,
    additive = mapping_additive,
    prune_zero_sinks = nothing,
    sink_opt_outs = nothing,
)
first_added_sink_ids = add_sinks_for_unmatched_metabolites!(fba_model, first_sink_specifications)
display(first(first_added_sink_ids, 10))

@info "Test case 1 optimization"
case_1_ct = case_1_constraint_tree(fba_model)
case_1_result = optimize_case_1(case_1_ct, case_1_ct.objective.value)
display(case_1_result)

include("src/UfbaSampler.jl")
using .UfbaSampler
include("src/FbaModelBuilder.jl")
using .FbaModelBuilder

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

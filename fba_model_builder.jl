using ArgParse

include("src/UfbaSampler.jl")
using .UfbaSampler
include("src/FbaModelBuilder.jl")
using .FbaModelBuilder

s = ArgParseSettings()
@add_arg_table! s begin
    "--nchains"
    help = "Number of chains during sampling"
    arg_type = Int64
    default = 10
end
n_chains = parse_args(s)["nchains"]

init_workers!()

base_rbc_gem = load_base_rbc_gem()
fba_model = create_fba_model(base_rbc_gem; exchanges = nothing)
result = fba(fba_model; n_chains = n_chains)

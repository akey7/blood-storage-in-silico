module ModelGraph

using COBREXA
import SBMLFBCModels
import AbstractFBCModels as A
import SBMLFBCModels as S
using JSON3
using Graphs
import Graphs.Parallel
using MetaGraphs
using DataFrames
using DataFramesMeta
using Arpack
using Clustering
using DataStructures
using LinearAlgebra
using ProgressMeter
using OrderedCollections

export model_to_dictionaries,
    metabolite_id_keys_reaction_id_values,
    metabolite_id_to_other_side_metabolite_ids,
    metabolite_ids_to_metabolite_names,
    reaction_ids_to_reaction_names,
    metabolite_names_to_reaction_names,
    make_graph,
    dfs_from_metabolite_id,
    betweenness_centrality_of_vertices,
    run_dfs_plan,
    summarize_dfs_plan_result,
    spectral_cluster_metabolite_graph,
    remove_metabolite_ids,
    prepare_spectral_clustering_for_yaml,
    find_isolated_vertices,
    load_ufba_models,
    make_graphs_for_ufba_models,
    run_all_dfs_plans,
    enrich_visited_reactions_df

"""
    load_ufba_models()

Loads uFBA models written during sampling from the `output/ufba_models` folder. Returns a dictionary with tuples of additive and final_time mapped to uFBA models.

# Returns
`Dict{Tuple{String,Int64},A.CanonicalModel.Model}`

Returns a dictionary mapping additives and final times to uFBA models.
"""
function load_ufba_models()
    ufba_model_folder = joinpath("output", "ufba_models")
    sbml_files = filter(
        f -> endswith(lowercase(f), ".xml") && isfile(joinpath(ufba_model_folder, f)),
        readdir(ufba_model_folder),
    )
    sbml_paths = joinpath.(ufba_model_folder, sbml_files)
    n_sbml_paths = length(sbml_paths)
    prog_sbml_paths = Progress(n_sbml_paths, "Loading uFBA models")
    ufba_models::Dict{Tuple{String,Int64},A.CanonicalModel.Model} = Dict()
    for sbml_path in sbml_paths
        bn = replace(basename(sbml_path), ".xml" => "", "uFBA " => "")
        additive, final_time_str = split(bn, "_")
        final_time = parse(Int64, final_time_str)
        ufba_model = load_model(S.SBMLFBCModel, sbml_path, A.CanonicalModel.Model)
        ufba_models[(additive, final_time)] = ufba_model
        next!(prog_sbml_paths)
    end
    return ufba_models
end

"""
    model_to_dictionaries(filename::String)

Load a given filename into an `A.CanonicalModel.Model`, which reveals `.reactions` and `.metabolites` properties, among other things, that are plain Julia dictionaries.

# Argument
1. `filename::String`: Path to the model file.

# Returns
`A.CanonicalModel.Model` containing the loaded model.
"""
model_to_dictionaries(filename::String) = load_model(filename, A.CanonicalModel.Model)

"""
    metabolite_id_keys_reaction_id_values(model::A.CanonicalModel.Model)

Create a dictionary with metabolite_ids as keys and vectors of reaction ids as values. This enable finding all reactions a metabolite participates in.

# Argument
1. `model::A.CanonicalModel.Model`: Model as loaded by `model_to_dictionaries()`.

# Returns
- `Dict{String,Vector{String}}`
Mapping metabolite_ids to the reactions they are involved in.
"""
function metabolite_id_keys_reaction_id_values(model::A.CanonicalModel.Model)
    reaction_stoichiometry =
        Dict(key => value.stoichiometry for (key, value) ∈ model.reactions)
    result =
        Dict(metabolite_string => Vector() for metabolite_string in keys(model.metabolites))
    for (key, value) ∈ reaction_stoichiometry
        for key0 ∈ keys(value)
            push!(result[key0], key)
        end
    end
    result
end

"""
    metabolite_id_to_other_side_metabolite_ids(model::A.CanonicalModel.Model, metabolite_id_to_reaction_ids::Dict{String, Vector{String}})

This function does two things:
1. For each `metabolite_id`, find each metabolite_id on the other side of the reactions, as found by the opposite sign of the stoichiometric coefficients (if oppostie side metabolites are present). These "other side" metabolites are further filtered by their stoichiometric coefficient to line up with the bounds of the reaction for which the edge is being specified. 
2. For reactions that have metabolites on both sides of the reaction (normal reactions and transporters), maps pairs of metabolite ids to the reaction that connects them.
3. For reactions that just have a metabolite on one side of the reactions (exchanges and sinks), maps single metabolites to the reaction id of the corresponding exchange/sink.

The first two operations are used by [`make_graph`](@ref BloodStorageInSilico.ModelGraph.make_graph) to make graphs of models. At this time, the results of the thirs operation are not used because that would lead to inconsistent defintions of vertices and edges. See [`make_graph`](@ref BloodStorageInSilico.ModelGraph.make_graph) for more information on this.

# Arguments
1. `model::A.CanonicalModel.Model`: Model to inspect
2. `metabolite_id_to_reaction_ids::Dict{String, Vector{Any}}`: Output of `metabolite_id_keys_reaction_id_values()`.

# Returns
`Dict{Symbol,Dict{String,Vector{Any}}}`

1. `:metabolite_to_metabolties`: metabolite_ids as keys and vectors of metabolite ids as values. The metabolite ids in the vectors are on the "other side" of the reactions according to stoichiometry in the reactions. For example, the metabolite id keys that are products will have lists of their reactants across all reactions in which they participate.
2. `:metabolite_pairs_to_reactions`: Metabolite pairs as keys and reaction ids as values. This can be used to label the edges in the metabolite graph. It shows what reaction pair metabolite ids together.
3. `:metabolites_to_exchanges_sinks`: This dictionary contains information about exchanges and sinks. Maps single metabolites ids as keys to single reaction ids as values.
"""
function metabolite_id_to_other_side_metabolite_ids(
    model::A.CanonicalModel.Model,
    metabolite_id_to_reaction_ids::Dict{String,Vector{Any}},
)
    metabolite_to_metabolites::Dict{String,Vector{String}} = Dict()
    metabolite_pairs_to_reactions::Dict{Tuple{String,String},String} = Dict()
    metabolites_to_exchanges_sinks::Dict{String,String} = Dict()
    for (metabolite_id, reaction_ids) ∈ metabolite_id_to_reaction_ids
        metabolite_to_metabolites[metabolite_id] = Vector()
        for reaction_id ∈ reaction_ids
            stoichiometry = model.reactions[reaction_id].stoichiometry
            lower_bound = model.reactions[reaction_id].lower_bound
            upper_bound = model.reactions[reaction_id].upper_bound
            other_side_sign = -1 * sign(stoichiometry[metabolite_id])
            if isapprox(lower_bound, 0.0) && upper_bound > 0.0
                other_metabolite_ids = [
                    other_metabolite_id for
                    (other_metabolite_id, coeff) ∈ stoichiometry if
                    coeff == other_side_sign && coeff > 0
                ]
                append!(metabolite_to_metabolites[metabolite_id], other_metabolite_ids)
            elseif lower_bound < 0.0 && isapprox(upper_bound, 0.0)
                other_metabolite_ids = [
                    other_metabolite_id for
                    (other_metabolite_id, coeff) ∈ stoichiometry if
                    coeff == other_side_sign && coeff < 0
                ]
                append!(metabolite_to_metabolites[metabolite_id], other_metabolite_ids)
            else
                other_metabolite_ids = [
                    other_metabolite_id for
                    (other_metabolite_id, coeff) ∈ stoichiometry if coeff == other_side_sign
                ]
                append!(metabolite_to_metabolites[metabolite_id], other_metabolite_ids)
            end
            if length(other_metabolite_ids) >= 1
                for other_metabolite_id ∈ other_metabolite_ids
                    metabolite_pair = (metabolite_id, other_metabolite_id)
                    metabolite_pairs_to_reactions[(metabolite_pair)] = reaction_id
                end
            else
                metabolites_to_exchanges_sinks[metabolite_id] = reaction_id
            end
        end
    end
    Dict(
        :metabolite_to_metabolites => metabolite_to_metabolites,
        :metabolite_pairs_to_reactions => metabolite_pairs_to_reactions,
        :metabolites_to_exchanges_sinks => metabolites_to_exchanges_sinks,
    )
end

"""
    metabolite_ids_to_metabolite_names(model::A.CanonicalModel.Model)

Maps metabolite ids to names and vice versa.

# Argument
1. `model::A.CanonicalModel.Model`: Model as loaded by `model_to_dictionaries()`.

# Returns
`Dict{Symbol,Dict{String,String}}` with two keys:
1. `:ids_to_names`, a 1:1 mapping of metabolite ids to names and
2. `:names_to_ids`, a 1:1 mapping of metabolite names to ids.
"""
function metabolite_ids_to_metabolite_names(model::A.CanonicalModel.Model)
    ids_to_names::Dict{String,String} = Dict(
        strip(metabolite_id) => strip(model.metabolites[metabolite_id].name) for
        metabolite_id ∈ keys(model.metabolites)
    )
    names_to_ids::Dict{String,String} = Dict(value => key for (key, value) ∈ ids_to_names)
    Dict(:ids_to_names => ids_to_names, :names_to_ids => names_to_ids)
end

"""
    reaction_ids_to_reaction_names(model::A.CanonicalModel.Model)

Maps reactions ids to names and vice versa.

# Argument
1. `model::A.CanonicalModel.Model`: Model as loaded by `model_to_dictionaries()`.

# Returns
`Dict{Symbol,Dict{String,String}}` with two keys:
1. `:ids_to_names`, a 1:1 mapping of reaction ids to names and
2. `:names_to_ids`, a 1:1 mapping of reaction names to ids.
"""
function reaction_ids_to_reaction_names(model::A.CanonicalModel.Model)
    ids_to_names::Dict{String,String} = Dict(
        strip(reaction_id) => strip(model.reactions[reaction_id].name) for
        reaction_id ∈ keys(model.reactions)
    )
    names_to_ids::Dict{String,String} = Dict(value => key for (key, value) ∈ ids_to_names)
    Dict(:ids_to_names => ids_to_names, :names_to_ids => names_to_ids)
end

"""
    function metabolite_names_to_reaction_names(metabolites_to_reactions::Dict{String,Vector{Any}}, metabolite_ids_to_names::Dict{String,String}, reaction_ids_to_names::Dict{String,String})

Creates a one-to-many map of metabolite ids to the reaction ids they participate in.

# Arguments
1. `metabolites_to_reactions::Dict{String,Vector{Any}}`: Map of metabolite ids to vectors of reaction ids.
2. `metabolite_ids_to_names::Dict{String,String}`: Map of metabolite ids to metabolite names.
3. `reaction_ids_to_names::Dict{String,String}`: Map of reaction ids to reaction names.

# Returns
`Dict{String,Vector{String}}`
Each key is a metabolite id, each value is a vector of reaction ids that the metabolite id is found in.
"""
function metabolite_names_to_reaction_names(
    metabolites_to_reactions::Dict{String,Vector{Any}},
    metabolite_ids_to_names::Dict{String,String},
    reaction_ids_to_names::Dict{String,String},
)
    result::Dict{String,Vector{String}} = Dict()
    for (metabolite_id, reaction_ids) ∈ metabolites_to_reactions
        key = metabolite_ids_to_names[metabolite_id]
        if !haskey(result, key)
            result[key] = Vector()
        end
        append!(
            result[key],
            [reaction_ids_to_names[reaction_id] for reaction_id ∈ reaction_ids],
        )
    end
    result
end

"""
    make_graph(model::A.CanonicalModel.Model)

Creates a graph and related data of the given model, with metabolite ids as the vertices and reactions as the edges.

Note: This function does not add edges for exchanges, though it does add edges for transports. The reason it doesn't add edges for edges for exchanges or sinks is because exchanges and sinks only include one metabolite and would therefore necessitate a different definition for an edge than the rest of the metabolites/reactions. So to maintain consistency sinks/exchanges are skipped.

# Argument
1. `model::A.CanonicalModel.Model`: Model loaded by `model_to_dictionaries()`

# Returns
`Dict{Symbol,Any}`
The dictionary contains the following keys:
1. `:metabolite_ids_to_ints`: Dictionary with metabolite id strings as keys and integer ids as values
2. `:ints_to_metabolite_ids`: Vector of strings, with each string being a metabolite id at the appropriate index.
3. `:adj_matrix`: NxN adjacency matrix, with N being number of metabolites ids.
4. `:mdg`: `MetaDiGraph` of the metaoblic network, useful for visualization and interoperatbility with the Graphs.jl ecosystem.
"""
function make_graph(model::A.CanonicalModel.Model)
    N = length(model.metabolites)
    metabolite_integer::Int64 = 1
    metabolite_ids_to_ints::Dict{String,Int64} = Dict()
    ints_to_metabolite_ids::Vector{String} = Vector{String}(undef, N)
    reaction_ids_to_vertices::Dict{String,Tuple{Int64,Int64}} = Dict()
    for metabolite_id ∈ sort(collect(keys(model.metabolites)))
        if haskey(metabolite_ids_to_ints, metabolite_id)
            continue
        end
        metabolite_ids_to_ints[metabolite_id] = metabolite_integer
        ints_to_metabolite_ids[metabolite_integer] = metabolite_id
        metabolite_integer += 1
    end
    metabolite_connections = metabolite_id_to_other_side_metabolite_ids(
        model,
        metabolite_id_keys_reaction_id_values(model),
    )
    adj_matrix = zeros(Int64, N, N)
    mdg = MetaDiGraph(N)
    set_indexing_prop!(mdg, :metabolite_id)
    metabolite_to_metabolites = metabolite_connections[:metabolite_to_metabolites]
    metabolite_pairs_to_reactions = metabolite_connections[:metabolite_pairs_to_reactions]
    for (vertex_metabolite_id, neighbor_metabolite_ids) ∈ metabolite_to_metabolites
        vertex_idx = metabolite_ids_to_ints[vertex_metabolite_id]
        set_prop!(mdg, vertex_idx, :metabolite_id, vertex_metabolite_id)
        for neighbor_metabolite_id ∈ neighbor_metabolite_ids
            neighbor_idx = metabolite_ids_to_ints[neighbor_metabolite_id]
            reaction_id = metabolite_pairs_to_reactions[(
                vertex_metabolite_id,
                neighbor_metabolite_id,
            )]
            if contains(lowercase(reaction_id), "r_ex")
                continue
            end
            adj_matrix[vertex_idx, neighbor_idx] = 1
            add_edge!(mdg, vertex_idx, neighbor_idx)
            set_prop!(mdg, Edge(vertex_idx, neighbor_idx), :reaction_id, reaction_id)
            reaction_ids_to_vertices[reaction_id] = (vertex_idx, neighbor_idx)
        end
    end
    Dict(
        :metabolite_ids_to_ints => metabolite_ids_to_ints,
        :ints_to_metabolite_ids => ints_to_metabolite_ids,
        :adj_matrix => adj_matrix,
        :N => N,
        :mdg => mdg,
        :model => model,
        :metabolite_pairs_to_reactions => metabolite_pairs_to_reactions,
        :reaction_ids_to_vertices => reaction_ids_to_vertices,
    )
end

"""
    make_graphs_for_ufba_models(ufba_models::Dict{Tuple{String,Int64},A.CanonicalModel.Model})

Makes graphs for all uFBA models provided as loaded by [`load_ufba_models`](@ref BloodStorageInSilico.ModelGraph.load_ufba_models)

# Arguments
1. `ufba_models::Dict{Tuple{String,Int64},A.CanonicalModel.Model}`: Loaded models.

# Returns
`Dict{Tuple{String,Int64},Dict{Symbol,Any}}`

Returns a dictionary mapping tuples of additive and final time to graph data from the [`make_graph`](@ref BloodStorageInSilico.ModelGraph.make_graph)
"""
function make_graphs_for_ufba_models(
    ufba_models::Dict{Tuple{String,Int64},A.CanonicalModel.Model},
)
    result::Dict{Tuple{String,Int64},Dict{Symbol,Any}} = Dict()
    n_models = length(keys(ufba_models))
    prog = Progress(n_models, "Creating graphs from uFBA models")
    for ((additive, final_time), ufba_model) in ufba_models
        result[(additive, final_time)] = make_graph(ufba_model)
        next!(prog)
    end
    return result
end

"""
    dfs_from_metabolite_id(graph_data::Dict{Symbol,Any}, metabolite_id::String, max_depth::Int64, metabolite_ids_to_skip::Union{Vector{String},Nothing})

Perform a depth-first search of the metabolic network starting from the given metabolite_id and limited to max_depth through the network. Prints a warning if a metabolite that does not exist is requested to be skipped.

# Returns
`Dict{Symbol,Any}`
A dictionary with the following keys:
1. `:metabolite_id`: The starting metabolite id.
2. `:max_depth`: The maximum depth of the search.
3. `:all_visited`: All visited metabolite ids.
4. `:all_visited_reactions`: All traversed reaction ids.
5. `:paths_metabolite_ids`: Metabolite ids encoutered over all paths.
6. `:paths_reaction_ids`: Reaction ids encoutnered over all paths.

# Arguments
1. `graph_data::Dict{Symbol,Any}`: Graph data as returned by `construct_adjacency_matrix()`.
2. `metabolite_id::String`: The metabolite_id to start from.
3. `max_depth::Int64`: Max number of reactions to traverse through the network.
4.  `metabolite_ids_to_skip::Union{Vector{String},Nothing}`: If specified, a list of metabolite ids not to traverse in the DFS. Defaults to `nothing`, which will traverse everything.
"""
function dfs_from_metabolite_id(
    graph_data::Dict{Symbol,Any},
    metabolite_id::String,
    max_depth::Int64,
    metabolite_ids_to_skip::Union{Vector{String},Nothing} = nothing,
)
    metabolite_ids_to_ints = graph_data[:metabolite_ids_to_ints]
    if !haskey(metabolite_ids_to_ints, metabolite_id)
        throw(KeyError("$metabolite_id not found in graph."))
    end
    metabolite_ids_to_skip0 =
        isnothing(metabolite_ids_to_skip) ? [] : metabolite_ids_to_skip
    for metabolite_id_to_skip ∈ metabolite_ids_to_skip0
        if !haskey(metabolite_ids_to_ints, metabolite_id_to_skip)
            @warn "$metabolite_id_to_skip was configured to skipped, but was not found in graph. Ignoring."
        end
    end
    skip_vec = [
        metabolite_ids_to_ints[metabolite_id] for
        metabolite_id ∈ metabolite_ids_to_skip0 if
        haskey(metabolite_ids_to_ints, metabolite_id)
    ]
    N = graph_data[:N]
    visited = zeros(Int64, N)
    adj_matrix = graph_data[:adj_matrix]
    paths::Vector{Vector{Int64}} = []
    start_vertex = metabolite_ids_to_ints[metabolite_id]

    function dfs(vertex::Int64, hop::Int64, current_path::Vector{Any})
        if hop > max_depth
            return
        end
        push!(current_path, vertex)
        visited[vertex] = 1
        for other_vertex ∈ 1:N
            if other_vertex ∉ skip_vec &&
               adj_matrix[vertex, other_vertex] == 1 &&
               visited[other_vertex] != 1
                dfs(other_vertex, hop + 1, deepcopy(current_path))
            end
        end
        if length(current_path) > 0
            push!(paths, current_path)
        end
    end

    dfs(start_vertex, 0, [])
    ints_to_metabolite_ids = graph_data[:ints_to_metabolite_ids]
    all_visited = [
        ints_to_metabolite_ids[index] for (index, value) ∈ enumerate(visited) if value == 1
    ]

    metabolite_pairs_to_reactions = graph_data[:metabolite_pairs_to_reactions]
    paths_metabolite_ids::Vector{Vector{String}} = []
    paths_reaction_ids::Vector{Vector{String}} = []
    for path ∈ paths
        path_metabolite_id = [ints_to_metabolite_ids[index] for index ∈ path]
        push!(paths_metabolite_ids, path_metabolite_id)
        path_reaction_id::Vector{String} = []
        n_path_metabolite_id = length(path_metabolite_id)
        if n_path_metabolite_id > 1
            for i ∈ 1:(length(path_metabolite_id)-1)
                pair = (path_metabolite_id[i], path_metabolite_id[i+1])
                push!(path_reaction_id, metabolite_pairs_to_reactions[pair])
            end
            push!(paths_reaction_ids, path_reaction_id)
        else
            # @warn "path_metabolite_id length $n_path_metabolite_id, contents $path_metabolite_id, this could cause a reaction to be skipped. Check ignored metabolite ids!"
        end
    end
    all_visited_reactions = unique(reduce(vcat, paths_reaction_ids))
    result = Dict(
        :metabolite_id => metabolite_id,
        :max_depth => max_depth,
        :all_visited => all_visited,
        :all_visited_reactions => all_visited_reactions,
        :paths_metabolite_ids => paths_metabolite_ids,
        :paths_reaction_ids => paths_reaction_ids,
    )
    return result
end

"""
    betweenness_centrality_of_vertices(graph_data::Dict{Symbol,Any})

Calculate normalized betweeness and closeness centrality and return in a DataFrame sorted by metabolite name.

# Returns
A `DataFrame`, sorted by metabolite id, with the following columns
1. `metabolite_name`: Full name of the metabolite (not unique because of cytosolic and extracellular metabolites)
2. `metabolite_id`: Unique metabolite id of every vertex.
3. `betweeness_centrality_normalized`: Normalized betweeness centrality of that vertex.
4. `closeness_centrality_normalized`: Normalized closness centrality of that vertex.

# Arugment
- `graph_data::Dict{Symbol,Any}`: Graph data from `make_graph()`
"""
function betweenness_centrality_of_vertices(graph_data::Dict{Symbol,Any})
    mdg = graph_data[:mdg]
    ints_to_metabolite_ids = graph_data[:ints_to_metabolite_ids]
    model = graph_data[:model]
    metabolite_names = metabolite_ids_to_metabolite_names(model)[:ids_to_names]
    bc = Parallel.betweenness_centrality(mdg; normalize = true)
    cc = Parallel.closeness_centrality(mdg; normalize = true)
    df = DataFrame(
        [
            (
                metabolite_names[metabolite_id],
                metabolite_id,
                betweenness_centrality_normalized,
                closeness_centrality_normalized,
            ) for (
                metabolite_id,
                betweenness_centrality_normalized,
                closeness_centrality_normalized,
            ) ∈ zip(ints_to_metabolite_ids, bc, cc)
        ],
        [
            :metabolite_name,
            :metabolite_id,
            :betweenness_centrality_normalized,
            :closeness_centrality_normalized,
        ],
    )
    sort(df, :metabolite_name)
end

"""
    run_dfs_plan(dfs_plan::DataFrame, graph_data::Dict{Symbol,Any}, metabolite_ids_to_skip::Union{Vector{String},Nothing})

Runs a "DFS plan" which is a list of of metabolite ids and their corresponding max depths to search. It then runs DFS on every metabolite id to the given depth and returns a data structure with the results of each DFS.

# Arguments
1. `dfs_plan::DataFrame`: A DataFrame with two columns: `metabolite_id` is the metabolite id to search and `max_depth` is the max depth to traverse.
2. `graph_data::Dict{Symbol,Any}`: The graph data as returned by `make_graph()`
3. `metabolite_ids_to_skip::Union{Vector{String},Nothing}`: A list of metabolite ids to skip traversal of. Defaults to `nothing`, which means it won't skip metabolite ids.

# Returns
`Vector{Dict{Symbol,Any}}`
Returns a vector of dictionaries. Each dictionary has the following keys:
1. `:metabolite_id`: The starting metabolite id.
2. `:max_depth`: The maximum depth of the search.
3. `:all_visited`: All visited metabolite ids.
4. `:all_visited_reactions`: All traversed reaction ids.
5. `:paths_metabolite_ids`: Metabolite ids encoutered over all paths.
6. `:paths_reaction_ids`: Reaction ids encoutnered over all paths.
"""
function run_dfs_plan(
    dfs_plan::DataFrame,
    graph_data::Dict{Symbol,Any},
    metabolite_ids_to_skip::Union{Vector{String},Nothing} = nothing,
)
    metabolite_ids = dfs_plan[!, :metabolite_id]
    max_depths = dfs_plan[!, :max_depth]
    results::Vector{Dict{Symbol,Any}} = []
    for (metabolite_id, max_depth) ∈ zip(metabolite_ids, max_depths)
        # println("run_dfs_plan(): Search from $metabolite_id for max_depth of $max_depth")
        result = dfs_from_metabolite_id(
            graph_data,
            String(metabolite_id),
            max_depth,
            metabolite_ids_to_skip,
        )
        push!(results, result)
    end
    return results
end

"""
    summarize_dfs_plan_result(results::Vector{Dict{Symbol,Any}})

Summarizes the results of running a DFS plan. It summarizes by extract all visited metabolite ids and all visited reactions.

# Arguments
1. `results::Vector{Dict{Symbol,Any}}`: DFS plan search results from `run_dfs_plan()`

# Returns
`Dict{String,Dict{Symbol,Vector{String}}}`
1. Returns a dictionary mapping keys of metabolite ids to values of summarized search results.
"""
function summarize_dfs_plan_result(results::Vector{Dict{Symbol,Any}})
    summary::Dict{String,Dict{Symbol,Vector{String}}} = Dict()
    for result ∈ results
        metabolite_id = result[:metabolite_id]
        summary[metabolite_id] = Dict(
            :all_visited => result[:all_visited],
            :all_visited_reactions => result[:all_visited_reactions],
        )
    end
    summary
end

"""
    run_all_dfs_plans(model_graphs::Dict{Tuple{String,Int64},Dict{Symbol,Any}}, dfs_plan::DataFrame, metabolite_ids_to_skip::Union{Vector{String},Nothing} = nothing)

Run DFS plan acorss all uFBA models specified in the arguments with [`dfs_from_metabolite_id`](@ref BloodStorageInSilico.ModelGraph.dfs_from_metabolite_id).

# Arguments
1. `model_graphs::Dict{Tuple{String,Int64},Dict{Symbol,Any}}`: uFBA model data from [`make_graphs_for_ufba_models`](@ref BloodStorageInSilico.ModelGraph.make_graphs_for_ufba_models).
2. `dfs_plan::DataFrame`: DataFrame of the DFS plan to execute on each model graph. Should have columns `metabolite_id` (metabolite id to start from) and `max_depth` (max number of hops to traverse).
3. `metabolite_ids_to_skip::Union{Vector{String},Nothing} = nothing`: Metabolite ids to skip in the traversal. This is used because some metabolites (like water) are highly connected and therefore might not be iteresting to traverse.

# Returns
`NamedTuple`

Returns a named tuple of two DataFrames.
1. `visited_metabolite_df`: A DataFrame of metabolites visited in the traversals. This DataFrame contains columns `additive`, `final_time`, `start_metabolite_id`, `path_idx`, `hops`, and `visited_metabolite_id`.
2. `visited_reaction_df`: A DataFrame of reactions visited in the traversals. This DataFrame contains the columns `additive`, `final_time`, `start_metabolite_id`, `path_idx`, `hops`, and `reaction_id`.
"""
function run_all_dfs_plans(
    model_graphs::Dict{Tuple{String,Int64},Dict{Symbol,Any}},
    dfs_plan::DataFrame,
    metabolite_ids_to_skip::Union{Vector{String},Nothing} = nothing,
)
    visited_metabolite_rows = []
    visited_reaction_rows = []
    n_model_graphs = length(keys(model_graphs))
    prog = Progress(n_model_graphs, "Executing DFS plan for uFBA models")
    for ((additive, final_time), graph_data) ∈ model_graphs
        dfs_plan_results = run_dfs_plan(dfs_plan, graph_data, metabolite_ids_to_skip)
        for dfs_plan_result in dfs_plan_results
            start_metabolite_id = dfs_plan_result[:metabolite_id]
            paths_metabolite_ids = dfs_plan_result[:paths_metabolite_ids]
            paths_reaction_ids = dfs_plan_result[:paths_reaction_ids]
            for (path_idx, path_metabolite_ids) in enumerate(paths_metabolite_ids)
                visited_metabolite_ids = path_metabolite_ids
                for (idx, visited_metabolite_id) in enumerate(visited_metabolite_ids)
                    hops = idx - 1
                    visited_metabolite_row = (
                        additive = additive,
                        final_time = final_time,
                        start_metabolite_id = start_metabolite_id,
                        path_idx = path_idx,
                        hops = hops,
                        visited_metabolite_id = visited_metabolite_id,
                    )
                    push!(visited_metabolite_rows, visited_metabolite_row)
                end
            end
            for (path_idx, path_reaction_ids) in enumerate(paths_reaction_ids)
                visited_reaction_ids = path_reaction_ids
                for (idx, visited_reaction_id) in enumerate(visited_reaction_ids)
                    hops = idx - 1
                    visited_reaction_row = (
                        additive = additive,
                        final_time = final_time,
                        start_metabolite_id = start_metabolite_id,
                        path_idx = path_idx,
                        hops = hops,
                        reaction_id = visited_reaction_id,
                    )
                    push!(visited_reaction_rows, visited_reaction_row)
                end
            end
        end
        next!(prog)
    end
    visited_metabolite_df = DataFrame(visited_metabolite_rows)
    visited_reaction_df = DataFrame(visited_reaction_rows)
    println("run_all_dfs_plans()")
    return (
        visited_metabolite_df = @orderby(
            visited_metabolite_df,
            :additive,
            :final_time,
            :start_metabolite_id,
            :path_idx,
            :hops,
        ),
        visited_reaction_df = @orderby(
            visited_reaction_df,
            :additive,
            :final_time,
            :start_metabolite_id,
            :path_idx,
            :hops,
        ),
    )
end

"""
    enrich_visited_reactions_df(visited_reactions_df::DataFrame, rxn_ids_to_strings::OrderedDict{String,Any}, median_df::DataFrame)

"Enrich" the visited reactions DataFrame returned by [`run_all_dfs_plans`](@ref BloodStorageInSilico.ModelGraph.run_all_dfs_plans) by adding columns with median flux rate and reaction string to the visited reactions DataFrame.

# Arguments
1. `visited_reactions_df::DataFrame`: Visited reactions DataFrame
2. `rxn_ids_to_strings::OrderedDict{String,Any}` The mapping of reaction ids to reactions returned by [`map_reaction_ids_to_reaction_strings`](@ref BloodStorageInSilico.UfbaSampler.map_reaction_ids_to_reaction_strings)
3. `median_df::DataFrame`: Median flux DataFrame returned by [`calc_median_flux_df`](@ref BloodStorageInSilico.UfbaSamplerAnalysis.calc_median_flux_df)

# Returns
`DataFrame`

Returns a DataFrame with the following columns: 

1. `:additive`: The additive of the model.
2. `:final_time`: The final time of the model.
3. `:start_metabolite_id`: The starting metabolite id of the traversed path.
4. `:path_idx`: The path index, an integer (starting at 1), that specifies which depth-first path from the starting metabolite id is traversed for this row.
5. `:hops`: The number of hops across the graph (starting at 0 for a reaction containing the metabolite id) that this reaction is found at on the depth first path starting from the starting metabolite id.
6. `:reaction_id`: Reaction id of the traversed edge.
7. `:median_flux`: The median flux of the reaction on the traversed edge.
8. `:reaction_string`: The reaction string of the traversed edge.
"""
function enrich_visited_reactions_df(
    visited_reactions_df::DataFrame,
    rxn_ids_to_strings::OrderedDict{String,Any},
    median_df::DataFrame,
)
    reaction_ids = []
    reaction_strings = []
    for (reaction_id, vals) in rxn_ids_to_strings
        push!(reaction_ids, reaction_id)
        push!(reaction_strings, vals["rxn_string"])
    end
    reaction_map_df =
        DataFrame(reaction_id = reaction_ids, reaction_string = reaction_strings)
    result_df = @chain visited_reactions_df begin
        innerjoin(reaction_map_df; on = :reaction_id)
        innerjoin(median_df; on = [:additive, :final_time, :reaction_id])
        @orderby(:additive, :final_time, :start_metabolite_id, :path_idx, :hops)
        @select(
            :additive,
            :final_time,
            :start_metabolite_id,
            :path_idx,
            :hops,
            :reaction_id,
            :median_flux,
            :reaction_string
        )
    end
    return result_df
end

"""
    find_isolated_vertices(mdg::MetaDiGraph)

Print the metabolite ids of the isolated vertices (vertices with no neighbors) in the metabolite graph `mdg`. This function is used to inform the list of vertices to exclude from clustering because isolated vertices put zeros in the degree matrix which prevent finding eigenvectros of the symmeterized Laplacian in the spectral clustering step.

# Argument
1. `mdg::MetaDiGraph`: MetaDiGraph with properties on each vertex of metabolite ids associated with each respective vertex.
"""
function find_isolated_vertices(mdg::MetaDiGraph)
    W = adjacency_matrix(mdg)
    D = Diagonal(sum(W, dims = 2)[:])
    isolated_vertices = findall(x -> x == 0, diag(D))
    isolated_metabolite_ids = [get_prop(mdg, v, :metabolite_id) for v ∈ isolated_vertices]
    println(isolated_metabolite_ids)
end

"""
    remove_metabolite_ids(mdg::MetaDiGraph, metabolite_ids::Vector{String})

Remove vertices referenced their associated `metabolite_ids` from `mdg`. Returns a new graph leaving the original unchanged.

This function can be used to remove clusters of metabolites with very high connectivity (such as O2) that would interfere with spectral clustering's ability to cluster the other metabolites.

# Arguments
1. `mdg::MetaDiGraph`: Undirected graph with metabolites as vertices connected by edges that represent reactions.
2. `metabolite_ids::Vector{String}`: Metabolite ids to remove.

# Returns
`MetaDiGraph`
Returns a new graph with the given metabolite ids removed.
"""
function remove_metabolite_ids(mdg::MetaDiGraph, metabolite_ids::Vector{String})
    new_mdg = deepcopy(mdg)
    for metabolite_id ∈ metabolite_ids
        for v ∈ 1:nv(new_mdg)
            if get_prop(new_mdg, v, :metabolite_id) == metabolite_id
                rem_vertex!(new_mdg, v)
                break
            end
        end
    end
    new_mdg
end

"""
    spectral_cluster_metabolite_graph(mdg::MetaDiGraph; k::Int64)

Perform a spectral clustering of the metabolite graph `mdg` into `k` clusters. Return a dictionary with a bunch of key value pairs that are the results of the clustering.

# Arguments
1. `mdg::MetaDiGraph`: Undirected metabolite graph with metabolites as vertices connected by edges representing reactions.
2. `k::Int64`: How many clusters to split the metabolite graph into.

# Returns
`Dict{Symbol,Any}`
Returns a dictionary with the following keys:
1. `:mdg`: The original metabolite graph on which this clustering was performed.
2. `:cluster_to_vertices`: `OrderedDict{Int64,Vector{Int64}}` Keys are cluster ids and values are vectors of vertex ids within that cluster.
3. `:cluster_ids`: `Vector{Int64}` a vector of all cluster ids.
"""
function spectral_cluster_metabolite_graph(mdg::MetaDiGraph; k::Int64)
    W = adjacency_matrix(mdg)
    D = Diagonal(sum(W, dims = 2)[:])
    D_inv = Diagonal(1 ./ diag(D))
    L = Array(laplacian_matrix(mdg))
    L_rw = D_inv * L
    _, eigvecs = eigs(L_rw, nev = k)
    max_imag = maximum(abs.(imag.(eigvecs)))
    println("INFO: Maximum imaginary part: ", max_imag)
    eigvecs_real = real(eigvecs)
    clustering_result = kmeans(eigvecs_real', k)
    cluster_assignments = clustering_result.assignments
    vertex_to_cluster =
        Dict(vertex => cluster for (vertex, cluster) in enumerate(cluster_assignments))
    cluster_to_vertices::Dict{Int64,Vector{Int64}} = Dict(i => [] for i ∈ 1:k)
    for (vertex, assignment) ∈ enumerate(cluster_assignments)
        push!(cluster_to_vertices[assignment], vertex)
    end
    num_clusters_multiple_members = 0
    for (cluster, vertices) ∈ cluster_to_vertices
        if length(vertices) > 1
            num_clusters_multiple_members += 1
        end
    end
    cluster_ids = collect(keys(cluster_to_vertices))
    Dict(
        :cluster_to_vertices => cluster_to_vertices,
        :cluster_ids => cluster_ids,
        :mdg => mdg,
    )
end

"""
    prepare_spectral_clustering_for_yaml(clustering)

Return a dictionary data structure specifically designed for easy readability of output YAML or JSON.

# Arguments
1. `clustering`: Output of `spectral_cluster_metabolite_graph()`.

# Returns
`Dict{Symbol,Any}`
Returns a dictionary with the following keys:
1. `:cluster_sizes`: `OrderedDict{Int64,Int64}` keys are cluster ids and values are the number of vertices in that cluster.
2. `:cluster_assignments`: `OrderedDict{Int64,OrderedDict{String,Int64}}` Each top-level key is a cluster id, and each value is a nested `OrderedDict{String,Int64}` where the keys are metabolite ids and the values are the number of neighbors that metabolite id has on the graph.
"""
function prepare_spectral_clustering_for_yaml(clustering)
    all_clusters_sorted = sort(clustering[:cluster_ids])
    mdg = clustering[:mdg]
    cluster_sizes::OrderedDict{Int64,Int64} = OrderedDict()
    for cluster ∈ all_clusters_sorted
        cluster_sizes[cluster] = length(keys(clustering[:cluster_to_vertices][cluster]))
    end
    cluster_assignments =
        OrderedDict(cluster => OrderedDict() for cluster ∈ all_clusters_sorted)
    for (cluster, vertices) ∈ clustering[:cluster_to_vertices]
        cluster_assignments[cluster] = Dict(
            get_prop(mdg, v, :metabolite_id) => length(all_neighbors(mdg, v)) for
            v ∈ vertices
        )
    end
    Dict(:cluster_sizes => cluster_sizes, :cluster_assignments => cluster_assignments)
end

end

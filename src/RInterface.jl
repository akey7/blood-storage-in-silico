module RInterface

using DataFrames
using CSV

export export_correlation_dict_for_r

"""
    export_correlation_dict_for_r(
        corr_dict::Dict{Tuple{String,Int64},DataFrame},
        out_dir::AbstractString;
        row_col::Symbol = :row_variable,
        symmetry_atol::Real = 1e-8,
        diagonal_atol::Real = 1e-8,
        check_diagonal_is_one::Bool = true,
        export_upper_triangle_only::Bool = true,
        export_diagonal_in_long::Bool = false,
        missingstring::AbstractString = "NA",
    )

Validate and export a dictionary of correlation matrices for use in R.

Input format
------------
`corr_dict` must map `(additive, timepoint)` to a DataFrame where:
- one column, named by `row_col`, contains row reaction IDs
- all remaining columns are reaction IDs
- the DataFrame represents a square correlation matrix

Validation
----------
For every matrix, the function checks:
- `row_col` exists
- row IDs are unique
- column reaction IDs are unique
- the matrix is square
- row reaction IDs match column reaction IDs
- all matrices share the same reaction ID set and ordering
- the matrix is symmetric within `symmetry_atol`
- optionally, diagonal entries are near 1.0 within `diagonal_atol`

Outputs
-------
Creates `out_dir` containing:
- `manifest.csv`
- `correlations_long.csv`
- `matrices/` with one CSV per matrix

Returns
-------
Named tuple with:
- `manifest_df`
- `long_df`
- `reaction_ids`
"""
function export_correlation_dict_for_r(
    corr_dict::Dict{Tuple{String,Int64},DataFrame},
    out_dir::AbstractString;
    row_col::String = "row_variable",
    symmetry_atol::Real = 1e-8,
    diagonal_atol::Real = 1e-8,
    check_diagonal_is_one::Bool = true,
    export_upper_triangle_only::Bool = true,
    export_diagonal_in_long::Bool = false,
    missingstring::AbstractString = "NA",
)
    isempty(corr_dict) && error("corr_dict is empty.")

    # ---------- small helpers ----------
    safe_filename_part(x::AbstractString) = begin
        s = strip(x)
        s = replace(s, r"\s+" => "_")
        s = replace(s, r"[^A-Za-z0-9_\-]" => "")
        isempty(s) && error("Could not derive a safe filename token from '$x'.")
        s
    end

    function normalize_reaction_ids(df, row_col)
        row_ids = string.(df[!, row_col])
        col_syms = filter(!=(row_col), names(df))
        col_ids = string.(col_syms)
        return row_ids, col_syms, col_ids
    end

    function matrix_from_df(df, col_syms)
        n = nrow(df)
        p = length(col_syms)
        n == p || error(
            "Matrix is not square: nrow(df) = $n but number of reaction columns = $p.",
        )

        mat = Matrix{Float64}(undef, n, p)
        for (j, col) in enumerate(col_syms)
            colvec = df[!, col]
            for i = 1:n
                val = colvec[i]
                if ismissing(val)
                    error("Missing value found in column $(col), row $i.")
                end
                mat[i, j] = Float64(val)
            end
        end
        return mat
    end

    # ---------- set up output dirs ----------
    matrices_dir = joinpath(out_dir, "correlation_matrices_1")
    mkpath(matrices_dir)

    manifest_df = DataFrame(additive = String[], timepoint = Int[], filename = String[])

    long_chunks = DataFrame[]

    # ---------- establish reference reaction set/order ----------
    sorted_entries = sort(collect(corr_dict); by = x -> x[1])

    first_key, first_df = first(sorted_entries)
    row_col in names(first_df) ||
        error("Expected column $(row_col) in all DataFrames. Broke at $(first_key)")

    ref_row_ids, ref_col_syms, ref_col_ids = normalize_reaction_ids(first_df, row_col)

    length(unique(ref_row_ids)) == length(ref_row_ids) ||
        error("Duplicate reaction IDs found in $(row_col) for key $first_key.")
    length(unique(ref_col_ids)) == length(ref_col_ids) ||
        error("Duplicate reaction column names found for key $first_key.")

    Set(ref_row_ids) == Set(ref_col_ids) || error(
        "For key $first_key, row reaction IDs and column reaction IDs do not match as sets.",
    )

    ref_row_ids == ref_col_ids || error(
        "For key $first_key, row reaction IDs and column reaction IDs do not have the same order.",
    )

    reaction_ids = copy(ref_row_ids)

    # ---------- main loop ----------
    for ((additive, timepoint), df) in sorted_entries
        row_col in names(df) ||
            error("Expected column $(row_col) for key ($(additive), $(timepoint)).")

        row_ids, col_syms, col_ids = normalize_reaction_ids(df, row_col)

        length(unique(row_ids)) == length(row_ids) ||
            error("Duplicate row reaction IDs for key ($(additive), $(timepoint)).")
        length(unique(col_ids)) == length(col_ids) ||
            error("Duplicate reaction columns for key ($(additive), $(timepoint)).")

        nrow(df) == length(col_ids) || error(
            "Non-square matrix for key ($(additive), $(timepoint)): $(nrow(df)) rows vs $(length(col_ids)) reaction columns.",
        )

        Set(row_ids) == Set(col_ids) ||
            error("Row and column reaction IDs differ for key ($(additive), $(timepoint)).")

        row_ids == col_ids || error(
            "Row and column reaction ID order differs for key ($(additive), $(timepoint)). Reorder before export.",
        )

        row_ids == reaction_ids || error(
            "Reaction ordering differs from reference for key ($(additive), $(timepoint)).",
        )

        mat = matrix_from_df(df, col_syms)

        # Symmetry check
        issymmetric_ok = isapprox(mat, transpose(mat); atol = symmetry_atol, rtol = 0.0)
        issymmetric_ok || error(
            "Matrix is not symmetric within atol=$(symmetry_atol) for key ($(additive), $(timepoint)).",
        )

        # Diagonal check
        if check_diagonal_is_one
            diagvals = [mat[i, i] for i = 1:size(mat, 1)]
            all(isapprox.(diagvals, 1.0; atol = diagonal_atol, rtol = 0.0)) || error(
                "Diagonal is not all ones within atol=$(diagonal_atol) for key ($(additive), $(timepoint)).",
            )
        end

        # Write matrix CSV exactly as provided
        additive_token = safe_filename_part(additive)
        filename = "corr__$(additive_token)__week_$(timepoint).csv"
        filepath = joinpath(matrices_dir, filename)
        CSV.write(filepath, df; missingstring = missingstring)

        push!(
            manifest_df,
            (additive = additive, timepoint = timepoint, filename = filename),
        )

        # Build long-form chunk
        long_df = DataFrame(
            additive = String[],
            timepoint = Int[],
            reaction_i = String[],
            reaction_j = String[],
            correlation = Float64[],
        )

        for (i, row) in enumerate(eachrow(df))
            reaction_i = String(row[row_col])

            if export_upper_triangle_only
                j_start = export_diagonal_in_long ? i : i + 1
                j_stop = length(col_syms)
            else
                j_start = 1
                j_stop = length(col_syms)
            end

            if j_start <= j_stop
                for j = j_start:j_stop
                    reaction_j = String(col_syms[j])
                    corr_value = Float64(row[col_syms[j]])
                    push!(
                        long_df,
                        (
                            additive = additive,
                            timepoint = timepoint,
                            reaction_i = reaction_i,
                            reaction_j = reaction_j,
                            correlation = corr_value,
                        ),
                    )
                end
            end
        end

        push!(long_chunks, long_df)
    end

    long_df = reduce(vcat, long_chunks)

    CSV.write(joinpath(out_dir, "manifest.csv"), manifest_df; missingstring = missingstring)
    CSV.write(
        joinpath(out_dir, "correlations_long.csv"),
        long_df;
        missingstring = missingstring,
    )

    return (manifest_df = manifest_df, long_df = long_df, reaction_ids = reaction_ids)
end

end

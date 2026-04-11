module UfbaSamplerViz3D

using DataFrames
using StatsBase
using Statistics
using PlotlyJS

export stacked_flux_histograms_3d, stacked_flux_histogram_steps_3d_colored

#####################################################################
# BEGIN: BLOCKY HISTOGRAMS WITH A LINE BEHIND THEM                  #
#####################################################################

"""
    cuboid_trace(x0, x1, y0, y1, z0, z1; color="royalblue", opacity=0.7, name="", showlegend=false)

Create a PlotlyJS mesh3d cuboid spanning:
- x in [x0, x1]
- y in [y0, y1]
- z in [z0, z1]
"""
function cuboid_trace(
    x0,
    x1,
    y0,
    y1,
    z0,
    z1;
    color = "royalblue",
    opacity = 0.7,
    name = "",
    showlegend = false,
)
    # 8 vertices
    x = [x0, x1, x1, x0, x0, x1, x1, x0]
    y = [y0, y0, y1, y1, y0, y0, y1, y1]
    z = [z0, z0, z0, z0, z1, z1, z1, z1]

    # 12 triangles
    i = Int[0, 0, 4, 4, 0, 0, 1, 1, 2, 2, 3, 3]
    j = Int[1, 2, 5, 6, 1, 5, 2, 6, 3, 7, 0, 4]
    k = Int[2, 3, 6, 7, 5, 4, 6, 5, 7, 6, 4, 7]

    return mesh3d(
        x = x,
        y = y,
        z = z,
        i = i,
        j = j,
        k = k,
        color = color,
        opacity = opacity,
        flatshading = true,
        hoverinfo = "skip",
        name = name,
        showlegend = showlegend,
        showscale = false,
    )
end

"""
    xy_plane_trace(xmin, xmax, ymax, z0; color="lightgray", opacity=0.15)

Horizontal XY plane at z = z0.
"""
function xy_plane_trace(xmin, xmax, ymax, z0; color = "lightgray", opacity = 0.15)
    return PlotlyJS.surface(
        x = [xmin, xmax],
        y = [0.0, ymax],
        z = fill(z0, 2, 2),
        colorscale = [[0.0, color], [1.0, color]],
        opacity = opacity,
        hoverinfo = "skip",
        showscale = false,
        name = "XY plane",
        showlegend = false,
    )
end

"""
    xz_plane_trace(xmin, xmax, zmin, zmax; color="gainsboro", opacity=0.18)

Vertical XZ plane at y = 0.
"""
function xz_plane_trace(xmin, xmax, zmin, zmax; color = "gainsboro", opacity = 0.18)
    return PlotlyJS.surface(
        x = [xmin, xmax],
        y = fill(0.0, 2, 2),
        z = [zmin zmin; zmax zmax],
        colorscale = [[0.0, color], [1.0, color]],
        opacity = opacity,
        hoverinfo = "skip",
        showscale = false,
        name = "XZ plane",
        showlegend = false,
    )
end

"""
    stacked_flux_histograms_3d(
        flux_df::DataFrame;
        edges=nothing,
        nbins::Int=25,
        z_thickness::Real=0.18,
        bar_color::AbstractString="steelblue",
        bar_opacity::Real=0.72,
        plane_opacity::Real=0.15,
        use_bin_midpoint_for_mode::Bool=true,
        tie_method::Symbol=:first,
        fig_title::AbstractString="3D Flux Histograms Across Weeks",
    )

Create a 3D stacked histogram plot from a long DataFrame with columns:
- :final_time  (integer week labels, e.g. 2,3,4,5,6)
- :flux        (Float64 flux samples)

Geometry:
- X axis = histogram bins
- Y axis = bin counts
- Z axis = week
- Mode trajectory = line in the y = 0 plane across weeks

Arguments
---------
edges:
    Optional shared histogram bin edges. Strongly recommended for reproducibility.
    If `nothing`, edges are computed from all fluxes using a uniform range.

nbins:
    Used only when `edges === nothing`.

tie_method:
    How to resolve multiple modal bins:
    - :first  => first max bin
    - :mean   => average center of all max bins
"""
function stacked_flux_histograms_3d(
    flux_df::DataFrame;
    edges = nothing,
    nbins::Int = 25,
    z_thickness::Real = 0.18,
    bar_color::AbstractString = "steelblue",
    bar_opacity::Real = 0.72,
    plane_opacity::Real = 0.15,
    use_bin_midpoint_for_mode::Bool = true,
    tie_method::Symbol = :first,
    fig_title::AbstractString = "3D Flux Histograms Across Weeks",
)
    required_cols = ["final_time", "flux"]
    missing_cols = setdiff(required_cols, names(flux_df))
    isempty(missing_cols) || error("flux_df is missing required columns: $(missing_cols)")

    nrow(flux_df) > 0 || error("flux_df is empty")

    df = dropmissing(flux_df, [:final_time, :flux])
    nrow(df) > 0 || error("flux_df has no non-missing rows in :final_time and :flux")

    weeks = sort(unique(df.final_time))
    issubset(Set(weeks), Set(df.final_time)) ||
        error("Unexpected issue while collecting weeks")

    all_flux = Float64.(df.flux)

    if edges === nothing
        flux_min = minimum(all_flux)
        flux_max = maximum(all_flux)

        if isapprox(flux_min, flux_max)
            # avoid zero-width bins
            δ = max(abs(flux_min) * 0.05, 1e-6)
            flux_min -= δ
            flux_max += δ
        end

        edges = collect(range(flux_min, flux_max; length = nbins + 1))
    else
        edges = collect(edges)
    end

    issorted(edges) || error("Histogram edges must be sorted")
    length(edges) >= 2 || error("Histogram edges must contain at least two values")

    bin_centers = (edges[1:(end-1)] .+ edges[2:end]) ./ 2

    traces = GenericTrace[]

    # Compute counts per week first so we can size the axes cleanly
    counts_by_week = Dict{Int,Vector{Int}}()
    max_count = 0

    for week in weeks
        week_flux = Float64.(df[df.final_time .== week, :flux])
        h = fit(Histogram, week_flux, edges)
        counts = Int.(h.weights)
        counts_by_week[week] = counts
        if !isempty(counts)
            max_count = max(max_count, maximum(counts))
        end
    end

    xmin = first(edges)
    xmax = last(edges)
    zmin = minimum(weeks) - 0.5
    zmax = maximum(weeks) + 0.5

    push!(
        traces,
        xy_plane_trace(
            xmin,
            xmax,
            max_count,
            minimum(weeks) - 0.5;
            opacity = plane_opacity,
        ),
    )
    push!(traces, xz_plane_trace(xmin, xmax, zmin, zmax; opacity = plane_opacity))

    mode_x = Float64[]
    mode_y = Float64[]
    mode_z = Float64[]

    first_bar = true

    for week in weeks
        counts = counts_by_week[week]

        # mode bin
        max_bins = findall(==(maximum(counts)), counts)

        mode_center = if tie_method == :first
            bin_centers[first(max_bins)]
        elseif tie_method == :mean
            mean(bin_centers[max_bins])
        else
            error("Unsupported tie_method: $tie_method. Use :first or :mean.")
        end

        push!(mode_x, mode_center)
        push!(mode_y, 0.0)
        push!(mode_z, float(week))

        for bin_idx in eachindex(counts)
            count = counts[bin_idx]
            count == 0 && continue

            x0 = edges[bin_idx]
            x1 = edges[bin_idx+1]
            y0 = 0.0
            y1 = count
            z0 = week - z_thickness
            z1 = week + z_thickness

            push!(
                traces,
                cuboid_trace(
                    x0,
                    x1,
                    y0,
                    y1,
                    z0,
                    z1;
                    color = bar_color,
                    opacity = bar_opacity,
                    name = "Histogram bins",
                    showlegend = first_bar,
                ),
            )
            first_bar = false
        end
    end

    push!(
        traces,
        scatter3d(
            x = mode_x,
            y = mode_y,
            z = mode_z,
            mode = "lines+markers",
            name = "Mode trajectory",
            line = attr(width = 6),
            marker = attr(size = 5),
            hovertemplate = "Week: %{z}<br>Mode bin center: %{x:.4f}<extra></extra>",
        ),
    )

    layout = Layout(
        title = fig_title,
        scene = attr(
            xaxis = attr(title = "Flux"),
            yaxis = attr(title = "Bin count"),
            zaxis = attr(
                title = "Week",
                tickmode = "array",
                tickvals = weeks,
                ticktext = ["Week $w" for w in weeks],
            ),
            aspectmode = "manual",
            aspectratio = attr(x = 1.6, y = 1.0, z = 1.0),
            camera = attr(eye = attr(x = 1.7, y = 1.4, z = 1.1)),
        ),
        showlegend = true,
    )

    return PlotlyJS.plot(traces, layout)
end

#####################################################################
# END: BLOCKY HISTOGRAMS WITH A LINE BEHIND THEM                    #
#####################################################################

#####################################################################
# BEGIN: SIMPLER SCATTER VERSION                                    #
#####################################################################

function stacked_flux_histogram_steps_3d_colored(
    flux_df::DataFrame;
    edges = nothing,
    nbins::Int = 25,
    line_width::Real = 5,
    plane_opacity::Real = 0.12,
    tie_method::Symbol = :first,
    fig_title::AbstractString = "3D Flux Histogram Step Plots Across Weeks",
)
    df = dropmissing(flux_df, [:final_time, :flux])
    weeks = sort(unique(df.final_time))
    all_flux = Float64.(df.flux)

    # ----------------------------
    # Bin edges
    # ----------------------------
    if edges === nothing
        edges = collect(range(minimum(all_flux), maximum(all_flux); length = nbins + 1))
    else
        edges = collect(edges)
    end

    bin_centers = (edges[1:(end-1)] .+ edges[2:end]) ./ 2

    # ----------------------------
    # Assign colors per week
    # ----------------------------
    palette = [
        "#1f77b4",
        "#d62728",
        "#2ca02c",
        "#ff7f0e",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
    ]

    week_colors =
        Dict(week => palette[mod1(i, length(palette))] for (i, week) in enumerate(weeks))

    # ----------------------------
    # Compute histograms
    # ----------------------------
    counts_by_week = Dict{Int,Vector{Int}}()
    max_count = 0

    for week in weeks
        vals = Float64.(df[df.final_time .== week, :flux])
        h = fit(Histogram, vals, edges)
        counts = Int.(h.weights)
        counts_by_week[week] = counts
        max_count = max(max_count, maximum(counts))
    end

    xmin, xmax = first(edges), last(edges)
    zmin, zmax = minimum(weeks) - 0.5, maximum(weeks) + 0.5

    traces = GenericTrace[]

    # ----------------------------
    # Reference planes
    # ----------------------------
    push!(
        traces,
        surface(
            x = [xmin, xmax],
            y = [0.0, max_count],
            z = fill(zmin, 2, 2),
            colorscale = [[0.0, "white"], [1.0, "white"]],
            opacity = plane_opacity,
            showscale = false,
            hoverinfo = "skip",
            showlegend = false,
        ),
    )

    push!(
        traces,
        surface(
            x = [xmin, xmax],
            y = fill(0.0, 2, 2),
            z = [zmin zmin; zmax zmax],
            colorscale = [[0.0, "white"], [1.0, "white"]],
            opacity = plane_opacity,
            showscale = false,
            hoverinfo = "skip",
            showlegend = false,
        ),
    )

    # ----------------------------
    # Histogram step lines
    # ----------------------------
    mode_x = Float64[]
    mode_y = Float64[]
    mode_z = Float64[]

    for week in weeks
        counts = counts_by_week[week]

        # Step coordinates
        x_step = Float64[]
        y_step = Float64[]

        push!(x_step, edges[1]);
        push!(y_step, 0.0)

        for i in eachindex(counts)
            left = edges[i]
            right = edges[i+1]
            c = counts[i]

            push!(x_step, left);
            push!(y_step, c)
            push!(x_step, right);
            push!(y_step, c)
        end

        push!(x_step, edges[end]);
        push!(y_step, 0.0)

        z_step = fill(float(week), length(x_step))

        push!(
            traces,
            scatter3d(
                x = x_step,
                y = y_step,
                z = z_step,
                mode = "lines",
                name = "Week $week",
                line = attr(width = line_width, color = week_colors[week]),
                hovertemplate = "Week $week<br>Flux: %{x:.4f}<br>Count: %{y}<extra></extra>",
            ),
        )

        # Mode
        max_bins = findall(==(maximum(counts)), counts)
        mode_center =
            tie_method == :mean ? mean(bin_centers[max_bins]) : bin_centers[first(max_bins)]

        push!(mode_x, mode_center)
        push!(mode_y, 0.0)
        push!(mode_z, float(week))
    end

    # ----------------------------
    # Mode trajectory
    # ----------------------------
    push!(
        traces,
        scatter3d(
            x = mode_x,
            y = mode_y,
            z = mode_z,
            mode = "lines+markers",
            name = "Mode trajectory",
            line = attr(width = 6, color = "black"),
            marker = attr(size = 5),
        ),
    )

    # ----------------------------
    # Layout (white background)
    # ----------------------------
    layout = Layout(
        title = fig_title,
        scene = attr(
            xaxis = attr(
                title = "Flux",
                backgroundcolor = "white",
                gridcolor = "lightgray",
            ),
            yaxis = attr(
                title = "Bin count",
                backgroundcolor = "white",
                gridcolor = "lightgray",
            ),
            zaxis = attr(
                title = "Week",
                tickmode = "array",
                tickvals = weeks,
                ticktext = ["Week $w" for w in weeks],
                backgroundcolor = "white",
                gridcolor = "lightgray",
            ),
            aspectmode = "manual",
            aspectratio = attr(x = 1.6, y = 1.0, z = 1.0),
        ),
        paper_bgcolor = "white",
        plot_bgcolor = "white",
    )

    return plot(traces, layout)
end

#####################################################################
# END: SIMPLER SCATTER VERSION                                      #
#####################################################################

end

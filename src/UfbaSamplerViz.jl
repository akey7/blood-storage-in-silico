module UfbaSamplerViz

using DataFrames
using StatsBase
using Statistics
using PlotlyJS
using KernelDensity

export stacked_flux_histogram_steps_3d_colored, stacked_flux_kde_3d

"""
    stacked_flux_histogram_steps_3d_colored(
        flux_df::DataFrame;
        edges = nothing,
        nbins::Int = 25,
        line_width::Real = 5,
        plane_opacity::Real = 0.12,
        tie_method::Symbol = :first,
        fig_title::AbstractString = "3D Flux Histogram Step Plots Across Weeks",
    )

Makes a 3D plot of a given flux distribution trajectory through time. Histograms are plotted up the z axis and a line connecting the mode bins of all the histograms is plotted behind these histograms to show the trajectory of flux over time.

# Arguments
1. `flux_df::DataFrame`: Long format DataFrame of flux sampling results.
2. `edges = nothing`: Edges for bins, leave as `nothing` to accept default assignment.
3. `nbins::Int = 25`: Number of bins.
4. `line_width::Real = 5`: Line width for the plot.
5. `plane_opacity::Real = 0.12`: Opactiy behind the lines
6. `tie_method::Symbol = :first`: How to resolve ties in the mode bin selection. In practice, ties shouldn't happen, but it is here just in case.
7. `z_spacing::Real = 0.45`: How far apart each histogram is up the z axis.
8. `fig_title::AbstractString = "3D Flux Histogram Step Plots Across Weeks"`: Title for the figure.

# Returns
`PlotlyJS.Plot`

A PlotlyJS plot to display or save.
"""
function stacked_flux_histogram_steps_3d_colored(
    df::DataFrame;
    edges = nothing,
    nbins::Int = 25,
    line_width::Real = 5,
    plane_opacity::Real = 0.12,
    tie_method::Symbol = :first,
    z_spacing::Real = 0.45,
    fig_title::AbstractString = "3D Flux Histogram Step Plots Across Weeks",
)
    weeks = sort(unique(df.final_time))
    length(weeks) > 0 || error("No weeks found in flux_df")
    z_spacing > 0 || error("z_spacing must be positive")
    all_flux = Float64.(df.flux)
    if edges === nothing
        flux_min = minimum(all_flux)
        flux_max = maximum(all_flux)
        if isapprox(flux_min, flux_max)
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
    week_z = Dict(week => (i - 1) * z_spacing for (i, week) in enumerate(weeks))
    zvals = [week_z[w] for w in weeks]
    counts_by_week = Dict{eltype(weeks),Vector{Int}}()
    max_count = 0
    for week in weeks
        vals = Float64.(df[df.final_time .== week, :flux])
        h = fit(Histogram, vals, edges)
        counts = Int.(h.weights)
        counts_by_week[week] = counts
        max_count = max(max_count, maximum(counts))
    end
    xmin, xmax = first(edges), last(edges)
    zmin = minimum(zvals) - 0.15
    zmax = maximum(zvals) + 0.15
    traces = GenericTrace[]
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
            name = "XY plane",
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
            name = "XZ plane",
        ),
    )
    mode_x = Float64[]
    mode_y = Float64[]
    mode_z = Float64[]
    for week in weeks
        counts = counts_by_week[week]
        zpos = week_z[week]
        x_step = Float64[]
        y_step = Float64[]
        push!(x_step, edges[1])
        push!(y_step, 0.0)
        for i in eachindex(counts)
            left = edges[i]
            right = edges[i+1]
            c = counts[i]
            push!(x_step, left)
            push!(y_step, c)
            push!(x_step, right)
            push!(y_step, c)
        end
        push!(x_step, edges[end])
        push!(y_step, 0.0)
        z_step = fill(zpos, length(x_step))
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
        max_bins = findall(==(maximum(counts)), counts)
        mode_center =
            tie_method == :mean ? mean(bin_centers[max_bins]) :
            tie_method == :first ? bin_centers[first(max_bins)] :
            error("Unsupported tie_method: $tie_method. Use :first or :mean.")
        push!(mode_x, mode_center)
        push!(mode_y, 0.0)
        push!(mode_z, zpos)
    end
    push!(
        traces,
        scatter3d(
            x = mode_x,
            y = mode_y,
            z = mode_z,
            mode = "lines+markers",
            name = "Mode trajectory",
            line = attr(width = 6, color = "black"),
            marker = attr(size = 5, color = "black"),
            hovertemplate = "Week: %{text}<br>Mode bin center: %{x:.4f}<extra></extra>",
            text = string.(weeks),
        ),
    )
    layout = Layout(
        title = fig_title,
        paper_bgcolor = "white",
        plot_bgcolor = "white",
        scene = attr(
            xaxis = attr(
                title = "Flux",
                backgroundcolor = "white",
                gridcolor = "lightgray",
                zerolinecolor = "lightgray",
            ),
            yaxis = attr(
                title = "Bin count",
                backgroundcolor = "white",
                gridcolor = "lightgray",
                zerolinecolor = "lightgray",
            ),
            zaxis = attr(
                title = "Week",
                tickmode = "array",
                tickvals = zvals,
                ticktext = ["Week $w" for w in weeks],
                backgroundcolor = "white",
                gridcolor = "lightgray",
                zerolinecolor = "lightgray",
            ),
            aspectmode = "manual",
            aspectratio = attr(x = 1.6, y = 1.0, z = 0.7),
            camera = attr(eye = attr(x = 1.7, y = 1.4, z = 1.1)),
        ),
        showlegend = true,
    )
    return plot(traces, layout)
end

"""
    stacked_flux_kde_3d(
        df::DataFrame;
        z_spacing::Real = 0.35,
        npoints::Int = 256,
        bandwidth = nothing,
        line_width::Real = 5,
        peak_projection_width::Real = 2,
        plane_opacity::Real = 0.10,
        fig_title::AbstractString = "3D KDE Curves Across Weeks",
    )

Makes a 3D plot of a given flux distribution trajectory through time. Kernel density estimation are plotted up the z axis and a line connecting the mode bins of all the KDEs is plotted behind these histograms to show the trajectory of flux over time. Lines from the peak of each curve to the trajectory line make the plot more readable.

# Arguments
1. `df::DataFrame`: Long format DataFrame of flux sampling results.
2. `z_spacing::Real = 0.35`: How far apart each histogram is up the z axis.
3. `npoints::Int = 256`: Number of points on the KDE curve
4. `bandwidth = nothing`: Bandwidth of KDE, if specified.
5. `line_width::Real = 5`: Line width of the plot.
6. `peak_projection_width::Real = 2`: Width of line between peaks of curves and the mode trajectory line.
7. `plane_opacity::Real = 0.10`: Opacity of reference planes.
8. `fig_title::AbstractString = "3D KDE Curves Across Weeks"`: Figure title if specified.

# Returns
`PlotlyJS.Plot`

A PlotlyJS plot to display or save.
"""
function stacked_flux_kde_3d(
    df::DataFrame;
    z_spacing::Real = 0.35,
    npoints::Int = 256,
    bandwidth = nothing,
    line_width::Real = 5,
    peak_projection_width::Real = 2,
    plane_opacity::Real = 0.10,
    fig_title::AbstractString = "3D KDE Curves Across Weeks",
)
    weeks = sort(unique(df.final_time))
    length(weeks) > 0 || error("No weeks found in flux_df")
    z_spacing > 0 || error("z_spacing must be positive")
    npoints >= 32 || error("npoints should be at least 32 for a smooth KDE curve")
    all_flux = Float64.(df.flux)
    xmin = minimum(all_flux)
    xmax = maximum(all_flux)
    if isapprox(xmin, xmax)
        δ = max(abs(xmin) * 0.05, 1e-6)
        xmin -= δ
        xmax += δ
    end
    week_z = Dict(week => (i - 1) * z_spacing for (i, week) in enumerate(weeks))
    zvals = [week_z[w] for w in weeks]
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
    kdes = Dict{eltype(weeks),Any}()
    ymax = 0.0

    for week in weeks
        vals = Float64.(df[df.final_time .== week, :flux])

        kd = if bandwidth === nothing
            kde(vals; boundary = (xmin, xmax), npoints = npoints)
        else
            kde(vals; boundary = (xmin, xmax), npoints = npoints, bandwidth = bandwidth)
        end

        kdes[week] = kd
        ymax = max(ymax, maximum(kd.density))
    end
    zmin = minimum(zvals) - 0.12
    zmax = maximum(zvals) + 0.12
    traces = GenericTrace[]
    push!(
        traces,
        surface(
            x = [xmin, xmax],
            y = [0.0, ymax],
            z = fill(zmin, 2, 2),
            colorscale = [[0.0, "white"], [1.0, "white"]],
            opacity = plane_opacity,
            showscale = false,
            hoverinfo = "skip",
            showlegend = false,
            name = "XY plane",
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
            name = "XZ plane",
        ),
    )
    peak_x = Float64[]
    peak_y = Float64[]
    peak_z = Float64[]
    for week in weeks
        kd = kdes[week]
        zpos = week_z[week]
        x_curve = Float64.(kd.x)
        y_curve = Float64.(kd.density)
        z_curve = fill(zpos, length(x_curve))
        push!(
            traces,
            scatter3d(
                x = x_curve,
                y = y_curve,
                z = z_curve,
                mode = "lines",
                name = "Week $week",
                line = attr(width = line_width, color = week_colors[week]),
                hovertemplate = "Week $week<br>Flux: %{x:.4f}<br>Density: %{y:.4f}<extra></extra>",
            ),
        )
        peak_idx = argmax(y_curve)
        px = x_curve[peak_idx]
        py = y_curve[peak_idx]
        push!(peak_x, px)
        push!(peak_y, 0.0)
        push!(peak_z, zpos)
        push!(
            traces,
            scatter3d(
                x = [px, px],
                y = [0.0, py],
                z = [zpos, zpos],
                mode = "lines",
                name = "Peak projection",
                showlegend = false,
                line = attr(
                    width = peak_projection_width,
                    color = week_colors[week],
                    dash = "dot",
                ),
                hoverinfo = "skip",
            ),
        )
    end
    push!(
        traces,
        scatter3d(
            x = peak_x,
            y = peak_y,
            z = peak_z,
            mode = "lines+markers",
            name = "KDE peak trajectory",
            line = attr(width = 6, color = "black"),
            marker = attr(size = 5, color = "black"),
            hovertemplate = "Week: %{text}<br>Peak flux: %{x:.4f}<extra></extra>",
            text = string.(weeks),
        ),
    )
    layout = Layout(
        title = fig_title,
        paper_bgcolor = "white",
        plot_bgcolor = "white",
        scene = attr(
            xaxis = attr(
                title = "Flux",
                backgroundcolor = "white",
                gridcolor = "lightgray",
                zerolinecolor = "lightgray",
            ),
            yaxis = attr(
                title = "Density",
                backgroundcolor = "white",
                gridcolor = "lightgray",
                zerolinecolor = "lightgray",
            ),
            zaxis = attr(
                title = "Week",
                tickmode = "array",
                tickvals = zvals,
                ticktext = ["Week $w" for w in weeks],
                backgroundcolor = "white",
                gridcolor = "lightgray",
                zerolinecolor = "lightgray",
            ),
            aspectmode = "manual",
            aspectratio = attr(x = 1.6, y = 1.0, z = 0.7),
            camera = attr(eye = attr(x = 1.7, y = 1.4, z = 1.1)),
        ),
        showlegend = true,
    )
    return plot(traces, layout)
end

end

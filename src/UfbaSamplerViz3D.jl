module UfbaSamplerViz3D

using DataFrames
using StatsBase
using Statistics
using PlotlyJS

export stacked_flux_histogram_steps_3d_colored

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
7. `fig_title::AbstractString = "3D Flux Histogram Step Plots Across Weeks"`: Title for the figure.

# Returns
`PlotlyJS.Plot`

A PlotlyJS plot to display or save.
"""
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

    # Bin edges
    if isnothing(edges)
        edges = collect(range(minimum(all_flux), maximum(all_flux); length = nbins + 1))
    else
        edges = collect(edges)
    end

    bin_centers = (edges[1:(end-1)] .+ edges[2:end]) ./ 2

    # Assign colors per week
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

    # Compute histograms
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

    # Reference planes
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

    # Histogram step lines
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

    # Mode trajectory
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

    # Layout (white background)
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

end

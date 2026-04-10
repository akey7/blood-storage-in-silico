module UfbaSamplerAnalysisAndViz

using Base.Iterators
using Statistics
using StatsBase
using CairoMakie
using AlgebraOfGraphics
using ColorSchemes
using DataFrames
using DataFramesMeta
using ProgressMeter
using HypothesisTests
using MultipleTesting
using Chain
using ThreadsX
using Random
using EffectSizes
using CategoricalArrays
using MixedModels
using MixedModels: likelihoodratiotest
using GLM
using PlotlyJS

export histograms_for_reaction_v2,
    plot_all_histograms_for_reactions,
    diagnose_flux_stats,
    pivot_sampling_df_long,
    net_sink_fluxes,
    net_flux_from_up_and_down,
    calc_median_flux_df,
    combine_and_clean_addititve_final_time,
    prepare_median_flux_vector_matrix,
    prepare_measurements_and_sinks_report_df,
    compare_flux_distributions,
    global_mixed_model_test,
    pivot_sampling_df_long_cat,
    per_reaction_additive_time_test,
    reaction_additive_heatmap,
    reaction_additive_timecourse_heatmap_dfs,
    stacked_flux_histograms_3d

"""
    histograms_for_reaction_v2(long_sampling_df, reaction_id, reaction_string; bins = 20)

Plots histograms for a single reaction, with time points as separate panels and additives layered on top of each other in different colors. Draws a thick black dashed vertical line at the 0 point on all rows.

# Arguments
1. `long_sampling_df`: Sampling DataFrame, pivoted long
2. `reaction_id`: The reaction id for which the samples are being plotted.
3. `reaction_string`: The human-readable reaction string to place as a subtitle on the plot.
4. `subsystem`: Human-readable susbsytem of the reaction
5. `bins`: Number of bins in the histograms.

# Returns
`Figure`

Returns a Makie `Figure` to display or save.
"""
function histograms_for_reaction_v2(
    long_sampling_df,
    reaction_id,
    reaction_string,
    subsystem;
    bins = 20,
)
    plt_df = @chain long_sampling_df begin
        @rsubset(:reaction_id == reaction_id)
        @rtransform(:time_span = "Week $(:final_time - 1) to $(:final_time)")
    end
    clean_reaction_id = replace(reaction_id, "R_" => "")
    title = "$clean_reaction_id ($subsystem)\n$reaction_string"
    additive_palette = [
        "01-Ctrl AS3" => :dodgerblue,
        "02-Adenosine" => :orange,
        "03-Glutamine" => :blueviolet,
        "04-Methionine" => :crimson,
        "07-NAC" => :brown,
        "08-Taurine" => :magenta,
    ]
    hist_layer =
        data(plt_df) *
        mapping(:flux; color = :additive, row = :time_span => nonnumeric) *
        histogram(bins = bins) *
        visual(alpha = 0.5)
    zero_line_layer =
        data((flux = [0],)) *
        mapping(:flux) *
        visual(VLines; color = :black, linestyle = :dash, linewidth = 3)
    plt = hist_layer + zero_line_layer
    return draw(
        plt,
        scales(
            Color = (; palette = additive_palette),
            X = (; label = "Flux (mM/week)"),
            Y = (; label = "Sample Count"),
        );
        facet = (; linkxaxes = :all, linkyaxes = :all),
        figure = (; title = title, size = (700, 700)),
    )
end

"""
    pivot_sampling_df_long(sampling_df)

Pivots the sampling_df longer.

# Arguments
1. `sampling_df`: The sampling DataFrame in long format. The DataFrame should have a column for each reaction sampled, along with `:additive` and `:final_time` columns.

# Returns
`DataFrame`

Returns a DataFrame pivoted to long with the following columns:
1. `:additive`: The additive
2. `:final_time`: Final time
3. `:reaction_id`: The reaction id
4. `:flux`: The flux through that reaction at that sample.
"""
function pivot_sampling_df_long(sampling_df)
    long_sampling_df = stack(
        sampling_df,
        Not([:additive, :final_time]),
        variable_name = :reaction_id,
        value_name = :flux,
    )
    return long_sampling_df
end

"""
    plot_all_histograms_for_reactions(sampling_df, rxn_ids_to_strings; bins = 20)

Plots version 2 of all histograms (with time points for all additives on the same figure). This function saves each figure as they are made to the `output/uFBA_histograms_v2` folder. Displays a progress meter as the plots are made.

# Arguments
1. `sampling_df`: Wide DataFrame of uFBA sampling results.
2. `rxn_ids_to_strings`: Dictionary mapping reaction ids to human readable strings for plot subtitles.
3. `bins`: Number of bins to put onto histograms.
"""
function plot_all_histograms_for_reactions(sampling_df, rxn_ids_to_strings; bins = 20)
    if nrow(sampling_df) == 0
        @info "uFBA: Nothing to plot"
    else
        @info "uFBA: Plotting histograms, version 2"
        long_sampling_df = pivot_sampling_df_long(sampling_df)
        reaction_ids = unique(long_sampling_df.reaction_id)
        n_reaction_ids = length(reaction_ids)
        prog = Progress(n_reaction_ids, desc = "Writing histograms, version 2")
        for reaction_id in reaction_ids
            reaction_string = rxn_ids_to_strings[reaction_id]["rxn_string"]
            subsystem = rxn_ids_to_strings[reaction_id]["subsystem"]
            # reaction_name = rxn_ids_to_strings[reaction_id]["name"]
            fig = histograms_for_reaction_v2(
                long_sampling_df,
                reaction_id,
                reaction_string,
                subsystem;
                bins = bins,
            )
            filename = joinpath("output", "uFBA_histograms_v2", "$reaction_id.png")
            save(filename, fig)
            next!(prog)
        end
        finish!(prog)
    end
end

"""
    diagnose_flux_stats(sampling_df)

Calculates diagnostic statistics for uFBA models. Expecially useful for finding fluxes that average zero flux.

# Arugments
1. `sampling_df`: Wide format DataFrame of sampling results.

# Returns
`DataFrame`

Returns a DataFrame with the following columns
1. `:additive`: Additive of the model
2. `:final_time`: Final time of the regression of rates of metabolite concentration change.
3. `:reaction_id`: Id of the flux.
4. `:mean_is_approx_zero`: True if the mean is approximately zero, false if the mean is non-zero.
5. `:mean`: Mean of the flux distribution
6. `:std`: Sample standard deviation of the flux distribution.
"""
function diagnose_flux_stats(sampling_df)
    @info "Calculating flux descriptive statistics"
    long_sampling_df = pivot_sampling_df_long(sampling_df)
    descriptions_df = @chain long_sampling_df begin
        @groupby(:additive, :final_time, :reaction_id)
        @combine(
            :mean_is_approx_zero = isapprox(mean(:flux), 0.0, atol = 1.0e-10),
            :mean = mean(:flux),
            :std = std(:flux),
        )
    end
    return descriptions_df
end

"""
    calc_median_flux_df(sampling_df)

Calculate median fluxes for all reactions in the given wide sampling DataFrame.

# Arguments
1. `sampling_df`: Wide sampling DataFrame returned by [`sample_fluxes`](@ref BloodStorageInSilico.UfbaSampler.sample_fluxes)

# Returns
`DataFrame`

Returns a DataFrame with the following columns: `additive`, `final_time`, `reaction_id`, and `median_flux` columns.
"""
function calc_median_flux_df(sampling_df)
    long_df = pivot_sampling_df_long(sampling_df)
    median_df = @chain long_df begin
        @groupby(:additive, :final_time, :reaction_id)
        @combine(:median_flux = median(:flux))
    end
    return median_df
end

"""
    combine_and_clean_addititve_final_time(additive, final_time)

Combine additive and final time specifications into a single lowercase string with `-` and ` ` substituted with `_`.

# Arguments
1. `additive`: The additive string
2. `final_time`: The final time integer

# Returns
`String`

Returns a string formatted in the way specified above.
"""
function combine_and_clean_addititve_final_time(additive, final_time)
    cleaned_additive = @chain additive begin
        lowercase()
        replace("-" => "_", " " => "_")
    end
    combined = "$(cleaned_additive)_$(final_time)"
    return combined
end

"""
    prepare_median_flux_vector_matrix(sampling_df)

Prepare a data matrix of the uFBA results. Each row is a reaction, each column is an additive at a time point, and each element is the median flux for that row and column.

# Arguments
1. `sampling_df`: The wide formatted sampling DataFrame

# Returns
`DataFrame`

Returns a data matrix in the form of a DataFrame as specified above.
"""
function prepare_median_flux_vector_matrix(sampling_df)
    median_df = calc_median_flux_df(sampling_df)
    transformed_df = @chain median_df begin
        @rtransform(
            :additive_final_time =
                combine_and_clean_addititve_final_time(:additive, :final_time)
        )
        @select(:additive_final_time, :reaction_id, :median_flux)
        @orderby(:additive_final_time, :reaction_id)
        unstack(:reaction_id, :additive_final_time, :median_flux)
    end
    return transformed_df
end

"""
    prepare_measurements_and_sinks_report_df(absolute_quant_long_df, fba_model_metabolites_df, ufba_optimized_sinks_df, sampling_df)

Prepares two DataFrames detailing the number of metabolites in each model and how many sinks there are.

The first DataFrame is longer than the second and contains one row per metabolite. On each row, the additive, time, and metabolite are listed, along with boolean columns with whether the metabolite has up/down sinks, or is measured.

The second DataFrame aggregates these per-metabolite rows into per-model rows, outlining overal counts of unique metabolites, how many metabolites are measured, and the number of up and down sinks.

# Arguments
1. `absolute_quant_long_df`: The absolute quant measurements from `output/absolute_quant_long.csv`
2. `fba_model_metabolites_df`: The FBA model metabolites from `output/fba_model_metabolites.csv`
3. `ufba_optimized_sinks_df`: The sinks added to models during the uFBA process from `output/ufba_optimized_sinks_df.csv`
4. `sampling_df`: The wide sampling DataFrame from all uFBA runs.

# Returns
`Tuple{DataFrame,DataFrame}`

Returns two DataFrames, one for each report. The first element is the per-metabolite report, and the second element is the per-model report.

The first DataFrame contains the following columns
1. `additive`: Additive of the model
2. `final_time`: Final time of the model
3. `fba_metabolite_id`: Metabolite id from the FBA model
4. `is_measured`: `true` if the metabolite was measured
5. `has_sink`: `true` if there is a sink for the metabolite

The second DataFrame contains the following columns:
1. `additive`: Additive of the model
2. `final_time`: Final time of the model
3. `n_fba_metabolites`: Number of metabolites in the model
4. `n_measured`: The number of metabolites that have absolute quant approximations.
5. `n_with_sink`: The number of metabolites that have a sink.
6. `n_no_measure_no_sink`: The number of metabolites that have neither an absolute quant approximation nor any sinks.
"""
function prepare_measurements_and_sinks_report_df(
    absolute_quant_long_df,
    fba_model_metabolites_df,
    ufba_optimized_sinks_df,
    sampling_df,
)
    display(first(ufba_optimized_sinks_df, 5))
    long_sampling_df = pivot_sampling_df_long(sampling_df)
    additives = sort(unique(long_sampling_df.additive))
    final_times = sort(unique(long_sampling_df.final_time))
    fba_metabolite_ids = sort(unique(fba_model_metabolites_df.metabolite_id))
    tasks = product(fba_metabolite_ids, additives, final_times)
    n_tasks = length(tasks)
    prog = Progress(n_tasks, "Matching FBA, measured, and sink metabolite ids")
    rows = map(tasks) do t
        fba_metabolite_id, additive, final_time = t
        measured_df = @rsubset(
            absolute_quant_long_df,
            :Additive == additive,
            :Time == final_time,
            :Metabolite == fba_metabolite_id
        )
        sink_df = @rsubset(
            ufba_optimized_sinks_df,
            :additive == additive,
            :final_time == final_time,
            :metabolite_id == fba_metabolite_id,
        )
        is_measured = nrow(measured_df) > 0
        has_sink = nrow(sink_df) > 0
        next!(prog)
        return (
            additive = additive,
            final_time = final_time,
            fba_metabolite_id = fba_metabolite_id,
            is_measured = is_measured,
            has_sink = has_sink,
        )
    end
    unsorted_report_df = DataFrame(rows)
    report_df = @orderby(unsorted_report_df, :additive, :final_time, :fba_metabolite_id)
    report_by_model_df = @chain report_df begin
        @rtransform(:no_measure_no_sink = !:is_measured && !:has_sink)
        @groupby(:additive, :final_time)
        @combine(
            :n_fba_metabolites = length(unique(:fba_metabolite_id)),
            :n_measured = sum(:is_measured),
            :n_with_sink = sum(:has_sink),
            :n_no_measure_no_sink = sum(:no_measure_no_sink),
        )
        @orderby(:additive, :final_time)
    end
    return report_df, report_by_model_df
end

"""
    abs_maximum(xs)

Returns the SIGNED value with the maximum absolute value in the given vector. In other words, looks for the maximum magnitude while preserving the sign. A helper function for [`compare_flux_distributions`](@ref BloodStorageInSilico.UfbaSamplerAnalysisAndViz.compare_flux_distributions)

# Arguments
1. `xs`: The vector to search through.

# Returns
`Float64`

Returns the value that has the maximum magnitude while preserving the sign.
"""
abs_maximum(xs) = xs[argmax(abs.(xs))]

"""
    compare_flux_distributions(
        sampling_df;
        control_additive = "01-Ctrl AS3",
        n_samples = nothing,
        alpha = 0.01,
        interesting_cohen_effect_z = 2.0,
    )

For each (non-control) additive, time point, and reaction, compare all additives to the control to find statistical differences that point to interesting histograms and reactions to investigate.

# Arguments
1. `sampling_df`: The wide formatted sampling DataFrame
2. `control_additive = "01-Ctrl AS3"`: The name of the additive to use as the "control".
3. `n_samples = nothing`: If specified, number of samples without replacement to take from the control and treatment fluxes. The use of this is to reduce the power of the statistical tests, because with thousands of samples, most of the adjusted p-values tend to be significant.
4. `alpha = 0.01`: Either the adjusted p-value considered significant or `1.0 - alpha` is the confidence interval for the Cohen's effect measurement.
5. `interesting_cohen_effect_z = 2.0`: Z-scores for the Cohen's effect sizes are computed per reaction across all additives and time points. For an effect size to be considered interesting, its z-score must be greater than mor equal to this value.

# Returns
`NamedTuple`

Returns a tuple of two DataFrames:
1. `interesting_df`: DataFrame with interesting additives/time points/reactions. The most important columns in this DataFrame are `treatment_additive`, `final_time`, `reaction_id`, `all_interesting`. If `all_interesting` is `true`, that row might be worth a look!
2. `interesting_vs_uninteresting_df`: An aggregated report of the number of rows that are `all_interesting` or not. Shows if the statistical test thresholds are too permissive or too tight.
3. `score_ranking_df`: Ranking reactions by their most influential treatment additive and time point.
4. `effects_wide_df`: Standardized Cohen's effect sizes in a wide format for plotting in a heatmap. Ordered in descending order of the maximum effect size across all additives per each reaction.
5. `significance_wide_df`: Minimum t-test p-values across all additives per reaction in a wide format for plotting in a heatmap. Ordered the same way as the wide signficance DataFrame.
6. `heatmap_rank_df`: The DataFrame used to order the wide effects and significance DataFrames.
"""
function compare_flux_distributions(
    sampling_df;
    control_additive = "01-Ctrl AS3",
    n_samples = nothing,
    alpha = 0.01,
    interesting_cohen_effect_z = 2.0,
)
    Random.seed!(123)
    ci_quantile = 1.0 - alpha
    long_sampling_df = pivot_sampling_df_long(sampling_df)
    final_times = sort(unique(long_sampling_df.final_time))
    reaction_ids = sort(unique(long_sampling_df.reaction_id))
    treatments_df = @rsubset(long_sampling_df, :additive != control_additive)
    treatment_additives = sort(unique(treatments_df.additive))
    control_df = @rsubset(long_sampling_df, :additive == control_additive)
    tasks = product(treatment_additives, reaction_ids, final_times)
    n_tasks = length(tasks)
    println("n_tasks: $n_tasks")
    test_rows = ThreadsX.map(tasks) do t
        treatment_additive, reaction_id, final_time = t
        control_reaction_df =
            @rsubset(control_df, :final_time == final_time, :reaction_id == reaction_id)
        treatment_reaction_df = @rsubset(
            treatments_df,
            :additive == treatment_additive,
            :final_time == final_time,
            :reaction_id == reaction_id
        )
        control_fluxes_0 = control_reaction_df.flux
        treatment_fluxes_0 = treatment_reaction_df.flux
        control_fluxes =
            isnothing(n_samples) ? control_fluxes_0 :
            sample(control_fluxes_0, n_samples, replace = false)
        treatment_fluxes =
            isnothing(n_samples) ? treatment_fluxes_0 :
            sample(treatment_fluxes_0, n_samples, replace = false)
        t_test = UnequalVarianceTTest(treatment_fluxes, control_fluxes)
        t_test_p = pvalue(t_test)
        mw_test = MannWhitneyUTest(treatment_fluxes, control_fluxes)
        mw_p = pvalue(mw_test)
        cohen_d = CohenD(treatment_fluxes, control_fluxes; quantile = ci_quantile)
        cohen_effect = effectsize(cohen_d)
        # cohen_effect_size_ci = confint(cohen_d)
        unadjusted_row = (
            treatment_additive = treatment_additive,
            reaction_id = reaction_id,
            final_time = final_time,
            t_test_p = t_test_p,
            mw_p = mw_p,
            cohen_effect = cohen_effect,
            # cohen_effect_low = lower(cohen_effect_size_ci),
            # cohen_effect_high = upper(cohen_effect_size_ci),
        )
        print(".")
        return unadjusted_row
    end
    println("done")
    test_df = @chain test_rows begin
        DataFrame()
        @transform(
            :adj_t_test_p = adjust(:t_test_p, BenjaminiHochberg()),
            :adj_mw_p = adjust(:mw_p, BenjaminiHochberg())
        )
    end
    reaction_cohen_effect_z_df = @chain test_df begin
        @groupby(:reaction_id)
        @transform(:reaction_cohen_effect_z = zscore(:cohen_effect))
        @select(:treatment_additive, :final_time, :reaction_id, :reaction_cohen_effect_z)
    end
    interesting_df = @chain test_df begin
        leftjoin(
            reaction_cohen_effect_z_df;
            on = [:treatment_additive, :final_time, :reaction_id],
        )
        @rtransform(
            :t_test_significant = :adj_t_test_p <= alpha,
            :mw_significant = :adj_mw_p <= alpha,
            :large_effect = abs(:reaction_cohen_effect_z) >= interesting_cohen_effect_z
        )
        @rtransform(
            :all_interesting = :t_test_significant && :mw_significant && :large_effect
        )
        @orderby(:treatment_additive, :final_time, :all_interesting, :reaction_id)
    end
    interesting_vs_uninteresting_df = @chain interesting_df begin
        @groupby(:all_interesting)
        DataFrames.combine(nrow => :count)
    end
    log_p_max = 2.0
    score_ranking_df = @chain test_df begin
        leftjoin(
            reaction_cohen_effect_z_df;
            on = [:treatment_additive, :final_time, :reaction_id],
        )
        @rtransform(
            :score =
                abs(:reaction_cohen_effect_z) * min(-log10(:adj_t_test_p), log_p_max)
        )
        @groupby(:reaction_id)
        @combine @astable begin
            idx = argmax(:score)
            :treatment_additive = :treatment_additive[idx]
            :final_time = :final_time[idx]
            :max_score = :score[idx]
            :reaction_cohen_effect_z = :reaction_cohen_effect_z[idx]
            :adj_t_test_p = :adj_t_test_p[idx]
        end
        @orderby(-:max_score)
    end
    heatmap_rank_df = @chain score_ranking_df begin
        @groupby(:reaction_id)
        @combine(:max_max_score = maximum(:max_score))
        @orderby(-:max_max_score)
    end
    effects_wide_df = @chain interesting_df begin
        @select(:reaction_id, :treatment_additive, :reaction_cohen_effect_z)
        unstack(
            :reaction_id,
            :treatment_additive,
            :reaction_cohen_effect_z;
            combine = abs_maximum,
        )
        innerjoin(heatmap_rank_df, on = :reaction_id)
        @orderby(-:max_max_score)
        @select(Not(:max_max_score))
    end
    significance_wide_df = @chain interesting_df begin
        @select(:reaction_id, :treatment_additive, :adj_t_test_p)
        unstack(:reaction_id, :treatment_additive, :adj_t_test_p; combine = minimum)
        innerjoin(heatmap_rank_df, on = :reaction_id)
        @orderby(-:max_max_score)
        @select(Not(:max_max_score))
    end
    result = (
        interesting_df = interesting_df,
        interesting_vs_uninteresting_df = interesting_vs_uninteresting_df,
        score_ranking_df = score_ranking_df,
        effects_wide_df = effects_wide_df,
        significance_wide_df = significance_wide_df,
        heatmap_rank_df = heatmap_rank_df,
    )
    return result
end

"""
    pivot_sampling_df_long_cat(sampling_df)

Pivots the sampling DataFrame long, with a twist: It transforms `additive` and `final_time` into categorical variables.

# Arguments
1. `sampling_df`: Wide-format sampling DataFrame.

# Returns
`DataFrame`

Returns the wide DataFrame pivoted long, with `additive` transformed to `additive_cat` and `final_time` transformed to `final_time_cat`.
"""
function pivot_sampling_df_long_cat(sampling_df)
    long_cat_df = @chain sampling_df begin
        stack(Not([:additive, :final_time]), variable_name = :reaction_id, value_name = :flux)
        @transform begin
            :additive_cat = categorical(:additive; levels = sort(unique(:additive)))
            :final_time_cat = categorical(
                :final_time;
                levels = sort(unique(:final_time)),
                ordered = true,
            )
        end
        @select(:reaction_id, :additive_cat, :final_time_cat, :flux)
    end
    return long_cat_df
end

"""
    global_mixed_model_test(sampling_df)

Performs a global test to answer a simple question: Does the additive treatment and time make any statistical difference whatsoever in any reactions? Since the outcome of this test is "yes", as supported by visual inspection of flux distributions, this function just prints the results of the test instead of gathering the results into a neater data structure.

Note: Assumes that `01-Ctrl AS3` is the control group, and that in a sort of additive names, it will be placed first.

# Arguments
1. `sampling_df`: Wide-format sampling DataFrame.

# Returns
`NamedTuple`

Returns a named tuple with the results of the tests.
"""
function global_mixed_model_test(sampling_df)
    long_cat_df = pivot_sampling_df_long_cat(sampling_df)
    m_null = fit(MixedModel, @formula(flux ~ 1 + (1 | reaction_id)), long_cat_df)
    m_additive =
        fit(MixedModel, @formula(flux ~ additive_cat + (1 | reaction_id)), long_cat_df)
    m_time =
        fit(MixedModel, @formula(flux ~ final_time_cat + (1 | reaction_id)), long_cat_df)
    m_additive_time = fit(
        MixedModel,
        @formula(flux ~ additive_cat + final_time_cat + (1 | reaction_id)),
        long_cat_df,
    )
    m_full = fit(
        MixedModel,
        @formula(flux ~ additive_cat * final_time_cat + (1 | reaction_id)),
        long_cat_df,
    )
    println("========== FULL MODEL ==========")
    println(m_full)

    println("\n========== HYPOTHESIS TESTS ==========")

    println("\nMain effect of additive (controlling for time):")
    println(likelihoodratiotest(m_time, m_additive_time))

    println("\nMain effect of time (controlling for additive):")
    println(likelihoodratiotest(m_additive, m_additive_time))

    println("\nAdditive x time interaction:")
    println(likelihoodratiotest(m_additive_time, m_full))

    result = (
        long_cat_df = long_cat_df,
        m_null = m_null,
        m_additive = m_additive,
        m_time = m_time,
        m_additive_time = m_additive_time,
        m_full = m_full,
    )

    return result
end

"""
    per_reaction_additive_time_test(sampling_df)

For each reaction, this function answers the question: are there any times and additives that make any difference on the fluxes for each individual reaction? This function tests time alone, additive alone, and time interacting with additive.

# Arguments
1. `sampling_df`: Wide-format sampling DataFrame.

# Returns
`DataFrame`

Returns a DataFrame with a row per reaction and the results of F-tests and adjusted p-values for each reaction.
"""
function per_reaction_additive_time_test(sampling_df)
    long_cat_df = pivot_sampling_df_long_cat(sampling_df)
    reaction_ids = sort(unique(long_cat_df.reaction_id))
    n_reaction_ids = length(reaction_ids)
    prog = Progress(n_reaction_ids, "Per reaction additive/time test")
    rows = map(reaction_ids) do reaction_id
        reaction_df = @rsubset(long_cat_df, :reaction_id == reaction_id)
        m_time = lm(@formula(flux ~ final_time_cat), reaction_df)
        m_additive = lm(@formula(flux ~ additive_cat), reaction_df)
        m_time_additive =
            lm(@formula(flux ~ additive_cat + final_time_cat), reaction_df)
        m_time_additive_interaction =
            lm(@formula(flux ~ additive_cat * final_time_cat), reaction_df)
        additive_ftest = GLM.ftest(m_time.model, m_time_additive.model)
        time_ftest = GLM.ftest(m_additive.model, m_time_additive.model)
        interaction_ftest =
            GLM.ftest(m_time_additive.model, m_time_additive_interaction.model)
        row = (
            reaction_id = reaction_id,
            additive_fstat = additive_ftest.fstat[2],
            additive_p = additive_ftest.pval[2],
            time_fstat = time_ftest.fstat[2],
            time_p = time_ftest.pval[2],
            interaction_fstat = interaction_ftest.fstat[2],
            interaction_p = interaction_ftest.pval[2],
        )
        next!(prog)
        return row
    end
    unsorted_df = DataFrame(rows)
    sorted_and_adjusted_df = @chain unsorted_df begin
        @transform begin
            :additive_adj_p = adjust(:additive_p, BenjaminiHochberg())
            :time_adj_p = adjust(:time_p, BenjaminiHochberg())
            :interaction_adj_p = adjust(:interaction_p, BenjaminiHochberg())
        end
        @select(
            :reaction_id,
            :additive_fstat,
            :additive_adj_p,
            :time_fstat,
            :time_adj_p,
            :interaction_fstat,
            :interaction_adj_p
        )
        @orderby(:reaction_id)
    end
    return sorted_and_adjusted_df
end

"""
    reaction_additive_timecourse_heatmap_dfs(
        sampling_df;
        control_additive = "01-Ctrl AS3",
    )

Splits the samples per reaction and additives into pairs of control and treatment groups. Then it fits models that (1) test the effect of time only vs (2) the effects of additive and time. It then does an f-test for the statistical difference between the models to determine if additive has additional explanatory power beyond just time alone.

TODO: These tests are ridiculously overpowered. In order to prevent taking -log10(0.0), which many adjusted p-values are, the minimum adjusted p-value is clamped at `eps(Float64)` on the low end. This is higher than even the maximum adjusted p-values. This means that the `significance_value` for all reactions is fixed at approximately ~15. Perhaps thinning of the samples could be done in the future to fix this? Or sorting reactions not by significance but by order of magnitude of the F-statistic? I am keeping this here in case it is useful in the future, but am not generating the plot based on this in the current release.

# Arguments
1. `sampling_df`: Wide-format sampling DataFrame.
2. `control_additive`: The additive that is considered the "control" group for all the tests.

# Returns
`NamedTuple`

Returns a named tuple with data suitable for (1) diagnostics and (2) plotting with [`reaction_additive_heatmap`](@ref BloodStorageInSilico.UfbaSamplerAnalysisAndViz.reaction_additive_heatmap).

Available fields are:
1. `effects_wide_df`: The wide format of the F-tests for each test. Reactions on rows, additives on columns.
2. `significance_wide_df`: The wide format of `-log10.(max.(results_long_df.adj_p_value, eps(Float64)))`, with reactions on rows and additives on columns. See the TODO caveat above.
3. `results_long_df`: Long format of the results of all tests.
4. `rank_df`: DataFrame that controls the ranking of additives.
"""
function reaction_additive_timecourse_heatmap_dfs(
    sampling_df;
    control_additive = "01-Ctrl AS3",
)
    analysis_df = pivot_sampling_df_long(sampling_df)
    reactions = unique(analysis_df.reaction_id)
    additives = sort(unique(analysis_df.additive))
    treatment_additives = [a for a in additives if a != control_additive]
    isempty(treatment_additives) && throw(ArgumentError("No non-control additives found."))
    jobs = [
        (reaction_id, additive) for reaction_id in reactions for
        additive in treatment_additives
    ]
    n_jobs = length(jobs)
    println("n_jobs: $n_jobs")
    result_rows = ThreadsX.map(jobs) do (reaction_id, additive)
        pair_df = @chain analysis_df begin
            @rsubset(:reaction_id == reaction_id)
            @rsubset(:additive == control_additive || :additive == additive)
            @rtransform(:group = :additive == control_additive ? "control" : "treatment")
        end
        pair_df.group = categorical(pair_df.group)
        pair_df.final_time = categorical(pair_df.final_time)
        try
            reduced_model = lm(@formula(flux ~ final_time), pair_df)
            full_model =
                lm(@formula(flux ~ final_time + group + final_time & group), pair_df)
            ft = GLM.ftest(reduced_model.model, full_model.model)
            statistic = Float64(ft.fstat[2])
            p_value = Float64(ft.pval[2])
            print(".")
            return (
                reaction_id = reaction_id,
                additive = additive,
                statistic = statistic,
                p_value = p_value,
                n_obs = nrow(pair_df),
                model_ok = true,
            )
        catch err
            @warn "ANOVA fit failed for reaction/additive pair" reaction_id additive exception =
                (err, catch_backtrace())
            print(".")
            return (
                reaction_id = reaction_id,
                additive = additive,
                statistic = NaN,
                p_value = NaN,
                n_obs = nrow(pair_df),
                model_ok = false,
            )
        end
    end
    println("done")
    results_long_df = DataFrame(result_rows)
    results_long_df.adj_p_value = fill(NaN, nrow(results_long_df))
    valid_idx = findall(x -> !isnan(x), results_long_df.p_value)
    if !isempty(valid_idx)
        results_long_df.adj_p_value[valid_idx] =
            adjust(results_long_df.p_value[valid_idx], BenjaminiHochberg())
    end
    results_long_df.significance_value =
        -log10.(max.(results_long_df.adj_p_value, eps(Float64)))

    rank_df = @chain results_long_df begin
        @rsubset(:adj_p_value < 0.05)
        @groupby(:reaction_id)
        @combine(:sort_order = maximum(:significance_value))
        @orderby(-:sort_order)
    end

    effects_wide_df = @chain results_long_df begin
        @select(:reaction_id, :additive, :statistic)
        unstack(:reaction_id, :additive, :statistic)
    end
    significance_wide_df = @chain results_long_df begin
        @select(:reaction_id, :additive, :significance_value)
        unstack(:reaction_id, :additive, :significance_value)
    end
    for additive in treatment_additives
        if !(additive in names(effects_wide_df))
            effects_wide_df[!, additive] = fill(NaN, nrow(effects_wide_df))
        end
        if !(additive in names(significance_wide_df))
            significance_wide_df[!, additive] = fill(NaN, nrow(significance_wide_df))
        end
    end
    effects_wide_df = select(effects_wide_df, :reaction_id, treatment_additives...)
    significance_wide_df =
        select(significance_wide_df, :reaction_id, treatment_additives...)
    effects_wide_sorted_df = @chain effects_wide_df begin
        innerjoin(rank_df, on = :reaction_id)
        @orderby(-:sort_order)
        @select(Not(:sort_order))
    end
    significance_wide_sorted_df = @chain significance_wide_df begin
        innerjoin(rank_df, on = :reaction_id)
        @orderby(-:sort_order)
        @select(Not(:sort_order))
    end
    return (
        effects_wide_df = effects_wide_sorted_df,
        significance_wide_df = significance_wide_sorted_df,
        results_long_df = results_long_df,
        rank_df = rank_df,
    )
end

"""
    reaction_additive_heatmap(
        effects_result;
        top_n = 20,
        fig_size = (800, 800),
    )

Plots a pair of heatmaps side-by-side, one with effect sizes and the other with significance values. Meant to be useful for a variety of tests.

# Arguments
1. `effects_result`: A named tuple with at least two fields `effects_wide_df` (the effects taken over time) and `significance_wide_df` (significance of each effect test). Both DataFrames need reactions on the rows and additives on the columns, and the reactions should be ordered in some way and the same in both DataFrames.
2. `top_n = 20`: Limit the plot to the top n reactions. Defaults to 20.
3. `fig_size = (800, 800)`: Size of the figure, to accomodate total vertical height and a width for both heatmaps and their color legends.
4. `include_significance = false`: If true, includes the significance heatmap.
5. `effect_title = "Heatmap"`: Plot title for the effect heatmap.
6. `effect_colorbar_label = "Legend"`: Title for the colorbar legend.

# Returns
`Figure`

Returns a figure suitable for display or plotting.
"""
function reaction_additive_heatmap(
    effects_result;
    top_n = 20,
    fig_size = (800, 800),
    include_significance = false,
    effect_title = "Heatmap",
    effect_colorbar_label = "Legend",
)
    effects_wide_df = effects_result.effects_wide_df
    effects_plot_df = reverse(first(effects_wide_df, top_n))
    effects_row_labels = effects_plot_df.reaction_id
    effects_col_labels = names(effects_plot_df)[2:end]
    effects_heatmap_mat = Matrix(effects_plot_df[:, 2:end])
    effects_clims = (-maximum(abs, effects_heatmap_mat), maximum(abs, effects_heatmap_mat))
    fig = Figure(size = fig_size)
    effects_ax = Axis(
        fig[1, 1],
        title = effect_title,
        xticks = (1:length(effects_col_labels), effects_col_labels),
        yticks = (1:length(effects_row_labels), effects_row_labels),
        xticklabelrotation = π/4,
    )
    effects_hm = heatmap!(
        effects_ax,
        effects_heatmap_mat';
        colormap = Reverse(:RdBu_9),
        colorrange = effects_clims,
    )
    Colorbar(fig[1, 2], effects_hm; label = effect_colorbar_label, labelsize = 14)
    if include_significance
        significance_wide_df = effects_result.significance_wide_df
        significance_plot_df = reverse(first(significance_wide_df, top_n))
        significance_row_labels = significance_plot_df.reaction_id
        significance_col_labels = names(significance_plot_df)[2:end]
        significance_heatmap_mat = Matrix(significance_plot_df[:, 2:end])
        significance_clims = (
            -maximum(abs, significance_heatmap_mat),
            maximum(abs, significance_heatmap_mat),
        )
        significance_ax = Axis(
            fig[1, 3],
            title = "Significance",
            xticks = (1:length(significance_col_labels), significance_col_labels),
            yticks = (1:length(significance_row_labels), significance_row_labels),
            xticklabelrotation = π/4,
        )
        significance_hm = heatmap!(
            significance_ax,
            significance_heatmap_mat';
            colormap = :Blues,
            colorrange = significance_clims,
        )
        Colorbar(fig[1, 4], significance_hm; label = "Significance", labelsize = 14)
    end
    return fig
end

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

end

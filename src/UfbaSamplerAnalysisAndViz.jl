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
    reaction_additive_across_time_df,
    reaction_additive_across_time_heatmap,
    reaction_additive_timecourse_heatmap_dfs

"""
    histograms_for_reaction_v2(long_sampling_df, reaction_id, reaction_string; bins = 20)

Plots histograms for a single reaction, with time points as separate panels and additives layered on top of each other in different colors. Draws a thick black dashed vertical line at the 0 point on all rows.

# Arguments
1. `long_sampling_df`: Sampling DataFrame, pivoted long
2. `reaction_id`: The reaction id for which the samples are being plotted.
3. `reaction_string`: The human-readable reaction string to place as a subtitle on the plot.
4: `bins`: Number of bins in the histograms.

# Returns
`Figure`

Returns a Makie `Figure` to display or save.
"""
function histograms_for_reaction_v2(
    long_sampling_df,
    reaction_id,
    reaction_string;
    bins = 20,
)
    plt_df = @chain long_sampling_df begin
        @rsubset(:reaction_id == reaction_id)
        @rtransform(:time_span = "Week $(:final_time - 1) to $(:final_time)")
    end
    title = "$reaction_id\n$reaction_string"
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
            reaction_string = rxn_ids_to_strings[reaction_id]
            fig = histograms_for_reaction_v2(
                long_sampling_df,
                reaction_id,
                reaction_string;
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
3. `ranked_df`: Ranking reactions by their most influential treatment additive and time point.
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
    ranked_df = @chain test_df begin
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
    result = (
        interesting_df = interesting_df,
        interesting_vs_uninteresting_df = interesting_vs_uninteresting_df,
        ranked_df = ranked_df,
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

function reaction_additive_across_time_df(sampling_df; reference_additive = "01-Ctrl AS3")
    long_cat_df = pivot_sampling_df_long_cat(sampling_df)
    all_additives = sort(unique(long_cat_df.additive_cat))
    non_reference_additives =
        [additive for additive in all_additives if additive != reference_additive]
    reaction_ids = sort(unique(long_cat_df.reaction_id))
    reactions_additives = vec(collect(product(reaction_ids, non_reference_additives)))
    n_reactions_additives = length(reactions_additives)
    println("n_reactions_additives: $n_reactions_additives")
    effects_rows =
        ThreadsX.map(reactions_additives) do (reaction_id, non_reference_additive)
            comparison_additives = [reference_additive, non_reference_additive]
            sub_df = @chain long_cat_df begin
                @rsubset(:reaction_id == reaction_id, :additive_cat in comparison_additives)
                @select(:additive_cat, :final_time_cat, :flux)
            end
            model = lm(@formula(flux ~ additive_cat + final_time_cat), sub_df)
            ct = coeftable(model)
            coef_df = DataFrame(
                term = String.(ct.rownms),
                estimate = ct.cols[1],
                p_value = ct.cols[4],
            )
            additive_term_df = @rsubset(coef_df, occursin("additive_cat", :term))
            if nrow(additive_term_df) != 1
                throw(
                    ArgumentError(
                        "Expected exactly one additive coefficient for reaction=$reaction_id additive=$additive_id, found $(nrow(additive_term_df))",
                    ),
                )
            end
            estimate = additive_term_df.estimate[1]
            p_value = additive_term_df.p_value[1]
            comparison_row = (
                reaction_id = string(reaction_id),
                reference_additive = string(reference_additive),
                additive = string(non_reference_additive),
                estimate = estimate,
                p_value = p_value,
            )
            print(".")
            return comparison_row
        end
    println("done")
    effects_adj_df = @chain effects_rows begin
        DataFrame()
        @transform(:adj_p_value = adjust(:p_value, BenjaminiHochberg()))
        @transform(:neg_log1_p = -log1p.(:adj_p_value))

        # TODO: Sorting on -log10(p) + pivoting may not sort the reactions
        # properly after the pivot.

        @orderby(:neg_log1_p)
    end
    effects_wide_df = unstack(effects_adj_df, :reaction_id, :additive, :estimate)
    significance_wide_df = unstack(effects_adj_df, :reaction_id, :additive, :neg_log1_p)
    result = (
        effects_adj_df = effects_adj_df,
        effects_wide_df = effects_wide_df,
        significance_wide_df = significance_wide_df,
    )
    return result
end

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
    result_rows = map(jobs) do (reaction_id, additive)
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
        clamp.(-log10.(max.(results_long_df.adj_p_value, eps())), 0.0, 10.0)

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

function reaction_additive_across_time_heatmap(
    effects_result;
    top_n = 20,
    fig_size = (800, 800),
)
    effects_wide_df = effects_result.effects_wide_df
    effects_plot_df = first(effects_wide_df, top_n)
    effects_row_labels = effects_plot_df.reaction_id
    effects_col_labels = names(effects_plot_df)[2:end]
    effects_heatmap_mat = Matrix(effects_plot_df[:, 2:end])
    effects_clims = (-maximum(abs, effects_heatmap_mat), maximum(abs, effects_heatmap_mat))
    significance_wide_df = effects_result.significance_wide_df
    significance_plot_df = first(significance_wide_df, top_n)
    significance_row_labels = significance_plot_df.reaction_id
    significance_col_labels = names(significance_plot_df)[2:end]
    significance_heatmap_mat = Matrix(significance_plot_df[:, 2:end])
    significance_clims =
        (-maximum(abs, significance_heatmap_mat), maximum(abs, significance_heatmap_mat))
    fig = Figure(size = fig_size)
    effects_ax = Axis(
        fig[1, 1],
        title = "Effect Estimate",
        xticks = (1:length(effects_col_labels), effects_col_labels),
        yticks = (1:length(effects_row_labels), effects_row_labels),
        xticklabelrotation = π/4,
    )
    effects_hm = heatmap!(
        effects_ax,
        effects_heatmap_mat';
        colormap = :RdBu,
        colorrange = effects_clims,
    )
    Colorbar(fig[1, 2], effects_hm; label = "Estimate", labelsize = 14)
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
    return fig
end

end

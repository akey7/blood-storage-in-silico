# Guide for Reviewers

This document outlines reccomended steps for code reviewers interested in examining this code's inputs, quality control outputs, and outputs for downstream data visualization.

## Quick Start: Suggested Review Path

### Inspecting the Code

For a scientific review, inspect the code in this order:

1. [`NEW_TO_JULIA.md`](NEW_TO_JULIA.md): Enough Julia to understand the broad strokes of this repo if you are coming from a Python or R background.
2. [`INSTALLATION.md`](INSTALLATION.md): How to install the dependencies in this repo.
1. [`EXECUTION.md`](EXECUTION.md) to understand the workflow order.
2. `absolute_quant.jl` and `src/AbsoluteQuant.jl` to understand how metabolomics data are converted into concentration change rates.
3. `ufba_sampler.jl` and the following modules to understand how the uFBA sampling is performed:
   - `src/FbaModelBuilder.jl`
   - `src/MetaboliteBounds.jl`
   - `src/PruningOptimizations.jl`
   - `src/UfbaSampler.jl`
4. Downstream flux summaries and statistical comparisons: 
   - `ufba_sampler_analysis.jl` and `src/UfbaSamplerAnalysis.jl`
   - `ufba_sampler_analysis_2.jl` and `src/UfbaSamplerAnalysis2.jl`
5. Diagnostic outputs in `output/`, especially:
   - `ufba_sampling_status.csv`
   - `ufba_pruning_overview.csv`
   - `ufba_fba_breaks.csv`
   - `ufba_prune_breaks.csv`
   - `measurements_and_sinks_report.csv`
   - `measurements_and_sinks_report_by_model.csv`
6. Outputs the drive data visualization in `blood-storage-in-silico-viz` repo:
   - `ufba_sampling.csv`
   - `ufba_sampling_status.csv` (also in item 5 above)
   - `viz_k_means_pca_distance.xlsx`

### Input Data

Not all data could be committed to this repo because some of it is proprietary. Below is a list of all inputs and whether they are included in this repo. Please see an exhaustive list of inputs in [`INSTALLATION.md`](INSTALLATION.md).

## Overview

For scientific review, the most important questions are:

* What data entered the pipeline?
* How were model constraints constructed?
* Were condition-specific models feasible?
* How were flux samples generated?
* Were expected samples produced?
* How were results summarized?
* Which files directly generated the manuscript outputs?
* What assumptions or caveats affect biological interpretation?

## Computational workflow from inputs to outputs

At a high level, this repository implements a computational pipeline for modeling red blood cell metabolism under storage conditions using constraint-based modeling and flux sampling.

The exact details are described elsewhere in the repository and manuscript. Conceptually, the workflow follows this structure:

```text
experimental metabolomics data
        ↓
condition-specific metabolite constraints
        ↓
red blood cell metabolic model
        ↓
feasibility checks and constraint construction
        ↓
optimization steps
        ↓
flux sampling
        ↓
summary tables of sampled flux distributions
        ↓
statistical analysis and visualization
```

The central idea is that experimental metabolomics measurements are used to define condition-specific constraints on a red blood cell metabolic model. The constrained model is then analyzed computationally to estimate feasible flux distributions under each condition.

### Inputs

The main inputs generally include:

* A red blood cell metabolic model.
* Metabolomics-derived constraints.
* Condition labels, such as additive or treatment condition.
* Time-point information.
* Configuration files or tables that define analysis settings.

### Constraint construction

The repository constructs mathematical constraints that describe which flux states are feasible under each experimental condition.

These constraints may include:

* Reaction bounds.
* Metabolite accumulation or depletion constraints.
* Exchange, sink, or feasibility-restoring reactions.
* Condition-specific bounds derived from experimental data.
* Additional constraints required for optimization or sampling.

### Feasibility and optimization

Before flux sampling, the model must be feasible under the relevant constraints. Optimization steps are used to evaluate feasibility, define valid solution spaces, and apply project-specific objective functions.

Depending on the analysis stage, optimization may be used to:

* Confirm that a condition-specific model can satisfy its constraints.
* Minimize or penalize feasibility-restoring sink reactions.
* Generate constrained models suitable for sampling.
* Check solver status and numerical validity.

### Flux sampling

Flux sampling estimates a distribution of feasible flux states rather than a single optimal solution.

This is important because metabolic models often admit many possible flux configurations consistent with the same constraints. Sampling provides a way to characterize this feasible solution space.

The sampled fluxes are then summarized and compared across conditions and time points.

### Outputs

The major outputs generally include:

* Sampled flux tables.
* Summary statistics by reaction, condition, and time point.
* Quality-control tables.
* Intermediate files documenting model feasibility or constraint status.
* Tables used by downstream statistical analysis or figure-generation workflows.

### Downstream analysis

Downstream analyses may include:

* Comparing treatment conditions to the AS3 control.
* Ranking additives or treatments by distance from control.
* Identifying reactions with large treatment-associated flux changes.
* Performing dimensionality reduction such as PCA.
* Performing clustering or correlation analysis.
* Generating publication figures or diagnostic plots.

Some downstream visualization or figure-generation work may occur outside this Julia repository, particularly in Python or R.

---

## How to reproduce the main analysis

The main analysis is comprised of a seven-step workflow. Please see the [installation](INSTALLATION.md) and [execution](EXECUTION.md) documents on how to setup the project, ensure inputs are in place, and execute the modules of the workflow.

---

## Quality-control and validation checks

This repository includes, or should be interpreted alongside, quality-control checks that help establish whether the computational pipeline behaved as expected.

For scientific peer review, these checks are often more important than whether the code is written in a particular programming style. They help answer whether the computational outputs are traceable, internally consistent, and biologically interpretable.

### Feasibility checks

Constraint-based metabolic models must be feasible before their outputs can be interpreted.

Relevant checks may include:

* Whether the solver found a feasible or optimal solution.
* Whether infeasible models were detected and documented.
* Whether feasibility-restoring reactions were introduced.
* Whether solver termination statuses were inspected rather than ignored.

Solver status should be checked explicitly after optimization. A numerical result should not be interpreted as scientifically meaningful unless the corresponding optimization problem solved successfully.

### Sample count checks

Flux sampling should produce the expected number of samples for each condition and time point.

Useful checks include:

* Expected sample counts.
* Observed sample counts.
* Whether any condition/time-point combinations are missing.
* Whether all reactions expected in the output are present.
* Whether sampling failures are logged or reported.

For example, if the expected number of samples is computed from the number of chains, warmup/start variables, and collection iterations, the observed output should be compared against that expectation.

### Control-condition checks

The AS3 control condition can be used as an internal reference for several analyses.

Useful checks include:

* AS3 compared against itself should have zero distance in distance-based analyses.
* Control rows should be retained in ranking tables when useful as quality-control references.
* Each time point should have the expected AS3 control entry.
* Treatment comparisons should clearly identify which control condition was used.

Keeping the control condition visible in some output tables can make the analysis easier to audit.

### Random seed checks

Some parts of the pipeline may involve stochastic algorithms, especially sampling or clustering.

Relevant checks include:

* Whether random seeds are set.
* Whether seeds are documented.
* Whether repeated runs give consistent scientific conclusions.
* Whether stochastic algorithms are used only where appropriate.

Exact equality across stochastic runs may not always be expected, but results should be stable enough to support the scientific conclusions.

### Numerical checks

Optimization and sampling workflows can be sensitive to numerical tolerances.

Useful checks include:

* Solver termination status.
* Objective values.
* Constraint violations.
* Flux bounds.
* Unexpected `NaN`, `Inf`, or `missing` values.
* Extremely large or implausible flux values.
* Reactions with zero variance across all samples or conditions.

Numerical checks help distinguish biological signals from computational artifacts.

### Data integrity checks

Tabular data should be checked before and after major transformations.

Relevant checks include:

* Expected row counts.
* Expected column names.
* Expected condition labels.
* Expected time points.
* Missing values.
* Duplicate rows where uniqueness is expected.
* Successful joins between reaction identifiers and reaction annotations.
* Preservation of key identifiers such as `additive`, `final_time`, and `reaction_id`.

These checks are especially important when data are reshaped between long and wide formats.

### Statistical analysis checks

For statistical outputs, reviewers should verify that the analysis design matches the structure of the data.

Useful checks include:

* Whether comparisons are made within the correct time point.
* Whether AS3 is used consistently as the control.
* Whether multiple-testing correction is applied where appropriate.
* Whether effect sizes are reported alongside p-values.
* Whether scaling or centering is performed within the intended grouping structure.
* Whether missing values are handled explicitly.
* Whether PCA or clustering is fit separately or globally, as intended.

Statistical code should make the comparison structure explicit.

### Testing

Automated tests are useful for validating core computational behavior.

Tests may include:

* Constraint construction tests.
* Solver-status tests.
* Expected behavior for toy models.
* Data transformation tests.
* Regression tests for known outputs.
* Tests confirming that expected failures are handled gracefully.

Tests do not prove that the scientific model is correct, but they can provide evidence that the implementation behaves as intended.

### Reviewer-oriented interpretation

A reviewer inspecting this repository should look for whether the code:

* Clearly defines inputs and outputs.
* Checks feasibility before interpreting model results.
* Documents assumptions and caveats.
* Makes control comparisons explicit.
* Preserves identifiers through data transformations.
* Reports unexpected or failed cases.
* Separates exploratory code from manuscript-generating code.

These checks are central to evaluating the computational reliability of the project.

## Known caveats

This repository uses Julia for computational modeling and analysis. Julia may be less familiar to some biomedical reviewers than Python or R.

Known caveats include the following.

### Julia familiarity

Reviewers unfamiliar with Julia may need extra time to interpret syntax, package conventions, and project structure. Please see [the "New to Julia" guide.](NEW_TO_JULIA.md)

### Numerical optimization

Constraint-based modeling depends on numerical optimization. Solver behavior can be affected by tolerances, model formulation, and constraint scaling.

For this reason, solver status, feasibility, and numerical validity should be checked explicitly.

### Flux sampling

Flux sampling characterizes feasible regions of a constrained model. Sampling-based outputs should be interpreted as distributions over feasible flux states, not as direct measurements of intracellular flux.

Sampling may depend on algorithm settings, random seeds, warmup behavior, and model constraints.

Condition-specific differences in feasibility-restoring reactions should be documented and treated as caveats when interpreting downstream results.

### Model assumptions

Constraint-based models are abstractions of biological systems. They depend on assumptions about reaction inclusion, reversibility, bounds, compartments, measured metabolites, unmeasured metabolites, and objective functions.

The outputs should be interpreted in the context of those assumptions.

### Experimental constraints

Metabolomics-derived constraints are only as complete and reliable as the available experimental data and preprocessing steps. Missing metabolites, measurement uncertainty, and differences among studies can affect the feasible solution space.

### Downstream workflows

Some downstream figure generation or exploratory visualization may occur in Python or R. The Julia repository should be interpreted as the core modeling and sampling workflow unless otherwise documented.

### Reproducibility

Reproducibility depends on using the intended Julia environment, input data, configuration files, solver versions, and random seeds. The `Project.toml` and `Manifest.toml` files are important for recreating the computational environment.

# New to Julia?

## Purpose of this document

This document is intended for readers who are reviewing this repository but are not familiar with Julia. It is not a full Julia tutorial. Instead, it provides enough context to help scientific reviewers understand the structure, intent, and major computational patterns used in this project.

The main goals are to help reviewers:

1. Understand why Julia was used for this project.
2. Recognize common Julia syntax used throughout the repository.
3. Understand the major packages used in the computational workflow.

Reviewers do **not** need to become Julia programmers to evaluate the scientific logic of this repository. The goal is to make the codebase legible enough to assess provenance, reproducibility, and computational validity.

---

## Why this project uses Julia

This project uses Julia for the computationally intensive parts of the analysis, especially model construction, optimization, constraint handling, and flux sampling.

Julia was chosen because it combines several features that are useful for scientific computing:

* High-performance numerical execution without requiring most project code to be written in a lower-level language.
* Syntax that remains relatively close to mathematical notation.
* Strong support for mathematical optimization through packages such as `JuMP.jl`.
* Access to constraint-based modeling tools through packages such as `COBREXA.jl`.
* Native support for tabular data analysis through `DataFrames.jl` and `DataFramesMeta.jl`.
* Reproducible project environments through Julia’s `Project.toml` and `Manifest.toml` files.

The purpose of using Julia here is not to argue that Julia is inherently preferable to Python or R for all scientific work. Rather, Julia was selected because the core workflow involves repeated optimization, feasibility checks, and sampling over constrained metabolic models. These are areas where Julia’s scientific computing and optimization ecosystem is well suited.

---

## What parts of the repository generate manuscript results

<!-- TODO: Manually enter repository-specific content here. -->

---

## Minimal Julia syntax needed to read this code

This section summarizes the Julia syntax patterns most likely to appear in this repository. It is intended as a quick reference for reviewers who are already familiar with another programming language such as Python, R, MATLAB, or C/C++.

### Function definitions

Julia functions are commonly written using the `function ... end` syntax:

```julia
function my_function(x, y)
    result = x + y
    return result
end
```

Julia also supports keyword arguments. Arguments after a semicolon are keyword arguments:

```julia
function run_analysis(input_file; seed = 1234, verbose = true)
    # function body
end
```

In this example, `input_file` is a positional argument, while `seed` and `verbose` are keyword arguments with default values.

A function call using keyword arguments may look like this:

```julia
run_analysis("data/input.csv"; seed = 2026, verbose = false)
```

### `end` closes code blocks

Julia uses `end` to close functions, loops, conditionals, and other blocks.

```julia
if value > 0
    println("positive")
else
    println("not positive")
end
```

This is analogous to indentation in Python or braces in C-like languages.

### Comments

Single-line comments begin with `#`:

```julia
# This is a comment
```

Multiline comments use `#=` and `=#`:

```julia
#=
This is a multiline comment.
It can span multiple lines.
=#
```

### Variables and assignment

Variables are assigned using `=`:

```julia
x = 10
name = "AS3"
```

Julia variables do not require explicit type declarations in most routine code.

### Symbols

Julia code frequently uses symbols, especially as column names in DataFrames:

```julia
:additive
:final_time
:reaction_id
```

A symbol is not the same as a string, although symbols and strings are often used in related ways. In this repository, symbols most commonly appear as column identifiers.

For example:

```julia
df.additive
df.final_time
```

or:

```julia
select(df, :additive, :final_time)
```

### Broadcasting with dots

A dot before a function call means “apply this operation elementwise.”

```julia
log.(x)
```

This means “apply `log` to each element of `x`.”

Other examples:

```julia
x .+ 1
x .* y
abs.(values)
```

These are elementwise addition, elementwise multiplication, and elementwise absolute value.

This is conceptually similar to vectorized operations in R, NumPy, or pandas.

### Arrays and indexing

Julia arrays use square brackets:

```julia
values = [1, 2, 3, 4]
```

Julia uses **1-based indexing**, like R and MATLAB:

```julia
values[1]
```

returns the first element.

This differs from Python, which uses 0-based indexing.

### Ranges

A range can be written with `:`:

```julia
1:5
```

This represents the sequence `1, 2, 3, 4, 5`.

Ranges are often used in loops:

```julia
for i in 1:5
    println(i)
end
```

### Loops

A `for` loop looks like this:

```julia
for reaction_id in reaction_ids
    println(reaction_id)
end
```

A `while` loop looks like this:

```julia
while condition
    # repeated work
end
```

### Conditional logic

Conditional logic uses `if`, `elseif`, `else`, and `end`:

```julia
if final_time == 2
    println("early time point")
elseif final_time == 6
    println("late time point")
else
    println("intermediate time point")
end
```

### Mutating functions and the `!` convention

Julia functions that end in `!` usually modify one or more of their inputs in place.

Examples:

```julia
sort!(values)
insertcols!(df, :new_column => values)
```

The exclamation point is a naming convention rather than a language requirement. It signals to the reader that the function may change an existing object rather than returning a completely new one.

### `missing` and `NaN`

Julia distinguishes between `missing` and `NaN`.

`missing` represents an unknown or absent value:

```julia
x = missing
```

`NaN` means “not a number” and usually arises from invalid floating-point calculations:

```julia
0.0 / 0.0
```

For data analysis, this distinction matters. A `missing` value usually indicates absent data, while `NaN` usually indicates a numerical result that could not be computed.

### `nothing`

`nothing` is Julia’s equivalent of a deliberate “no value” result. It is conceptually similar to `None` in Python or `NULL` in R, although the exact semantics differ.

```julia
result = nothing
```

### DataFrame column access

Columns in a Julia `DataFrame` can often be accessed using dot syntax:

```julia
df.additive
df.final_time
```

They can also be referred to using symbols:

```julia
select(df, :additive, :final_time)
```

### Chained DataFrame workflows

This repository may use `@chain` to express a sequence of transformations:

```julia
@chain df begin
    @rsubset(:final_time == 3)
    @transform(:scaled_flux = (:flux .- mean(:flux)) ./ std(:flux))
end
```

This is similar in spirit to a tidyverse pipeline in R.

The input DataFrame is passed through each step in order. This style is used to make data transformations easier to read from top to bottom.

`@rtransform` is the row-wise counterpart to `@transform`: instead of evaluating its expression once as a vectorized operation over an entire column, it evaluates the expression once per row. Because `mean(:flux)` and `std(:flux)` need the whole column, they have to be computed beforehand; inside `@rtransform` the remaining arithmetic is written without broadcasting dots, since each evaluation only ever sees one row's values.

```julia
mean_flux = mean(df.flux)
std_flux = std(df.flux)

@chain df begin
    @rsubset(:final_time == 3)
    @rtransform(:scaled_flux = (:flux - mean_flux) / std_flux)
end
```

### DataFramesMeta syntax

`DataFramesMeta.jl` provides macros such as `@rsubset`, `@transform`, `@rtransform`, and `@combine`.

Common examples include:

```julia
@rsubset(df, :final_time == 3)
```

which keeps rows where `final_time == 3`.

```julia
@transform(df, :new_column = :old_column .* 2)
```

which creates a new column. `:old_column` refers to the whole column, so the elementwise multiplication needs the broadcasting dot in `.*`.

`@rtransform` does the same job row-by-row:

```julia
@rtransform(df, :new_column = :old_column * 2)
```

Here `:old_column` refers to a single row's value each time the expression is evaluated, so a plain `*` is enough — no broadcasting dot is needed.

```julia
@combine(groupby(df, :additive), :mean_flux = mean(:flux))
```

which summarizes grouped data.

The `:` before a column name indicates that the name refers to a DataFrame column.

### Multiple dispatch

Julia functions can have multiple definitions, called methods, depending on the types of their inputs.

For example, the same function name may be used for different kinds of objects:

```julia
process(x::DataFrame) = ...
process(x::Vector) = ...
```

For reviewers, the main implication is that understanding a function sometimes requires checking which type of object is being passed into it. In this repository, multiple dispatch may appear in code that handles models, constraints, reactions, or tabular data.

### Package imports

Julia packages are usually loaded with `using`:

```julia
using DataFrames
using CSV
using JuMP
```

This is broadly analogous to `import` in Python or `library()` in R.

### Reproducible environments

Julia projects commonly use two files to define their software environment:

* `Project.toml`
* `Manifest.toml`

`Project.toml` lists the main project dependencies.

`Manifest.toml` records the exact resolved versions of all dependencies. This makes it possible to recreate the same computational environment more precisely.

---

## Key packages and what they do

This repository uses several Julia packages. The table below, while not exhaustive, summarizes the major packages and their roles.

| Package                                                               | Role in this repository                                                                                              |
| --------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| `COBREXA.jl`                                                          | Constraint-based reconstruction and analysis. Used for working with metabolic models and constraint-based workflows. |
| `AbstractFBCModels.jl`                                                | Provides abstractions for flux-balance constraint models.                                                            |
| `SBMLFBCModels.jl`                                                    | Supports reading and working with SBML models that include flux-balance constraints.                                 |
| `ConstraintTrees.jl`                                                  | Supports structured representation of model constraints.                                                             |
| `JuMP.jl`                                                             | Mathematical optimization modeling language. Used to define optimization problems in a solver-independent way.       |
| `HiGHS.jl`                                                            | Interface to the HiGHS optimization solver. Used for linear and mixed-integer optimization problems.                 |
| `MathOptInterface.jl`                                                 | Common interface between JuMP and optimization solvers.                                                              |
| `DataFrames.jl`                                                       | Tabular data structures similar to pandas DataFrames or R data frames.                                               |
| `DataFramesMeta.jl`                                                   | Macro-based syntax for readable DataFrame transformations. Similar in spirit to tidyverse-style workflows.           |
| `CSV.jl`                                                              | Reading and writing CSV files.                                                                                       |
| `Chain.jl`                                                            | Provides the `@chain` macro for readable data transformation pipelines.                                              |
| `Statistics`                                                          | Standard Julia library for common statistical functions such as `mean` and `std`.                                    |
| `StatsBase.jl`                                                        | Additional statistical utilities.                                                                                    |
| `HypothesisTests.jl`                                                  | Statistical hypothesis tests.                                                                                        |
| `EffectSizes.jl`                                                      | Effect size calculations.                                                                                            |
| `MultivariateStats.jl`                                                | Multivariate methods such as PCA.                                                                                    |
| `Clustering.jl`                                                       | Clustering algorithms such as k-means.                                                                               |
| `CairoMakie.jl`, `AlgebraOfGraphics.jl`, and related plotting packages | Visualization and diagnostic plotting..                                                   |

Not every package is equally central to the final scientific outputs. Some packages are used for core model construction and optimization, while others support data wrangling, statistical analysis, plotting, testing, or development utilities.

---

## Notes for Python and R users

Many reviewers may be more familiar with Python or R than Julia. The table below gives rough analogies between Julia concepts and similar concepts in Python or R.

These analogies are meant to aid orientation. They are not exact equivalences.

| Julia concept                    | Rough Python analogy                          | Rough R analogy                                                |    |
| -------------------------------- | --------------------------------------------- | -------------------------------------------------------------- | -- |
| `DataFrame` from `DataFrames.jl` | `pandas.DataFrame`                            | `data.frame` or `tibble`                                       |    |
| `CSV.read`                       | `pandas.read_csv`                             | `readr::read_csv` or `data.table::fread`                       |    |
| `@chain`                         | pandas method chaining or pipe-like workflows | `%>%` or `                                                     | >` |
| `@rsubset`                       | row filtering with boolean masks              | `filter()`                                                     |    |
| `@transform`                     | `assign()` or column mutation                 | `mutate()`                                                     |    |
| `@rtransform`                    | row-wise `apply()` (no broadcasting needed)   | `rowwise()` + `mutate()`                                       |    |
| `@combine`                       | `groupby().agg()`                             | `summarize()`                                                  |    |
| `Project.toml`                   | project dependency metadata                   | `DESCRIPTION` or `renv` metadata                               |    |
| `Manifest.toml`                  | lockfile with exact dependency versions       | `renv.lock`                                                    |    |
| `missing`                        | `pd.NA` or missing values                     | `NA`                                                           |    |
| `NaN`                            | `numpy.nan`                                   | `NaN`                                                          |    |
| `nothing`                        | `None`                                        | `NULL`                                                         |    |
| `Dict`                           | `dict`                                        | named list or environment                                      |    |
| `Vector`                         | list or NumPy array                           | vector                                                         |    |
| `Matrix`                         | 2D NumPy array                                | matrix                                                         |    |
| `using PackageName`              | `import package`                              | `library(package)`                                             |    |
| functions ending in `!`          | in-place mutation by convention               | functions that modify an object by reference, where applicable |    |

### Important differences from Python

Julia uses 1-based indexing:

```julia
x[1]
```

returns the first element.

Python uses 0-based indexing:

```python
x[0]
```

returns the first element.

Julia also uses `end` to close blocks:

```julia
if x > 0
    println("positive")
end
```

Python uses indentation instead.

### Important differences from R

Julia distinguishes between `missing`, `NaN`, and `nothing`. These should not be treated as interchangeable.

Julia also has explicit package environments. The `Project.toml` and `Manifest.toml` files are important for reproducing the computational environment.

### DataFrame workflows

A Julia DataFrame workflow using `DataFramesMeta.jl` may look similar to tidyverse code:

```julia
@chain df begin
    @rsubset(:final_time == 3)
    @transform(:scaled_flux = (:flux .- mean(:flux)) ./ std(:flux))
end
```

This is broadly similar to:

```r
df |>
    filter(final_time == 3) |>
    mutate(scaled_flux = (flux - mean(flux)) / sd(flux))
```

or to a pandas workflow involving filtering and `assign`.

If you see `@rtransform` instead of `@transform` in the code, it's the row-wise version of the same idea — the expression is evaluated once per row rather than as a single vectorized operation over the whole column:

```julia
mean_flux = mean(df.flux)
std_flux = std(df.flux)

@chain df begin
    @rsubset(:final_time == 3)
    @rtransform(:scaled_flux = (:flux - mean_flux) / std_flux)
end
```

This is closer to `rowwise()` combined with `mutate()` in R, or a row-by-row `apply` in pandas, since neither `mean_flux` nor `std_flux` can be recomputed from inside a single row.

### Optimization workflows

The optimization code may be less familiar to readers coming from general Python or R data analysis.

Julia’s `JuMP.jl` package is a modeling language for optimization problems. It is used to define variables, constraints, and objectives before passing the problem to a solver such as HiGHS.

A simplified JuMP model may look like this:

```julia
model = Model(HiGHS.Optimizer)

@variable(model, x >= 0)
@objective(model, Min, x)
@constraint(model, x >= 1)

optimize!(model)
```

This code defines an optimization problem, solves it, and then allows the code to inspect the solution and solver status.

In this repository, similar concepts are applied to metabolic model constraints and flux variables.

---

## Code not used for manuscript results

Some code in a scientific repository may support exploration, debugging, diagnostics, or development. Such code can be useful, but it should be clearly distinguished from code that directly generates manuscript results.

When reviewing this repository, code should be interpreted according to its documented purpose.

Possible categories include:

| Category                   | Meaning                                                                                          |
| -------------------------- | ------------------------------------------------------------------------------------------------ |
| Manuscript-generating code | Code that directly produces tables, outputs, or intermediate files used for manuscript results.  |
| Supporting code            | Helper functions required by the manuscript-generating workflow.                                 |
| Quality-control code       | Code used to check feasibility, sampling completeness, numerical behavior, or data integrity.    |
| Diagnostic code            | Code used to inspect model behavior or debug unexpected results.                                 |
| Exploratory code           | Code used during development but not used to generate final results.                             |
| Deprecated code            | Code retained only for historical context or compatibility and not used in the current workflow. |

If exploratory, diagnostic, or deprecated code remains in the repository, it should not be interpreted as part of the final scientific analysis unless explicitly documented as such.

The manuscript-generating entry points should be listed elsewhere in this repository. Reviewers should prioritize those files when evaluating how the reported results were produced.

---

## Where to look next

Readers who want to inspect the repository in more detail should next consult the main [`README.md`](README.md) file, which will list many other places to look for documentation.

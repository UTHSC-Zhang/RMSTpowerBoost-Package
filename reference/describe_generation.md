# Summarize a generated dataset and its simulation mechanism

Given one element returned by
[`load_recipe_sets()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/load_recipe_sets.md)
(a list with `data` and a rich `meta` list), this builds tidy tables
describing the data-generating process: model family, baseline
parameters, linear predictor coefficients (intercept / treatment /
covariates), treatment assignment, and censoring.

## Usage

``` r
describe_generation(set)
```

## Arguments

- set:

  A single element from
  [`load_recipe_sets()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/load_recipe_sets.md),
  i.e. `list(data=<df>, meta=<list>)`.

## Value

A list of data frames:

- `header`: `n`, `L` (analysis horizon; stored as attr `"tau"`),
  `model`, `event_rate`, `achieved_censoring`.

- `baseline`: flattened baseline parameters.

- `effects`: coefficients for intercept/treatment and covariates (and
  formula terms if used).

- `treatment`: assignment mechanism and knobs.

- `censoring`: censoring mode/target/admin time.

- `covariates`: each generated covariate with its distribution and
  parameters.

- `files`: paths to on-disk files (csv/rds/rdata) when available.

## Examples

``` r
if (FALSE) { # \dontrun{
sets <- load_recipe_sets("checks/manifest.rds")
spec <- describe_generation(sets[[1]])
spec$header
spec$baseline
spec$effects
spec$treatment
spec$censoring
spec$covariates
spec$files
} # }
```

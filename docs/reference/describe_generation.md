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
covs <- list(covar_continuous("x", mean = 0, sd = 1))
rec <- recipe_quick_aft(
  n = 20, model = "aft_lognormal",
  baseline = list(mu = 2.7, sigma = 0.6),
  treat_effect = -0.2, covariates = covs,
  target_censoring = 0.2, seed = 123
)
out <- file.path(tempdir(), "rmst_describe_generation")
generate_recipe_sets(rec, out_dir = out, formats = "rds", n_reps = 1)
sets <- load_recipe_sets(file.path(out, "manifest.rds"))
spec <- describe_generation(sets[[1]])
spec$header
#>    n  L          model event_rate achieved_censoring
#> 1 20 NA AFT log-normal        0.7                0.3
spec$baseline
#>            baseline
#> 1 mu=2.7; sigma=0.6
spec$effects
#>          term coefficient
#> 1 (Intercept)         0.0
#> 2   treatment        -0.2
spec$treatment
#>      assignment allocation
#> 1 randomization        1:1
spec$censoring
#>             mode target admin_time
#> 1 target_overall    0.2        Inf
spec$covariates
#>   name       type   dist   parameters transform
#> 1    x continuous normal mean=0; sd=1      <NA>
spec$files
#>    csv
#> 1 <NA>
#>                                                                                      rds
#> 1 C:\\Users\\arnab\\AppData\\Local\\Temp\\RtmpcPUVbt/rmst_describe_generation/sc1_r1.rds
#>   rdata
#> 1  <NA>
```

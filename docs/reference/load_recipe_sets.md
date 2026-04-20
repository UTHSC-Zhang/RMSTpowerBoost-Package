# Load datasets from a recipe-sets manifest

Reads a `manifest.rds` created by
[`generate_recipe_sets()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/generate_recipe_sets.md),
loads one dataset per row (preferring `rds` â†’ `rdata` â†’ `csv` â†’
`txt`), restores attribute `"achieved_censoring"`, and returns a named
list of `list(data = <data.frame>, meta = <list>)`.

## Usage

``` r
load_recipe_sets(manifest_path)
```

## Arguments

- manifest_path:

  Path to `manifest.rds`.

## Value

A named list where each element is `list(data=..., meta=...)`.

## Examples

``` r
covs <- list(covar_continuous("x", mean = 0, sd = 1))
rec <- recipe_quick_aft(
  n = 20, model = "aft_lognormal",
  baseline = list(mu = 2.7, sigma = 0.6),
  treat_effect = -0.2, covariates = covs,
  target_censoring = 0.2, seed = 123
)
out <- file.path(tempdir(), "rmst_load_recipe_sets")
generate_recipe_sets(rec, out_dir = out, formats = "rds", n_reps = 1)
sets <- load_recipe_sets(file.path(out, "manifest.rds"))
names(sets)
#> [1] "sc001_r01"
str(sets[[1]]$meta)
#> List of 19
#>  $ dataset_id        : chr "sc001_r01"
#>  $ scenario_id       : int 1
#>  $ rep               : int 1
#>  $ seed_used         : int NA
#>  $ n                 : int 20
#>  $ n_treat           : int 5
#>  $ n_control         : int 15
#>  $ event_rate        : num 0.7
#>  $ achieved_censoring: num 0.3
#>  $ model             : chr "aft_lognormal"
#>  $ baseline          :List of 2
#>   ..$ mu   : num 2.7
#>   ..$ sigma: num 0.6
#>  $ effects           :List of 3
#>   ..$ intercept : num 0
#>   ..$ treatment : num -0.2
#>   ..$ covariates: NULL
#>  $ treatment         :List of 2
#>   ..$ assignment: chr "randomization"
#>   ..$ allocation: chr "1:1"
#>  $ censoring         :List of 3
#>   ..$ mode      : chr "target_overall"
#>   ..$ target    : num 0.2
#>   ..$ admin_time: num Inf
#>  $ covariates        :List of 1
#>   ..$ :List of 4
#>   .. ..$ name  : chr "x"
#>   .. ..$ type  : chr "continuous"
#>   .. ..$ dist  : chr "normal"
#>   .. ..$ params:List of 2
#>   .. .. ..$ mean: num 0
#>   .. .. ..$ sd  : num 1
#>  $ allocation        : chr "1:1"
#>  $ params            : list()
#>  $ files             :List of 4
#>   ..$ txt  : chr NA
#>   ..$ csv  : chr NA
#>   ..$ rds  : chr "C:\\Users\\arnab\\AppData\\Local\\Temp\\RtmpcPUVbt/rmst_load_recipe_sets/sc1_r1.rds"
#>   ..$ rdata: chr NA
#>  $ created_at        : chr "2026-04-19 22:23:20.839811"
```

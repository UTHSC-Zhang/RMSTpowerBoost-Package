# Simulate a dataset from a validated recipe (list-only)

Simulate a dataset from a validated recipe (list-only)

## Usage

``` r
simulate_from_recipe(recipe, seed = NULL)
```

## Arguments

- recipe:

  A validated recipe list (use
  [`validate_recipe()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/validate_recipe.md)).

- seed:

  Optional integer seed to override recipe\$seed.

## Value

A data.frame with columns `time`, `status`, `arm` (if treatment
present), plus covariates. Attribute `"achieved_censoring"` is attached.

## Examples

``` r
covs <- list(list(name="x", type="continuous", dist="normal", params=list(mean=0, sd=1)))
rec <- recipe_quick_aft(120, "aft_lognormal",
        baseline=list(mu=2.2, sigma=0.5), treat_effect=-0.2,
        covariates=covs, target_censoring=0.25)
dat <- simulate_from_recipe(rec, seed = 11)
```

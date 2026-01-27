# Validate a simulation recipe (list-only schema)

Checks/normalizes the simulation recipe and fills reasonable defaults.
This function does **not** use or require any analysis horizon (L/tau).

## Usage

``` r
validate_recipe(recipe)
```

## Arguments

- recipe:

  A named list defining n, covariates, treatment (optional), event_time
  (model, baseline, effects, optional frailty), and censoring.

## Value

A validated recipe list.

## Examples

``` r
r <- recipe_quick_aft(
  n = 100,
  model = "aft_lognormal",
  baseline = list(mu = 2.2, sigma = 0.5),
  treat_effect = -0.2,
  covariates = list(list(name="x", type="continuous", dist="normal",
                         params=list(mean=0, sd=1))),
  target_censoring = 0.25, allocation = "1:1"
)
r2 <- validate_recipe(r)
```

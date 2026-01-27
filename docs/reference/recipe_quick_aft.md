# Quick AFT recipe builder (list-only, no L/tau)

Quick AFT recipe builder (list-only, no L/tau)

## Usage

``` r
recipe_quick_aft(
  n,
  model = c("aft_lognormal", "aft_weibull"),
  baseline,
  treat_effect,
  covariates,
  target_censoring = 0.25,
  allocation = "1:1",
  seed = NULL
)
```

## Arguments

- n:

  Sample size.

- model:

  One of `"aft_lognormal"` or `"aft_weibull"`.

- baseline:

  Baseline parameter list (see model).

- treat_effect:

  Numeric treatment coefficient (on log-time scale).

- covariates:

  Covariate definitions (list of defs).

- target_censoring:

  Target overall censoring fraction (0-1).

- allocation:

  Allocation ratio string (e.g., "1:1").

- seed:

  Optional seed.

## Value

A recipe list suitable for
[`simulate_from_recipe`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/simulate_from_recipe.md).

## Examples

``` r
covs <- list(list(name="x", type="continuous", dist="normal", params=list(mean=0, sd=1)))
r <- recipe_quick_aft(120, "aft_lognormal",
       baseline=list(mu=2.3, sigma=0.5), treat_effect=-0.2,
       covariates=covs, target_censoring=0.25, allocation="1:1")
dat <- simulate_from_recipe(r, seed = 1)
```

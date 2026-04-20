# Quick PH recipe builder

Convenience wrapper that builds a recipe list for proportional-hazard
(PH) models, parallel to
[`recipe_quick_aft`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/recipe_quick_aft.md)
for AFT models. The returned recipe is validated and ready for
[`simulate_from_recipe`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/simulate_from_recipe.md).

## Usage

``` r
recipe_quick_ph(
  n,
  model = c("ph_exponential", "ph_weibull", "ph_gompertz", "cox_pwexp"),
  baseline,
  treat_effect,
  covariates = list(),
  target_censoring = 0.25,
  allocation = "1:1",
  seed = NULL
)
```

## Arguments

- n:

  Total sample size (integer).

- model:

  PH model. One of `"ph_exponential"`, `"ph_weibull"`, `"ph_gompertz"`,
  `"cox_pwexp"`.

- baseline:

  Named list of baseline hazard parameters:

  ph_exponential

  :   `list(rate = ...)`

  ph_weibull

  :   `list(shape = ..., scale = ...)`

  ph_gompertz

  :   `list(shape = ..., rate = ...)`

  cox_pwexp

  :   `list(rates = c(...), cuts = c(...))`

- treat_effect:

  Numeric log-hazard ratio for the treatment arm.

- covariates:

  List of covariate definitions. Use
  [`covar_continuous`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/covar_continuous.md),
  [`covar_binary`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/covar_binary.md),
  or
  [`covar_categorical`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/covar_categorical.md).
  Default [`list()`](https://rdrr.io/r/base/list.html).

- target_censoring:

  Target overall censoring fraction (0â€“1). Default `0.25`.

- allocation:

  Treatment allocation ratio string, e.g. `"1:1"` (default) or `"1:2"`.

- seed:

  Optional integer seed.

## Value

A recipe list suitable for
[`simulate_from_recipe`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/simulate_from_recipe.md).

## See also

[`recipe_quick_aft`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/recipe_quick_aft.md),
[`rmst.sim`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.sim.md),
[`simulate_from_recipe`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/simulate_from_recipe.md)

## Examples

``` r
r <- recipe_quick_ph(100, "ph_weibull",
       baseline = list(shape = 1.5, scale = 10),
       treat_effect = -0.5,
       covariates = list(covar_continuous("age")),
       target_censoring = 0.30)
dat <- simulate_from_recipe(r, seed = 1)
```

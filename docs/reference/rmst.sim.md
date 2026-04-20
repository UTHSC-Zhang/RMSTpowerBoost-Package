# Simulate survival data for RMST analysis

A unified, single-call wrapper for generating survival data suitable for
use as reference/pilot data in
[`rmst.power`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.power.md)
and
[`rmst.ss`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.ss.md).
Supports all seven built-in event-time models (AFT and PH families).

## Usage

``` r
rmst.sim(
  n,
  model = "aft_lognormal",
  baseline,
  treat_effect = 0,
  covariates = list(),
  target_censoring = 0.25,
  allocation = "1:1",
  L = NULL,
  seed = NULL
)
```

## Arguments

- n:

  Total sample size (split by `allocation` ratio).

- model:

  Event-time model. One of: `"aft_lognormal"` (default),
  `"aft_weibull"`, `"aft_loglogistic"`, `"ph_exponential"`,
  `"ph_weibull"`, `"ph_gompertz"`, `"cox_pwexp"`.

- baseline:

  Named list of baseline parameters (model-specific; see
  [`recipe_quick_aft`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/recipe_quick_aft.md)
  and
  [`recipe_quick_ph`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/recipe_quick_ph.md)).

- treat_effect:

  Numeric treatment coefficient. Log-time scale for AFT models;
  log-hazard ratio for PH models. Default `0`.

- covariates:

  List of covariate definitions. Elements can be created with
  [`covar_continuous`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/covar_continuous.md),
  [`covar_binary`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/covar_binary.md),
  or
  [`covar_categorical`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/covar_categorical.md),
  or supplied as raw named lists in the existing recipe format. Default
  [`list()`](https://rdrr.io/r/base/list.html) (no covariates).

- target_censoring:

  Target overall censoring fraction (0â€“1). Default `0.25`.

- allocation:

  Treatment allocation ratio string, e.g. `"1:1"` (default).

- L:

  Optional numeric truncation time. Stored as an attribute on the
  returned object for downstream use by
  [`rmst.power`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.power.md)
  /
  [`rmst.ss`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.ss.md).
  Does not affect data generation.

- seed:

  Optional integer seed for reproducibility.

## Value

A `data.frame` of class `c("rmst_simdata", "data.frame")` with columns
`time`, `status`, `arm` (when treatment is present), and one column per
covariate. Attributes:

- recipe:

  The validated recipe list used for generation.

- L:

  The truncation time if supplied, else `NULL`.

- achieved_censoring:

  Actual censoring fraction achieved.

## Details

Internally routes to
[`recipe_quick_aft`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/recipe_quick_aft.md)
for AFT models and
[`recipe_quick_ph`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/recipe_quick_ph.md)
for PH models, then calls
[`simulate_from_recipe`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/simulate_from_recipe.md).

## See also

[`rmst.power`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.power.md),
[`rmst.ss`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.ss.md),
[`recipe_quick_ph`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/recipe_quick_ph.md),
[`recipe_quick_aft`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/recipe_quick_aft.md),
[`covar_continuous`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/covar_continuous.md)

## Examples

``` r
df <- rmst.sim(
  n            = 150,
  model        = "aft_lognormal",
  baseline     = list(mu = 2.2, sigma = 0.5),
  treat_effect = -0.3,
  covariates   = list(covar_continuous("age"), covar_binary("female")),
  L            = 12,
  seed         = 42
)
print(df)
#> ── Simulated RMST Dataset ────────────────────────────────────
#>   Model          : AFT Log-Normal 
#>   N (total)      : 150  |  Arm 0: 81  |  Arm 1: 69
#>   Events         : 110 (73.3%)
#>   Censored       : 40 (26.7%)
#>   Truncation time: 12  (stored as attribute)
#> 
#> First 6 rows:
#>        time status arm        age female
#> 1  9.004186      1   0  1.3709584      0
#> 2 13.198736      1   0 -0.5646982      0
#> 3  6.817518      1   1  0.3631284      0
#> 4  9.655581      1   1  0.6328626      0
#> 5  3.960148      0   0  0.4042683      1
#> 6  4.942770      0   0 -0.1061245      1
s <- summary(df)
```

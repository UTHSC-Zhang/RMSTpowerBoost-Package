# Power analysis for RMST-based models via formula interface

Routes a formula-based call to the matching RMST power routine.

## Usage

``` r
rmst.power(
  formula,
  data,
  arm,
  sample_sizes,
  L,
  strata = NULL,
  strata_type = c("additive", "multiplicative"),
  dep_cens = FALSE,
  type = c("analytical", "boot"),
  alpha = 0.05,
  n_sim = 1000L,
  parallel.cores = 1L
)
```

## Arguments

- formula:

  A formula of the form `Surv(time, status) ~ cov1 + cov2`. Use `s()`
  wrapping (mgcv style) for smooth terms: `Surv(time, status) ~ s(age)`.

- data:

  A `data.frame` containing the reference (pilot) data.

- arm:

  Character string naming the treatment arm column (binary 0/1).

- sample_sizes:

  Integer vector of per-arm (or per-stratum) sample sizes to evaluate.

- L:

  Numeric truncation time for RMST.

- strata:

  Character column name, one-sided formula (`~col`), or `NULL`
  (default). Ignored when `dep_cens = TRUE`.

- strata_type:

  One of `"additive"` (default) or `"multiplicative"`. Only used when
  `strata` is non-`NULL` and `dep_cens = FALSE`.

- dep_cens:

  Logical; use dependent-censoring model? Default `FALSE`.

- type:

  One of `"analytical"` or `"boot"`. Auto-switched with a message when
  the requested type is unavailable for the chosen model.

- alpha:

  Significance level. Default `0.05`.

- n_sim:

  Number of bootstrap replicates (boot methods only). Default `1000`.

- parallel.cores:

  Number of cores for parallel processing. Default `1`.

## Value

An object of class `c("rmst_power", "list")` with elements
`results_data`, `results_plot`, `results_summary`, `model_output`, and
`.meta`.

## See also

[`rmst.ss`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.ss.md),
[`print.rmst_power`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst_power_methods.md),
[`summary.rmst_power`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst_power_methods.md),
[`plot.rmst_power`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst_power_methods.md)

## Examples

``` r
data(aft_lognormal_L12_n150, package = "RMSTpowerBoost")
r <- rmst.power(Surv(time, status) ~ age,
                data = aft_lognormal_L12_n150,
                arm  = "arm",
                sample_sizes = c(50, 100, 150),
                L    = 12)
#> --- Estimating parameters from pilot data for analytic calculation... ---
#> Model: Y_rmst ~ factor(arm) + age
#> --- Calculating asymptotic variance... ---
#> --- Calculating power for specified sample sizes... ---
print(r)
#> ── RMST Power Analysis ──────────────────────────
#>   Model           : Linear IPCW (Analytical) 
#>   Formula         : Surv(time, status) ~ age 
#>   Arm             : arm 
#>   Truncation time : 12 
#>   Alpha           : 0.05 
#> 
#> Power Results:
#>  N_per_Arm      Power
#>         50 0.08365565
#>        100 0.12692736
#>        150 0.16928696
#> 
#> Treatment Effect (from reference data):
#>   Estimand : RMST Difference 
#>   Estimate : -0.3007 
s <- summary(r)
plot(r)

```

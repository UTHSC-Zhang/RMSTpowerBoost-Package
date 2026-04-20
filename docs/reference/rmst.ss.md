# Sample size estimation for RMST-based models via formula interface

Routes a formula-based call to the matching RMST sample-size routine.

## Usage

``` r
rmst.ss(
  formula,
  data,
  arm,
  target_power,
  L,
  strata = NULL,
  strata_type = c("additive", "multiplicative"),
  dep_cens = FALSE,
  type = c("analytical", "boot"),
  alpha = 0.05,
  n_sim = 1000L,
  parallel.cores = 1L,
  n_start = 50L,
  n_step = 25L,
  max_n = 2000L,
  patience = 5L
)
```

## Arguments

- formula:

  A formula of the form `Surv(time, status) ~ cov1 + cov2`.

- data:

  A `data.frame` containing the reference (pilot) data.

- arm:

  Character string naming the treatment arm column (binary 0/1).

- target_power:

  Numeric target power (e.g., `0.80`).

- L:

  Numeric truncation time for RMST.

- strata:

  Character column name, one-sided formula (`~col`), or `NULL`.

- strata_type:

  One of `"additive"` (default) or `"multiplicative"`.

- dep_cens:

  Logical; use dependent-censoring model? Default `FALSE`.

- type:

  One of `"analytical"` or `"boot"`.

- alpha:

  Significance level. Default `0.05`.

- n_sim:

  Number of bootstrap replicates. Default `1000`.

- parallel.cores:

  Number of cores for parallel processing. Default `1`.

- n_start:

  Starting sample size for the search. Default `50`.

- n_step:

  Search increment. Default `25`.

- max_n:

  Maximum sample size to try. Default `2000`.

- patience:

  Number of consecutive non-improving steps before stopping. Default
  `5`.

## Value

An object of class `c("rmst_ss", "list")` with elements `results_data`,
`results_plot`, `results_summary`, `model_output`, and `.meta`.

## See also

[`rmst.power`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.power.md),
[`print.rmst_ss`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst_ss_methods.md),
[`summary.rmst_ss`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst_ss_methods.md),
[`plot.rmst_ss`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst_ss_methods.md)

## Examples

``` r
data(aft_lognormal_L12_n150, package = "RMSTpowerBoost")
r <- rmst.ss(Surv(time, status) ~ age,
             data         = aft_lognormal_L12_n150,
             arm          = "arm",
             target_power = 0.80,
             L            = 12)
#> --- Estimating parameters from pilot data for analytic search... ---
#> Model: Y_rmst ~ factor(arm) + age
#> --- Searching for Sample Size (Method: Analytic) ---
#>   N = 50/arm, Calculated Power = 0.084
#>   N = 75/arm, Calculated Power = 0.106
#>   N = 100/arm, Calculated Power = 0.127
#>   N = 125/arm, Calculated Power = 0.148
#>   N = 150/arm, Calculated Power = 0.169
#>   N = 175/arm, Calculated Power = 0.19
#>   N = 200/arm, Calculated Power = 0.211
#>   N = 225/arm, Calculated Power = 0.232
#>   N = 250/arm, Calculated Power = 0.253
#>   N = 275/arm, Calculated Power = 0.274
#>   N = 300/arm, Calculated Power = 0.294
#>   N = 325/arm, Calculated Power = 0.314
#>   N = 350/arm, Calculated Power = 0.334
#>   N = 375/arm, Calculated Power = 0.354
#>   N = 400/arm, Calculated Power = 0.374
#>   N = 425/arm, Calculated Power = 0.393
#>   N = 450/arm, Calculated Power = 0.412
#>   N = 475/arm, Calculated Power = 0.43
#>   N = 500/arm, Calculated Power = 0.449
#>   N = 525/arm, Calculated Power = 0.467
#>   N = 550/arm, Calculated Power = 0.484
#>   N = 575/arm, Calculated Power = 0.501
#>   N = 600/arm, Calculated Power = 0.518
#>   N = 625/arm, Calculated Power = 0.535
#>   N = 650/arm, Calculated Power = 0.551
#>   N = 675/arm, Calculated Power = 0.567
#>   N = 700/arm, Calculated Power = 0.582
#>   N = 725/arm, Calculated Power = 0.597
#>   N = 750/arm, Calculated Power = 0.611
#>   N = 775/arm, Calculated Power = 0.625
#>   N = 800/arm, Calculated Power = 0.639
#>   N = 825/arm, Calculated Power = 0.653
#>   N = 850/arm, Calculated Power = 0.666
#>   N = 875/arm, Calculated Power = 0.678
#>   N = 900/arm, Calculated Power = 0.69
#>   N = 925/arm, Calculated Power = 0.702
#>   N = 950/arm, Calculated Power = 0.714
#>   N = 975/arm, Calculated Power = 0.725
#>   N = 1000/arm, Calculated Power = 0.736
#>   N = 1025/arm, Calculated Power = 0.746
#>   N = 1050/arm, Calculated Power = 0.756
#>   N = 1075/arm, Calculated Power = 0.766
#>   N = 1100/arm, Calculated Power = 0.775
#>   N = 1125/arm, Calculated Power = 0.784
#>   N = 1150/arm, Calculated Power = 0.793
#>   N = 1175/arm, Calculated Power = 0.802
#> 
#> --- Calculation Summary ---
#> 
#> 
#> Table: Required Sample Size
#> 
#> | Target_Power| Required_N_per_Arm|
#> |------------:|------------------:|
#> |          0.8|               1175|
print(r)
#> ── RMST Sample Size Estimation ────────────────
#>   Model           : Linear IPCW (Analytical) 
#>   Formula         : Surv(time, status) ~ age 
#>   Arm             : arm 
#>   Truncation time : 12 
#>   Alpha           : 0.05 
#> 
#> Sample Size Result:
#>  Target_Power Required_N_per_Arm
#>           0.8               1175
#> 
#> Treatment Effect (from reference data):
#>   Estimand : RMST Difference 
#>   Estimate : -0.3007 
```

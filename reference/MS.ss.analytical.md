# Find Sample Size for a Multiplicative Stratified RMST Model (Analytic)

Calculates the required sample size for a target power using the
analytic (approximate) method from Wang et al. (2019).

## Usage

``` r
MS.ss.analytical(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  strata_var,
  target_power,
  linear_terms = NULL,
  L,
  alpha = 0.05,
  n_start = 50,
  n_step = 25,
  max_n_per_arm = 2000
)
```

## Arguments

- pilot_data:

  A `data.frame` containing pilot study data.

- time_var:

  A character string for the time-to-event variable.

- status_var:

  A character string for the event status variable (1=event,
  0=censored).

- arm_var:

  A character string for the treatment arm variable (1=treatment,
  0=control).

- strata_var:

  A character string for the stratification variable.

- target_power:

  A single numeric value for the desired power.

- linear_terms:

  An optional character vector of other covariate names.

- L:

  The numeric value for the RMST truncation time.

- alpha:

  The significance level (Type I error rate).

- n_start:

  The starting sample size *per stratum* for the search.

- n_step:

  The increment in sample size at each step of the search.

- max_n_per_arm:

  The maximum sample size *per stratum* to search up to.

## Value

A `list` containing:

- results_data:

  A `data.frame` with the target power and required sample size.

- results_plot:

  A `ggplot` object visualizing the search path.

- results_summary:

  A `data.frame` summarizing the estimated log(RMST Ratio).

## Details

This function performs an iterative search for the sample size required
to achieve a specified `target_power`. It uses the same underlying
theory and log-linear approximation as `MS.power.analytical`. It
performs a one-time estimation of the log-RMST ratio and its asymptotic
variance from the pilot data, then uses these parameters in an analytic
formula to efficiently find the required sample size.

## Examples

``` r
set.seed(456)
pilot_df_strat_effect <- data.frame(
 time = c(rexp(60, 0.15), rexp(60, 0.08)), # Effect
 status = rbinom(120, 1, 0.7),
 arm = rep(0:1, each = 60),
 region = factor(rep(c("A", "B"), each = 60))
)

ss_results <- MS.ss.analytical(
 pilot_data = pilot_df_strat_effect,
 time_var = "time", status_var = "status", arm_var = "arm", strata_var = "region",
 target_power = 0.80, L = 10,
 n_start = 100, n_step = 50
)
#> --- Estimating parameters from pilot data (log-linear approximation)... ---
#> --- Searching for Sample Size (Method: Analytic/Approximation) ---
#>   N = 100/stratum, Calculated Power = 1
#> 
#> --- Calculation Summary ---
#> 
#> 
#> Table: Required Sample Size
#> 
#> | Target_Power| Required_N_per_Stratum|
#> |------------:|----------------------:|
#> |          0.8|                    100|
print(ss_results$results_data)
#>   Target_Power Required_N_per_Stratum
#> 1          0.8                    100
```

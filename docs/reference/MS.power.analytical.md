# Analyze Power for a Multiplicative Stratified RMST Model (Analytic)

Performs power analysis for a multiplicative, stratified RMST model
using an analytic method based on the work of Wang et al. (2019).

## Usage

``` r
MS.power.analytical(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  strata_var,
  sample_sizes,
  linear_terms = NULL,
  L,
  alpha = 0.05
)
```

## Arguments

- pilot_data:

  A `data.frame` with pilot study data.

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

- sample_sizes:

  A numeric vector of sample sizes *per stratum* to calculate power for.

- linear_terms:

  An optional character vector of other covariate names.

- L:

  The numeric value for the RMST truncation time.

- alpha:

  The significance level (Type I error rate).

## Value

A `list` containing:

- results_data:

  A `data.frame` with sample sizes and corresponding powers.

- results_plot:

  A `ggplot` object visualizing the power curve.

## Details

This function is based on the method for evaluating center (stratum)
effects using a multiplicative model for RMST: \\\mu\_{ij} = \mu\_{0j}
\exp\\\beta'Z_i\\\\. The method uses IPCW with a stratified Cox model
for the censoring distribution.

Formal estimation of \\\beta\\ requires an iterative solver for the
estimating equation given in Equation (8) of Wang et al. (2019). Because
this is computationally intensive, this implementation uses a log-linear
working model fitted by weighted least squares to pseudo-observations
(`lm(log(Y_rmst) ~ ...)`) as a tractable approximation. The
approximation is consistent when the log-linear mean structure is
well-specified but may differ from the formal estimator under strong
misspecification.

The power calculation relies on the asymptotic variance of the log-RMST
ratio estimator, \\\hat{\tau}\\. The variance is derived from the robust
variance-covariance matrix of the `lm` fit, which serves as a proxy for
the formal sandwich estimator \\A^{-1}B(A^{-1})'\\ described in Theorem
1 of Wang et al. (2019).

## Examples

``` r
set.seed(123)
pilot_df_strat <- data.frame(
 time = rexp(120, 0.15),
 status = rbinom(120, 1, 0.6),
 arm = rep(0:1, each = 60),
 region = factor(rep(c("A", "B", "C"), each = 40))
)
pilot_df_strat$time[pilot_df_strat$arm == 1] <- pilot_df_strat$time[pilot_df_strat$arm == 1] * 1.5

power_results <- MS.power.analytical(
 pilot_data = pilot_df_strat,
 time_var = "time", status_var = "status", arm_var = "arm", strata_var = "region",
 sample_sizes = c(50, 75, 100),
 L = 10, alpha = 0.05
)
#> --- Estimating parameters from pilot data (log-linear approximation)... ---
#> Approximation Model: log(Y_rmst) ~ arm + region
#> --- Calculating power for specified sample sizes... ---
print(power_results$results_data)
#>   N_per_Stratum     Power
#> 1            50 0.4415527
#> 2            75 0.6027274
#> 3           100 0.7270471
```

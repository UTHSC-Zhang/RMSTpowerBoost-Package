# Analyze Power for a Multiplicative Stratified RMST Model via Simulation

Performs power analysis based on a multiplicative model for RMST for
stratified trials, using a bootstrap simulation approach.

## Usage

``` r
MS.power.boot(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  strata_var,
  sample_sizes,
  linear_terms = NULL,
  L,
  n_sim = 1000,
  alpha = 0.05,
  parallel.cores = 1
)
```

## Arguments

- pilot_data:

  A `data.frame` with pilot study data.

- time_var:

  A character string for the time-to-event variable.

- status_var:

  A character string for the event status variable.

- arm_var:

  A character string for the treatment arm variable.

- strata_var:

  A character string for the stratification variable.

- sample_sizes:

  A numeric vector of sample sizes *per stratum* to calculate power for.

- linear_terms:

  Optional character vector of covariates for the model.

- L:

  The numeric truncation time for RMST.

- n_sim:

  Number of bootstrap simulations.

- alpha:

  The significance level.

- parallel.cores:

  Number of cores for parallel processing.

## Value

A `list` containing:

- results_data:

  A `data.frame` of sample sizes and corresponding powers.

- results_plot:

  A `ggplot` object visualizing the power curve.

- results_summary:

  A `data.frame` with the estimated treatment effect (RMST Ratio).

## Details

This function estimates power through bootstrap simulation, resampling
within each stratum defined by `strata_var`. In each of the `n_sim`
iterations:

1.  Jackknife pseudo-observations for the RMST are calculated.

2.  A log-linear model ([`stats::lm`](https://rdrr.io/r/stats/lm.html))
    is fitted to the `log(pseudo_obs)`. This models the multiplicative
    relationship on the original RMST scale, i.e., \\\mu\_{ij} =
    \mu\_{0j} \exp\\\beta'Z_i\\\\. The model formula includes
    stratum-specific intercepts and interactions with the treatment arm.

3.  The p-value for the treatment effect is extracted from the model
    summary.

Power is determined as the proportion of simulations where the p-value
is less than `alpha`.

## Note

`status_var` should be `1`/`0`. `arm_var` should be `1`/`0`.
`strata_var` is a mandatory argument.

## Examples

``` r
if (FALSE) { # \dontrun{
pilot_df_strat <- data.frame(
 time = rexp(120, 0.15),
 status = rbinom(120, 1, 0.6),
 arm = rep(0:1, each = 60),
 region = factor(rep(c("A", "B", "C"), each = 40))
)
pilot_df_strat$time[pilot_df_strat$arm == 1] <- pilot_df_strat$time[pilot_df_strat$arm == 1] * 1.4

power_results <- MS.power.boot(
 pilot_data = pilot_df_strat,
 time_var = "time", status_var = "status", arm_var = "arm", strata_var = "region",
 sample_sizes = c(50, 75),
 L = 10,
 n_sim = 100 # Low n_sim for example
)
print(power_results$results_data)
} # }
```

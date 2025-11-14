# Estimate Sample Size for a Multiplicative Stratified RMST Model via Simulation

Performs sample size estimation based on a multiplicative model for RMST
for stratified trials, using iterative bootstrap simulations.

## Usage

``` r
MS.ss.boot(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  strata_var,
  target_power,
  linear_terms = NULL,
  L,
  n_sim = 1000,
  alpha = 0.05,
  parallel.cores = 1,
  patience = 5,
  n_start = 50,
  n_step = 25,
  max_n_per_arm = 2000
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

- target_power:

  A single numeric value for the target power (e.g., 0.80).

- linear_terms:

  Optional vector of covariates for the model.

- L:

  The numeric truncation time for RMST.

- n_sim:

  Number of bootstrap simulations per search step.

- alpha:

  The significance level.

- parallel.cores:

  Number of cores for parallel processing.

- patience:

  Number of consecutive non-improving steps in the search before
  terminating.

- n_start:

  Starting sample size per stratum for the search.

- n_step:

  Increment for the sample size search.

- max_n_per_arm:

  Maximum sample size per stratum to try.

## Value

A `list` containing:

- results_data:

  A `data.frame` with the target power and required N.

- results_plot:

  A `ggplot` object showing the search path.

- results_summary:

  A `data.frame` with the estimated treatment effect.

## Details

This function iteratively searches for the sample size required to
achieve a `target_power`. At each step of the search, it runs a full
bootstrap simulation (as described in `MS.power.boot`) to estimate the
power for the current sample size. The search proceeds until the target
power is met or other stopping criteria are satisfied. This process can
be very computationally intensive.

## Note

`status_var` should be `1`/`0`. `arm_var` should be `1`/`0`.
`strata_var` is a mandatory argument.

## Examples

``` r
if (FALSE) { # \dontrun{
pilot_df_strat_effect <- data.frame(
 time = c(rexp(60, 0.15), rexp(60, 0.08)), # Effect
 status = rbinom(120, 1, 0.7),
 arm = rep(0:1, each = 60),
 region = factor(rep(c("A", "B", "C"), each = 40))
)
ss_results <- MS.ss.boot(
 pilot_data = pilot_df_strat_effect,
 time_var = "time", status_var = "status", arm_var = "arm", strata_var = "region",
 target_power = 0.80, L = 10,
 n_sim = 100, # Low n_sim for example
 n_start = 100,
 n_step = 50, patience = 2
)
print(ss_results$results_data)
} # }
```

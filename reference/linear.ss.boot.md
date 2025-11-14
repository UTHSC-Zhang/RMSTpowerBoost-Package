# Find Sample Size for a Linear RMST Model via Simulation

Performs an iterative sample size search to achieve a target power based
on the direct linear regression model for RMST, using bootstrap
simulation.

## Usage

``` r
linear.ss.boot(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  target_power,
  linear_terms = NULL,
  L,
  n_sim = 1000,
  alpha = 0.05,
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

- target_power:

  A single numeric value for the target power (e.g., 0.80).

- linear_terms:

  Optional character vector of other covariates for the linear model.

- L:

  The numeric truncation time for RMST.

- n_sim:

  The number of bootstrap simulations per search step.

- alpha:

  The significance level.

- patience:

  The number of consecutive non-improving steps in the search before
  terminating.

- n_start:

  The starting sample size *per arm* for the search.

- n_step:

  The increment in sample size at each step of the search.

- max_n_per_arm:

  The maximum sample size *per arm* to search up to.

## Value

A `list` containing:

- results_data:

  A `data.frame` with the target power and the final required sample
  size per arm.

- results_plot:

  A `ggplot` object showing the search path.

- results_summary:

  A `data.frame` with summary statistics for the estimated treatment
  effect from the final simulation.

## Details

This function iteratively searches for the required sample size to
achieve a specified `target_power`. At each step of the search, it runs
a full bootstrap simulation (`n_sim` iterations), as described in
`linear.power.boot`, to estimate the power for the current sample size.
The search stops when the target power is achieved or other stopping
criteria (e.g., `patience`) are met. Due to the nested simulation
structure, this function can be very time-consuming.

## Note

`status_var` should be `1` for an event, `0` for censored. `arm_var`
should be `1` for treatment, `0` for control.

## Examples

``` r
if (FALSE) { # \dontrun{
pilot_df_effect <- data.frame(
  time = c(rexp(50, 0.1), rexp(50, 0.05)), # Effect present
  status = rbinom(100, 1, 0.8),
  arm = rep(0:1, each = 50)
)
ss_results <- linear.ss.boot(
  pilot_data = pilot_df_effect,
  time_var = "time",
  status_var = "status",
  arm_var = "arm",
  target_power = 0.80,
  L = 10,
  n_sim = 200, # Low n_sim for example
  patience = 2,
  n_start = 100,
  n_step = 50,
  max_n_per_arm = 500
)
print(ss_results$results_data)
print(ss_results$results_plot)
} # }
```

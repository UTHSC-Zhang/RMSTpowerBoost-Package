# Find Sample Size for a Semiparametric Additive RMST Model via Simulation

Performs an iterative sample size search to achieve a target power using
a flexible, semiparametric additive model for the RMST.

## Usage

``` r
GAM.ss.boot(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  strata_var = NULL,
  target_power,
  linear_terms = NULL,
  smooth_terms = NULL,
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

  An optional string for a stratification variable.

- target_power:

  A single numeric value for the target power.

- linear_terms:

  Optional character vector of covariates with a linear effect.

- smooth_terms:

  Optional character vector of covariates with a non-linear effect.

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

  The starting sample size per arm/stratum for the search.

- n_step:

  The increment in sample size at each step of the search.

- max_n_per_arm:

  The maximum sample size per arm/stratum to search up to.

## Value

A `list` containing:

- results_data:

  A `data.frame` with the target power and required sample size.

- results_plot:

  A `ggplot` object of the search path.

- results_summary:

  A `data.frame` summarizing the estimated treatment effect.

## Details

This function iteratively searches for the sample size required to
achieve a `target_power`. At each step, it runs a full bootstrap
simulation (as described in `GAM.power.boot`) to estimate the power for
the current sample size. The search proceeds until the target power is
met or other stopping criteria are satisfied.

## Note

This function's methodology is bootstrap-based.

## Examples

``` r
if (FALSE) { # \dontrun{
pilot_df_effect <- data.frame(
  time = c(stats::rexp(50, 0.1), stats::rexp(50, 0.04)), # Effect
  status = stats::rbinom(100, 1, 0.9),
  arm = rep(0:1, each = 50)
)

ss_results <- GAM.ss.boot(
  pilot_data = pilot_df_effect,
  time_var = "time",
  status_var = "status",
  arm_var = "arm",
  target_power = 0.80,
  L = 15,
  n_sim = 100,      # Low n_sim for example
  n_start = 100,
  n_step = 50,
  patience = 2,
  parallel.cores = 2
)
print(ss_results$results_data)
print(ss_results$results_plot)
} # }
```

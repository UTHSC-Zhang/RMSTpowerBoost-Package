# Estimate Sample Size for a Multiplicative Stratified RMST Model via Simulation

Performs sample size estimation based on a multiplicative model for RMST
for stratified trials, using iterative bootstrap simulations.

## Usage

``` r
MS.ss.boot.app(
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

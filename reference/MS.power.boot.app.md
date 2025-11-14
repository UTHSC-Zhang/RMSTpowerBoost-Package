# Analyze Power for a Multiplicative Stratified RMST Model via Simulation

Performs power analysis based on a multiplicative model for RMST for
stratified trials, using a bootstrap simulation approach. \#' @param
pilot_data A `data.frame` with pilot study data.

## Usage

``` r
MS.power.boot.app(
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

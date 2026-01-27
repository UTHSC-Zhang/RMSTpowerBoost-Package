# Analyze Power for a Linear RMST Model via Simulation

Performs a power analysis for given sample sizes based on the direct
linear regression model for RMST, using a bootstrap simulation approach.

## Usage

``` r
linear.power.boot.app(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  sample_sizes,
  linear_terms = NULL,
  L,
  n_sim = 1000,
  alpha = 0.05,
  parallel.cores
)
```

## Arguments

- sample_sizes:

  A numeric vector of sample sizes *per arm* to calculate power for.

## Value

A `list` containing results.

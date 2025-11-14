# Calculate Power for a Semiparametric Additive RMST Model via Simulation

Performs a power analysis for given sample sizes using a flexible,
semiparametric additive model for the RMST based on pseudo-observations.

## Usage

``` r
additive.power.boot.app(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  strata_var = NULL,
  sample_sizes,
  linear_terms = NULL,
  smooth_terms = NULL,
  L,
  n_sim = 1000,
  alpha = 0.05,
  parallel.cores = 1
)
```

## Arguments

- sample_sizes:

  A numeric vector of sample sizes per arm/stratum.

## Value

A list containing results.

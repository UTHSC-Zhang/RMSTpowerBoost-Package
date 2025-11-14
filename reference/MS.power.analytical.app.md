# Analyze Power for a Multiplicative Stratified RMST Model (Analytic)

Performs power analysis for a multiplicative, stratified RMST model
using an analytic method based on the work of Wang et al. (2019).

## Usage

``` r
MS.power.analytical.app(
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

- sample_sizes:

  A numeric vector of sample sizes *per stratum* to calculate power for.

- alpha:

  The significance level (Type I error rate).

## Value

A list containing results.

# Analyze Power for a Stratified Additive RMST Model (Analytic)

Performs power analysis for a stratified, additive RMST model.

## Usage

``` r
additive.power.analytical.app(
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

A list containing the results.

# Analyze Power for a Linear RMST Model (Analytic)

Performs power analysis using a direct formula based on the asymptotic
variance estimator for the linear RMST model.

## Usage

``` r
linear.power.analytical.app(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  sample_sizes,
  linear_terms = NULL,
  L,
  alpha = 0.05
)
```

## Arguments

- sample_sizes:

  A numeric vector of sample sizes *per arm* to calculate power for.

- alpha:

  The significance level for the power calculation (Type I error rate).

## Value

A `list` containing results.

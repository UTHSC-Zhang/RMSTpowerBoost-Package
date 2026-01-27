# Internal Helper to Estimate Multiplicative Stratified Model Parameters

This internal function contains the common logic for estimating the
log-RMST ratio and its variance from pilot data using a log-linear
approximation.

## Usage

``` r
.estimate_multiplicative_stratified_params(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  strata_var,
  linear_terms,
  L
)
```

## Value

A list containing `beta_effect`, `se_beta_n1`, and `n_strata`.

# Internal Helper to Estimate Additive Stratified Model Parameters

This internal function contains the common logic for estimating the
treatment effect and its variance from pilot data for the additive
stratified model. It is called by both the power and sample size
calculation functions to avoid code duplication.

## Usage

``` r
.estimate_additive_stratified_params(
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

A list containing `beta_effect` (the estimated treatment effect) and
`se_beta_n1` (the standard error for a sample size of 1).

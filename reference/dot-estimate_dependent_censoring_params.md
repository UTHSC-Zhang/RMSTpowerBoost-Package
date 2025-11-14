# Internal Helper to Estimate Covariate-Dependent Censoring Model Parameters

Internal helper for estimating the treatment effect and its variance
from pilot data under a single censoring mechanism that depends on
covariates.

This version uses one Cox model for censoring:
`Surv(time, status==0) ~ linear_terms` (treatment excluded by default),
builds IPCW \\w_i = 1/\hat G(Y_i\mid X_i)\\ at \\Y_i=\min(T_i,L)\\, fits
a weighted RMST regression, and computes a sandwich variance that
ignores uncertainty in \\\hat G\\.

## Usage

``` r
.estimate_dependent_censoring_params(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  dep_cens_status_var,
  linear_terms,
  L
)
```

## Value

A list with `beta_effect` (arm effect) and `se_beta_n1` (SE scaled to
N=1).

## Note

`dep_cens_status_var` is accepted for API compatibility but ignored
here.

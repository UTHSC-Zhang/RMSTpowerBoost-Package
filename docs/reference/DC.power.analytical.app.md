# Analyze Power for RMST with Covariate-Dependent Censoring (Analytic)

Analytic power for RMST regression when censoring depends on covariates
(single mechanism; no competing risks). Uses IPCW and a sandwich
variance.

## Usage

``` r
DC.power.analytical.app(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  dep_cens_status_var,
  sample_sizes,
  linear_terms = NULL,
  L,
  alpha = 0.05
)
```

## Arguments

- sample_sizes:

  Numeric vector of per-arm sample sizes for power.

- alpha:

  Two-sided Type I error (default 0.05).

## Value

A list: `results_data` (data.frame), `results_plot` (ggplot),
`results_summary` (data.frame with pilot effect).

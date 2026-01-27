# Analyze Power for RMST Model with Covariate-Dependent Censoring (Analytic)

Performs power analysis for an RMST model when the censoring mechanism
depends on observed covariates. Competing risks are not modeled here.

## Usage

``` r
DC.power.analytical(
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

- pilot_data:

  A `data.frame` with pilot data.

- time_var:

  Name of the time-to-event variable.

- status_var:

  Name of the primary event indicator (1=event, 0=censored).

- arm_var:

  Name of the treatment indicator (1=treatment, 0=control).

- sample_sizes:

  Numeric vector of per-arm sample sizes for power.

- linear_terms:

  Optional character vector of additional covariate names (used in both
  models).

- L:

  RMST truncation time.

- alpha:

  Two-sided Type I error.

## Value

A `list` with:

- results_data:

  A data.frame with `N_per_Arm` and `Power`.

- results_plot:

  A ggplot object of the power curve.

## Details

This function assumes a single censoring process whose hazard depends on
covariates (but not on treatment by default). It fits a Cox model for
censoring \$\$\Pr(\text{censoring by } t \mid X) = 1 - G(t \mid X)\$\$
using `Surv(time, status==0) ~ linear\_terms`, then forms
inverse-probability-of-censoring weights (IPCW) \\w_i = 1/\hat G(Y_i\mid
X_i)\\ evaluated at \\Y_i=\min(T_i,L)\\. The RMST regression \\E\[Y_i
\mid A_i,X_i\]\\ is then fit by weighted least squares, and power is
derived from a sandwich variance that ignores uncertainty from
estimating \\\hat G\\.

Note: `dep_cens_status_var` is accepted for API compatibility but is
ignored under this setting (no competing risks are modeled).

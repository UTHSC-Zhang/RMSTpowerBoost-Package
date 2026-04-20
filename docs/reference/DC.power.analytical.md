# Analyze Power for RMST Model with Covariate-Dependent Censoring (Analytic)

Performs power analysis for an RMST model when the censoring mechanism
depends on observed covariates.

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

  A `data.frame` with pilot study data.

- time_var:

  A character string for the time-to-event variable.

- status_var:

  A character string for the event status variable (1=event,
  0=censored).

- arm_var:

  A character string for the treatment arm variable (1=treatment,
  0=control).

- sample_sizes:

  A numeric vector of sample sizes *per arm* to calculate power for.

- linear_terms:

  Optional character vector of additional covariate names.

- L:

  The numeric value for the RMST truncation time.

- alpha:

  The significance level (Type I error rate).

## Value

A `list` with:

- results_data:

  A data.frame with `N_per_Arm` and `Power`.

- results_plot:

  A ggplot object of the power curve.

- results_summary:

  A data.frame summarizing the pilot arm effect.

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

Note: This implementation models a single censoring process and does not
handle competing risks.

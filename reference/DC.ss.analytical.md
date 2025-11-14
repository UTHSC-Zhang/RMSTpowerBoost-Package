# Find Sample Size for RMST with Covariate-Dependent Censoring (Analytic)

Iteratively finds required per-arm sample size for a target power, using
the same IPCW-based analytic variance as `DC.power.analytical`.

## Usage

``` r
DC.ss.analytical(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  target_power,
  linear_terms = NULL,
  L,
  alpha = 0.05,
  n_start = 50,
  n_step = 25,
  max_n_per_arm = 2000
)
```

## Arguments

- pilot_data:

  A `data.frame` containing pilot study data.

- time_var:

  Time variable.

- status_var:

  Event indicator (1=event, 0=censored).

- arm_var:

  Treatment indicator (1/0).

- target_power:

  Desired power.

- linear_terms:

  Optional covariates used in both models.

- L:

  RMST truncation time.

- alpha:

  Two-sided Type I error.

- n_start:

  Starting per-arm N.

- n_step:

  Step size for search.

- max_n_per_arm:

  Maximum per-arm N to search.

## Value

A `list` with:

- results_data:

  data.frame with `Target_Power` and `Required_N_per_Arm`.

- results_plot:

  ggplot showing the search path.

- results_summary:

  data.frame summarizing the pilot arm effect.

## Details

Uses a single censoring Cox model `Surv(time, status==0) ~ linear_terms`
to form IPCW and fits a weighted RMST regression. Treatment is excluded
from the censoring model by default. Competing risks are not modeled.
Variance ignores uncertainty in \\\hat G\\.

Note: `dep_cens_status_var` is accepted for API compatibility but
ignored here.

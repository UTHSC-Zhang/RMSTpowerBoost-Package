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

  A character string for the time-to-event variable.

- status_var:

  A character string for the event status variable (1=event,
  0=censored).

- arm_var:

  A character string for the treatment arm variable (1=treatment,
  0=control).

- target_power:

  A single numeric value for the desired power.

- linear_terms:

  Optional character vector of additional covariate names.

- L:

  The numeric value for the RMST truncation time.

- alpha:

  The significance level (Type I error rate).

- n_start:

  The starting sample size *per arm* for the search.

- n_step:

  The increment in sample size at each step of the search.

- max_n_per_arm:

  The maximum sample size *per arm* to search.

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

Note: This implementation models a single censoring process and does not
handle competing risks.

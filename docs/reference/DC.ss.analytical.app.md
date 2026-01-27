# Find Sample Size for RMST with Covariate-Dependent Censoring (Analytic)

Iterative per-arm N to reach target power using IPCW-based analytic
variance.

## Usage

``` r
DC.ss.analytical.app(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  dep_cens_status_var,
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

- target_power:

  Desired power.

- alpha:

  Two-sided Type I error (default 0.05).

- n_start:

  Starting per-arm N.

- n_step:

  Step size for search.

- max_n_per_arm:

  Maximum per-arm N to search.

## Value

A list: `results_data`, `results_plot`, `results_summary`.

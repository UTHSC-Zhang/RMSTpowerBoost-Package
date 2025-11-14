# Find Sample Size for a Linear RMST Model via Simulation

Performs an iterative sample size search to achieve a target power.

## Usage

``` r
linear.ss.boot.app(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  target_power,
  linear_terms = NULL,
  L,
  n_sim = 1000,
  alpha = 0.05,
  patience = 5,
  n_start = 50,
  n_step = 25,
  max_n_per_arm = 2000,
  parallel.cores
)
```

## Arguments

- target_power:

  A single numeric value for the target power (e.g., 0.80).

- patience:

  The number of consecutive non-improving steps before terminating.

- n_start:

  The starting sample size *per arm* for the search.

- n_step:

  The increment in sample size at each step of the search.

- max_n_per_arm:

  The maximum sample size *per arm* to search up to.

## Value

A `list` containing results.

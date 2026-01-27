# Find Sample Size for a Linear RMST Model (Analytic)

Calculates the required sample size for a target power using an analytic
formula based on the methods of Tian et al. (2014).

## Usage

``` r
linear.ss.analytical.app(
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

- target_power:

  A single numeric value for the desired power (e.g., 0.80 or 0.90).

- alpha:

  The significance level (Type I error rate).

- n_start:

  The starting sample size *per arm* for the search.

- n_step:

  The increment in sample size at each step of the search.

- max_n_per_arm:

  The maximum sample size *per arm* to search up to.

## Value

A `list` containing results.

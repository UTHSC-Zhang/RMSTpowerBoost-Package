# Find Sample Size for a Stratified Additive RMST Model (Analytic)

Calculates the required sample size for a target power.

## Usage

``` r
additive.ss.analytical.app(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  strata_var,
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

  A single numeric value for the desired power.

- alpha:

  The significance level (Type I error rate).

- n_start:

  The starting sample size *per stratum* for the search.

- n_step:

  The increment in sample size at each step of the search.

- max_n_per_arm:

  The maximum sample size *per stratum* to search up to.

## Value

A list containing the results.

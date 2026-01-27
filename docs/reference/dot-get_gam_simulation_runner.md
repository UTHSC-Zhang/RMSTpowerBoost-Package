# Internal Factory for Additive GAM RMST Simulation

This internal function prepares and returns a configured simulation
function.

## Usage

``` r
.get_gam_simulation_runner(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  strata_var,
  linear_terms,
  smooth_terms,
  L,
  alpha,
  n_sim,
  parallel.cores
)
```

## Value

A function that takes `n_per_group` and runs the simulation.

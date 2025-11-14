# Internal Factory for Linear IPCW RMST Simulation

This internal function prepares and returns a configured simulation
function.

## Usage

``` r
.get_linear_ipcw_simulation_runner(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  linear_terms,
  L,
  alpha,
  n_sim,
  parallel.cores
)
```

## Value

A function that takes `n_per_arm` and runs the simulation.

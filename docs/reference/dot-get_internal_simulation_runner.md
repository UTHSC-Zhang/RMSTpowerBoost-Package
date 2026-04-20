# Internal Factory for RMST Simulation Functions (Bootstrap)

This is a non-exported helper function that sets up and returns another
function responsible for running the core bootstrap simulation logic.
This factory pattern helps reduce code duplication between the
public-facing power and sample size functions. It prepares the model
formula and pilot data and returns the `run_power_sim` function which is
then called repeatedly.

## Usage

``` r
.get_internal_simulation_runner(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  strata_var,
  linear_terms,
  L,
  alpha,
  n_sim,
  parallel.cores = 10
)
```

## Value

A function that takes `n_per_stratum` and runs the simulation.

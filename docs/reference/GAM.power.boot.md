# Calculate Power for a Semiparametric Additive RMST Model via Simulation

Performs a power analysis for given sample sizes using a flexible,
semiparametric additive model for the RMST based on pseudo-observations.

## Usage

``` r
GAM.power.boot(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  strata_var = NULL,
  sample_sizes,
  linear_terms = NULL,
  smooth_terms = NULL,
  L,
  n_sim = 1000,
  alpha = 0.05,
  parallel.cores = 1
)
```

## Arguments

- pilot_data:

  A `data.frame` with pilot study data.

- time_var:

  A character string for the time-to-event variable.

- status_var:

  A character string for the event status variable.

- arm_var:

  A character string for the treatment arm variable.

- strata_var:

  An optional character string for a stratification variable.

- sample_sizes:

  A numeric vector of sample sizes per arm/stratum.

- linear_terms:

  Optional character vector of covariates with a linear effect.

- smooth_terms:

  Optional character vector of covariates with a non-linear effect.

- L:

  The numeric truncation time for RMST.

- n_sim:

  Number of bootstrap simulations.

- alpha:

  The significance level.

- parallel.cores:

  Number of cores for parallel processing.

## Value

A `list` containing:

- results_data:

  A `data.frame` of sample sizes and corresponding estimated power.

- results_plot:

  A `ggplot` object visualizing the power curve.

- results_summary:

  A `data.frame` with a summary of the estimated treatment effect.

## Details

This function estimates power by bootstrap resampling the pilot data,
computing RMST pseudo-observations, fitting a semiparametric additive
model, and recording whether the treatment effect is significant. If
`strata_var` is supplied, resampling is carried out within strata. The
model formula can include linear terms, non-linear smooth terms (`s()`),
and interactions.

Power is estimated as the proportion of simulations with a
treatment-effect p-value below `alpha`. This pseudo-observation approach
models RMST directly without using the IPCW estimating equations used
elsewhere in the package.

## Note

`status_var` should be `1` for an event, `0` for censored. `arm_var`
should be `1` for treatment, `0` for control.

## Examples

``` r
# \donttest{
pilot_df <- data.frame(
  time = rexp(100, 0.08),
  status = rbinom(100, 1, 0.7),
  arm = rep(0:1, each = 50),
  age = rnorm(100, 60, 10)
)
# Add a treatment effect
pilot_df$time[pilot_df$arm == 1] <- pilot_df$time[pilot_df$arm == 1] * 1.3

power_results <- GAM.power.boot(
  pilot_data = pilot_df,
  time_var = "time",
  status_var = "status",
  arm_var = "arm",
  sample_sizes = c(100, 150),
  linear_terms = "age",
  L = 15,
  n_sim = 100, # Use more sims in practice, e.g., 1000
  parallel.cores = 2
)
#> --- Calculating Power (Method: Additive GAM for RMST) ---
#> Model: pseudo_obs ~ arm + age
#> Simulating for n = 100 /arm ...
#> Simulating for n = 150 /arm ...
#> Total simulation time: 45.56 seconds
#> 
#> --- Simulation Summary ---
#> 
#> 
#> Table: Estimated Treatment Effect (RMST Difference)
#> 
#> |Statistic            |      Value|
#> |:--------------------|----------:|
#> |Mean RMST Difference | -1.6634763|
#> |Mean Standard Error  |  0.6459969|
#> |95% CI Lower         | -3.2260579|
#> |95% CI Upper         | -0.1008947|
print(power_results$results_data)
#>   N_per_Group Power
#> 1         100  0.64
#> 2         150  0.72
print(power_results$results_plot)

# }
```

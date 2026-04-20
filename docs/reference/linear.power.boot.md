# Analyze Power for a Linear RMST Model via Simulation

Performs a power analysis for given sample sizes based on the direct
linear regression model for RMST, using a bootstrap simulation approach.

## Usage

``` r
linear.power.boot(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  sample_sizes,
  linear_terms = NULL,
  L,
  n_sim = 1000,
  alpha = 0.05
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

- sample_sizes:

  A numeric vector of sample sizes *per arm* to calculate power for.

- linear_terms:

  Optional character vector of other covariates for the linear model.

- L:

  The numeric truncation time for RMST.

- n_sim:

  The number of bootstrap simulations to run for each sample size.

- alpha:

  The significance level (Type I error rate).

## Value

A `list` containing:

- results_data:

  A `data.frame` of sample sizes and corresponding estimated power.

- results_plot:

  A `ggplot` object visualizing the power curve.

- results_summary:

  A `data.frame` with summary statistics for the estimated treatment
  effect from the largest sample size simulation.

## Details

This function estimates power by generating a number of bootstrap
samples (`n_sim`) from the provided pilot data by resampling with
replacement. For each bootstrap sample, it performs the following steps:

1.  Estimates the censoring distribution using the Kaplan-Meier method
    ([`survival::survfit`](https://rdrr.io/pkg/survival/man/survfit.html)).

2.  Calculates Inverse Probability of Censoring Weights (IPCW) for each
    observation.

3.  Fits a weighted linear model
    ([`stats::lm`](https://rdrr.io/r/stats/lm.html)) to the RMST of the
    uncensored subjects.

4.  Extracts the p-value for the treatment `arm_var` coefficient.

The final power for a given sample size is the proportion of the `n_sim`
simulations where this p-value is less than the significance level
`alpha`. This simulation-based approach can be useful when analytic
approximations are less reliable, but it can be computationally
intensive.

## Note

`status_var` should be `1` for an event, `0` for censored. `arm_var`
should be `1` for treatment, `0` for control.

## Examples

``` r
pilot_df <- data.frame(
  time = rexp(100, 0.1),
  status = rbinom(100, 1, 0.7),
  arm = rep(0:1, each = 50),
  age = rnorm(100, 60, 8)
)
# Introduce a treatment effect for a more interesting example
pilot_df$time[pilot_df$arm == 1] <- pilot_df$time[pilot_df$arm == 1] * 1.5

power_results <- linear.power.boot(
  pilot_data = pilot_df,
  time_var = "time",
  status_var = "status",
  arm_var = "arm",
  linear_terms = "age",
  sample_sizes = c(100, 150, 200),
  L = 10,
  n_sim = 200 # Use more simulations in practice (e.g., 1000)
)
#> --- Calculating Power (Method: Linear RMST with IPCW) ---
#> Model: Y_rmst ~ arm + age
#> Simulating for n = 100 per arm...
#> Simulating for n = 150 per arm...
#> Simulating for n = 200 per arm...
#> Total simulation time: 2.03 seconds
#> 
#> --- Simulation Summary ---
#> 
#> 
#> Table: Estimated Treatment Effect (RMST Difference)
#> 
#> |Statistic            |     Value|
#> |:--------------------|---------:|
#> |Mean RMST Difference | 1.5554654|
#> |Mean Standard Error  | 0.3889167|
#> |95% CI Lower         | 0.8396978|
#> |95% CI Upper         | 2.2712330|
print(power_results$results_data)
#>   N_per_Arm Power
#> 1       100 0.830
#> 2       150 0.940
#> 3       200 0.985
print(power_results$results_plot)
```

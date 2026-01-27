# Find Sample Size for a Linear RMST Model (Analytic)

Calculates the required sample size for a target power using an analytic
formula based on the methods of Tian et al. (2014).

## Usage

``` r
linear.ss.analytical(
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

  A character string specifying the name of the time-to-event variable.

- status_var:

  A character string specifying the name of the event status variable
  (1=event, 0=censored).

- arm_var:

  A character string specifying the name of the treatment arm variable
  (1=treatment, 0=control).

- target_power:

  A single numeric value for the desired power (e.g., 0.80 or 0.90).

- linear_terms:

  An optional character vector of other covariate names to include in
  the model.

- L:

  The numeric value for the RMST truncation time.

- alpha:

  The significance level (Type I error rate).

- n_start:

  The starting sample size *per arm* for the search.

- n_step:

  The increment in sample size at each step of the search.

- max_n_per_arm:

  The maximum sample size *per arm* to search up to.

## Value

A `list` containing:

- results_data:

  A `data.frame` with the target power and the required sample size per
  arm.

- results_plot:

  A `ggplot` object visualizing the sample size search path.

- results_summary:

  A `data.frame` summarizing the treatment effect from the pilot data
  used for the calculation.

## Details

This function performs an iterative search to find the sample size
needed to achieve a specified `target_power`. It uses the same
underlying theory as `linear.power.analytical`. First, it estimates the
treatment effect size and its asymptotic variance from the pilot data.
Then, it iteratively calculates the power for increasing sample sizes
using the analytic formula until the target power is achieved.

## Examples

``` r
pilot_df <- data.frame(
  time = c(rexp(50, 0.1), rexp(50, 0.07)), # Introduce an effect
  status = rbinom(100, 1, 0.8),
  arm = rep(0:1, each = 50),
  age = rnorm(100, 55, 10)
)
ss_results <- linear.ss.analytical(
  pilot_data = pilot_df,
  time_var = "time",
  status_var = "status",
  arm_var = "arm",
  target_power = 0.80,
  L = 10
)
#> --- Estimating parameters from pilot data for analytic search... ---
#> Model: Y_rmst ~ factor(arm)
#> --- Searching for Sample Size (Method: Analytic) ---
#>   N = 50/arm, Calculated Power = 0.11
#>   N = 75/arm, Calculated Power = 0.144
#>   N = 100/arm, Calculated Power = 0.178
#>   N = 125/arm, Calculated Power = 0.212
#>   N = 150/arm, Calculated Power = 0.245
#>   N = 175/arm, Calculated Power = 0.278
#>   N = 200/arm, Calculated Power = 0.311
#>   N = 225/arm, Calculated Power = 0.343
#>   N = 250/arm, Calculated Power = 0.375
#>   N = 275/arm, Calculated Power = 0.405
#>   N = 300/arm, Calculated Power = 0.435
#>   N = 325/arm, Calculated Power = 0.464
#>   N = 350/arm, Calculated Power = 0.492
#>   N = 375/arm, Calculated Power = 0.519
#>   N = 400/arm, Calculated Power = 0.546
#>   N = 425/arm, Calculated Power = 0.571
#>   N = 450/arm, Calculated Power = 0.595
#>   N = 475/arm, Calculated Power = 0.618
#>   N = 500/arm, Calculated Power = 0.64
#>   N = 525/arm, Calculated Power = 0.662
#>   N = 550/arm, Calculated Power = 0.682
#>   N = 575/arm, Calculated Power = 0.701
#>   N = 600/arm, Calculated Power = 0.719
#>   N = 625/arm, Calculated Power = 0.737
#>   N = 650/arm, Calculated Power = 0.753
#>   N = 675/arm, Calculated Power = 0.769
#>   N = 700/arm, Calculated Power = 0.784
#>   N = 725/arm, Calculated Power = 0.798
#>   N = 750/arm, Calculated Power = 0.811
#> 
#> --- Calculation Summary ---
#> 
#> 
#> Table: Required Sample Size
#> 
#> | Target_Power| Required_N_per_Arm|
#> |------------:|------------------:|
#> |          0.8|                750|
print(ss_results$results_data)
#>   Target_Power Required_N_per_Arm
#> 1          0.8                750
print(ss_results$results_plot)
```

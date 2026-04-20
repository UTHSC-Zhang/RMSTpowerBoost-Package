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
#>   N = 50/arm, Calculated Power = 0.067
#>   N = 75/arm, Calculated Power = 0.081
#>   N = 100/arm, Calculated Power = 0.095
#>   N = 125/arm, Calculated Power = 0.108
#>   N = 150/arm, Calculated Power = 0.122
#>   N = 175/arm, Calculated Power = 0.135
#>   N = 200/arm, Calculated Power = 0.148
#>   N = 225/arm, Calculated Power = 0.162
#>   N = 250/arm, Calculated Power = 0.175
#>   N = 275/arm, Calculated Power = 0.188
#>   N = 300/arm, Calculated Power = 0.201
#>   N = 325/arm, Calculated Power = 0.214
#>   N = 350/arm, Calculated Power = 0.227
#>   N = 375/arm, Calculated Power = 0.24
#>   N = 400/arm, Calculated Power = 0.253
#>   N = 425/arm, Calculated Power = 0.266
#>   N = 450/arm, Calculated Power = 0.279
#>   N = 475/arm, Calculated Power = 0.292
#>   N = 500/arm, Calculated Power = 0.305
#>   N = 525/arm, Calculated Power = 0.317
#>   N = 550/arm, Calculated Power = 0.33
#>   N = 575/arm, Calculated Power = 0.343
#>   N = 600/arm, Calculated Power = 0.355
#>   N = 625/arm, Calculated Power = 0.367
#>   N = 650/arm, Calculated Power = 0.379
#>   N = 675/arm, Calculated Power = 0.391
#>   N = 700/arm, Calculated Power = 0.403
#>   N = 725/arm, Calculated Power = 0.415
#>   N = 750/arm, Calculated Power = 0.427
#>   N = 775/arm, Calculated Power = 0.438
#>   N = 800/arm, Calculated Power = 0.45
#>   N = 825/arm, Calculated Power = 0.461
#>   N = 850/arm, Calculated Power = 0.472
#>   N = 875/arm, Calculated Power = 0.483
#>   N = 900/arm, Calculated Power = 0.494
#>   N = 925/arm, Calculated Power = 0.505
#>   N = 950/arm, Calculated Power = 0.515
#>   N = 975/arm, Calculated Power = 0.526
#>   N = 1000/arm, Calculated Power = 0.536
#>   N = 1025/arm, Calculated Power = 0.546
#>   N = 1050/arm, Calculated Power = 0.556
#>   N = 1075/arm, Calculated Power = 0.566
#>   N = 1100/arm, Calculated Power = 0.575
#>   N = 1125/arm, Calculated Power = 0.585
#>   N = 1150/arm, Calculated Power = 0.594
#>   N = 1175/arm, Calculated Power = 0.603
#>   N = 1200/arm, Calculated Power = 0.612
#>   N = 1225/arm, Calculated Power = 0.621
#>   N = 1250/arm, Calculated Power = 0.63
#>   N = 1275/arm, Calculated Power = 0.639
#>   N = 1300/arm, Calculated Power = 0.647
#>   N = 1325/arm, Calculated Power = 0.655
#>   N = 1350/arm, Calculated Power = 0.663
#>   N = 1375/arm, Calculated Power = 0.671
#>   N = 1400/arm, Calculated Power = 0.679
#>   N = 1425/arm, Calculated Power = 0.687
#>   N = 1450/arm, Calculated Power = 0.694
#>   N = 1475/arm, Calculated Power = 0.702
#>   N = 1500/arm, Calculated Power = 0.709
#>   N = 1525/arm, Calculated Power = 0.716
#>   N = 1550/arm, Calculated Power = 0.723
#>   N = 1575/arm, Calculated Power = 0.73
#>   N = 1600/arm, Calculated Power = 0.737
#>   N = 1625/arm, Calculated Power = 0.743
#>   N = 1650/arm, Calculated Power = 0.75
#>   N = 1675/arm, Calculated Power = 0.756
#>   N = 1700/arm, Calculated Power = 0.762
#>   N = 1725/arm, Calculated Power = 0.768
#>   N = 1750/arm, Calculated Power = 0.774
#>   N = 1775/arm, Calculated Power = 0.78
#>   N = 1800/arm, Calculated Power = 0.785
#>   N = 1825/arm, Calculated Power = 0.791
#>   N = 1850/arm, Calculated Power = 0.796
#>   N = 1875/arm, Calculated Power = 0.801
#> 
#> --- Calculation Summary ---
#> 
#> 
#> Table: Required Sample Size
#> 
#> | Target_Power| Required_N_per_Arm|
#> |------------:|------------------:|
#> |          0.8|               1875|
print(ss_results$results_data)
#>   Target_Power Required_N_per_Arm
#> 1          0.8               1875
print(ss_results$results_plot)
```

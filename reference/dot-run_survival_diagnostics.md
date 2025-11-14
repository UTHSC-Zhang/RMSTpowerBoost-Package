# Run Survival Diagnostics for Pilot Data

This function performs survival diagnostics on pilot data, including a
log-rank test and Kaplan-Meier plot.

## Usage

``` r
.run_survival_diagnostics(
  pilot_data,
  time_var,
  status_var,
  arm_var,
  strata_var = NULL,
  alpha = 0.05
)
```

## Arguments

- pilot_data:

  A data frame containing the pilot data with survival information.

- time_var:

  A string specifying the name of the time variable in the pilot data.

- status_var:

  A string specifying the name of the status variable in the pilot data
  (1 for event, 0 for censored).

- arm_var:

  A string specifying the name of the treatment arm variable in the
  pilot data.

- strata_var:

  An optional string specifying the name of the stratification variable
  in the pilot data.

- alpha:

  The significance level for calculating the confidence interval
  (default is 0.05).

## Details

Run Survival Diagnostics for Pilot Data

# Simulated dataset: AFT log-normal, L = 12, n = 150

Columns:

- time: observed time

- status: event indicator (1 = event, 0 = censored)

- arm: treatment arm (0 = control, 1 = treatment)

- age: baseline age

- gender: factor with two levels

## Format

A data frame with 150 rows and 5 variables.

## Examples

``` r
if (FALSE) { # \dontrun{
data(aft_lognormal_L12_n150)
table(aft_lognormal_L12_n150$arm)
} # }
```

# Example dataset: small study (n = 50)

Columns:

- time: observed time

- status: event indicator (1 = event, 0 = censored)

- arm: treatment arm (0 = control, 1 = treatment)

- age: baseline age

- gender: factor with two levels

## Format

A data frame with 50 rows and 5 variables.

## Examples

``` r
if (FALSE) { # \dontrun{
data(rmst_small_50)
tapply(rmst_small_50$time, rmst_small_50$arm, median)
} # }
```

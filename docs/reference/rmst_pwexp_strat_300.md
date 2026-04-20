# Example dataset: stratified PH piecewise-exponential (n = 300)

Columns:

- time: observed time

- status: event indicator (1 = event, 0 = censored)

- arm: treatment arm (0 = control, 1 = treatment)

- age: baseline age

- gender: factor with two levels

## Format

A data frame with 300 rows and 5 variables.

## Examples

``` r
data(rmst_pwexp_strat_300)
#> Warning: data set 'rmst_pwexp_strat_300' not found
table(rmst_pwexp_strat_300$gender, rmst_pwexp_strat_300$arm)
#> Error: object 'rmst_pwexp_strat_300' not found
```

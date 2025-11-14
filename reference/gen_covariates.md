# Generate covariate matrix/data frame from a recipe

Generate covariate matrix/data frame from a recipe

## Usage

``` r
gen_covariates(n, covariates)
```

## Arguments

- n:

  sample size

- covariates:

  list(defs = list(...))

## Value

data.frame of covariates

## Examples

``` r
if (FALSE) { # \dontrun{
defs <- list(
  list(name="x", type="continuous", dist="normal", params=list(mean=0, sd=1)),
  list(name="z", type="categorical", dist="categorical",
       params=list(prob=c(0.3,0.7), labels=c("A","B")))
)
X <- gen_covariates(10, list(defs = defs))
} # }
```

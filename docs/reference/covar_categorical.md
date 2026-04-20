# Categorical covariate definition

Builds a multi-level categorical covariate definition list.

## Usage

``` r
covar_categorical(name, probs, labels = NULL)
```

## Arguments

- name:

  Column name in the generated data.

- probs:

  Numeric vector of category probabilities (must sum to 1).

- labels:

  Optional character vector of level labels. Length must equal
  `length(probs)`. Defaults to `"1"`, `"2"`, ...

## Value

A named list suitable as one element of the `covariates` argument.

## See also

[`covar_continuous`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/covar_continuous.md),
[`covar_binary`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/covar_binary.md),
[`rmst.sim`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.sim.md)

## Examples

``` r
covar_categorical("stage", probs = c(0.4, 0.35, 0.25),
                  labels = c("I", "II", "III"))
#> $name
#> [1] "stage"
#> 
#> $type
#> [1] "categorical"
#> 
#> $dist
#> [1] "categorical"
#> 
#> $params
#> $params$prob
#> [1] 0.40 0.35 0.25
#> 
#> $params$labels
#> [1] "I"   "II"  "III"
#> 
#> 
```

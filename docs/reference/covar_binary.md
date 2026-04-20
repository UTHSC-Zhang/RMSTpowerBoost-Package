# Binary covariate definition

Builds a binary (Bernoulli) covariate definition list.

## Usage

``` r
covar_binary(name, p = 0.5)
```

## Arguments

- name:

  Column name in the generated data.

- p:

  Probability of value 1. Default `0.5`.

## Value

A named list suitable as one element of the `covariates` argument.

## See also

[`covar_continuous`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/covar_continuous.md),
[`covar_categorical`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/covar_categorical.md),
[`rmst.sim`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.sim.md)

## Examples

``` r
covar_binary("female", p = 0.52)
#> $name
#> [1] "female"
#> 
#> $type
#> [1] "categorical"
#> 
#> $dist
#> [1] "bernoulli"
#> 
#> $params
#> $params$p
#> [1] 0.52
#> 
#> 
```

# Continuous covariate definition

Builds a covariate definition list for use in
[`rmst.sim`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.sim.md),
[`recipe_quick_aft`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/recipe_quick_aft.md),
or
[`recipe_quick_ph`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/recipe_quick_ph.md).

## Usage

``` r
covar_continuous(name, dist = "normal", ...)
```

## Arguments

- name:

  Column name in the generated data.

- dist:

  Distribution name. One of `"normal"` (default), `"lognormal"`,
  `"gamma"`, `"weibull"`, `"uniform"`, `"beta"`, `"t"`.

- ...:

  Named distribution parameters passed as the `params` list (e.g.,
  `mean = 0, sd = 1` for `dist = "normal"`).

## Value

A named list suitable as one element of the `covariates` argument.

## See also

[`covar_binary`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/covar_binary.md),
[`covar_categorical`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/covar_categorical.md),
[`rmst.sim`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.sim.md)

## Examples

``` r
covar_continuous("age", dist = "normal", mean = 50, sd = 10)
#> $name
#> [1] "age"
#> 
#> $type
#> [1] "continuous"
#> 
#> $dist
#> [1] "normal"
#> 
#> $params
#> $params$mean
#> [1] 50
#> 
#> $params$sd
#> [1] 10
#> 
#> 
covar_continuous("bmi", dist = "lognormal", meanlog = 3.2, sdlog = 0.2)
#> $name
#> [1] "bmi"
#> 
#> $type
#> [1] "continuous"
#> 
#> $dist
#> [1] "lognormal"
#> 
#> $params
#> $params$meanlog
#> [1] 3.2
#> 
#> $params$sdlog
#> [1] 0.2
#> 
#> 
```

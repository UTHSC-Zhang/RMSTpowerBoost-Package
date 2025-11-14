# Expand a recipe over a grid of values (list-only)

Expand a recipe over a grid of values (list-only)

## Usage

``` r
recipe_grid(base, vary)
```

## Arguments

- base:

  A recipe list.

- vary:

  Named list of vectors; keys are dotted paths.

## Value

A list of recipe lists.

## Examples

``` r
r <- recipe_quick_aft(60, "aft_lognormal", baseline=list(mu=2.7, sigma=0.6),
       treat_effect=-0.2, covariates=list(list(name="x", type="continuous", dist="normal",
       params=list(mean=0, sd=1))))
recipe_grid(r, list(n=c(60,80), "event_time.effects.treatment"=c(-0.2,-0.4)))
#> [[1]]
#> [[1]]$n
#> [1] 60
#> 
#> [[1]]$covariates
#> [[1]]$covariates$defs
#> [[1]]$covariates$defs[[1]]
#> [[1]]$covariates$defs[[1]]$name
#> [1] "x"
#> 
#> [[1]]$covariates$defs[[1]]$type
#> [1] "continuous"
#> 
#> [[1]]$covariates$defs[[1]]$dist
#> [1] "normal"
#> 
#> [[1]]$covariates$defs[[1]]$params
#> [[1]]$covariates$defs[[1]]$params$mean
#> [1] 0
#> 
#> [[1]]$covariates$defs[[1]]$params$sd
#> [1] 1
#> 
#> 
#> 
#> 
#> 
#> [[1]]$treatment
#> [[1]]$treatment$assignment
#> [1] "randomization"
#> 
#> [[1]]$treatment$allocation
#> [1] "1:1"
#> 
#> 
#> [[1]]$event_time
#> [[1]]$event_time$model
#> [1] "aft_lognormal"
#> 
#> [[1]]$event_time$baseline
#> [[1]]$event_time$baseline$mu
#> [1] 2.7
#> 
#> [[1]]$event_time$baseline$sigma
#> [1] 0.6
#> 
#> 
#> [[1]]$event_time$effects
#> [[1]]$event_time$effects$intercept
#> [1] 0
#> 
#> [[1]]$event_time$effects$treatment
#> [1] -0.2
#> 
#> [[1]]$event_time$effects$covariates
#> NULL
#> 
#> 
#> 
#> [[1]]$censoring
#> [[1]]$censoring$mode
#> [1] "target_overall"
#> 
#> [[1]]$censoring$target
#> [1] 0.25
#> 
#> [[1]]$censoring$admin_time
#> [1] Inf
#> 
#> 
#> [[1]]$seed
#> NULL
#> 
#> 
#> [[2]]
#> [[2]]$n
#> [1] 80
#> 
#> [[2]]$covariates
#> [[2]]$covariates$defs
#> [[2]]$covariates$defs[[1]]
#> [[2]]$covariates$defs[[1]]$name
#> [1] "x"
#> 
#> [[2]]$covariates$defs[[1]]$type
#> [1] "continuous"
#> 
#> [[2]]$covariates$defs[[1]]$dist
#> [1] "normal"
#> 
#> [[2]]$covariates$defs[[1]]$params
#> [[2]]$covariates$defs[[1]]$params$mean
#> [1] 0
#> 
#> [[2]]$covariates$defs[[1]]$params$sd
#> [1] 1
#> 
#> 
#> 
#> 
#> 
#> [[2]]$treatment
#> [[2]]$treatment$assignment
#> [1] "randomization"
#> 
#> [[2]]$treatment$allocation
#> [1] "1:1"
#> 
#> 
#> [[2]]$event_time
#> [[2]]$event_time$model
#> [1] "aft_lognormal"
#> 
#> [[2]]$event_time$baseline
#> [[2]]$event_time$baseline$mu
#> [1] 2.7
#> 
#> [[2]]$event_time$baseline$sigma
#> [1] 0.6
#> 
#> 
#> [[2]]$event_time$effects
#> [[2]]$event_time$effects$intercept
#> [1] 0
#> 
#> [[2]]$event_time$effects$treatment
#> [1] -0.2
#> 
#> [[2]]$event_time$effects$covariates
#> NULL
#> 
#> 
#> 
#> [[2]]$censoring
#> [[2]]$censoring$mode
#> [1] "target_overall"
#> 
#> [[2]]$censoring$target
#> [1] 0.25
#> 
#> [[2]]$censoring$admin_time
#> [1] Inf
#> 
#> 
#> [[2]]$seed
#> NULL
#> 
#> 
#> [[3]]
#> [[3]]$n
#> [1] 60
#> 
#> [[3]]$covariates
#> [[3]]$covariates$defs
#> [[3]]$covariates$defs[[1]]
#> [[3]]$covariates$defs[[1]]$name
#> [1] "x"
#> 
#> [[3]]$covariates$defs[[1]]$type
#> [1] "continuous"
#> 
#> [[3]]$covariates$defs[[1]]$dist
#> [1] "normal"
#> 
#> [[3]]$covariates$defs[[1]]$params
#> [[3]]$covariates$defs[[1]]$params$mean
#> [1] 0
#> 
#> [[3]]$covariates$defs[[1]]$params$sd
#> [1] 1
#> 
#> 
#> 
#> 
#> 
#> [[3]]$treatment
#> [[3]]$treatment$assignment
#> [1] "randomization"
#> 
#> [[3]]$treatment$allocation
#> [1] "1:1"
#> 
#> 
#> [[3]]$event_time
#> [[3]]$event_time$model
#> [1] "aft_lognormal"
#> 
#> [[3]]$event_time$baseline
#> [[3]]$event_time$baseline$mu
#> [1] 2.7
#> 
#> [[3]]$event_time$baseline$sigma
#> [1] 0.6
#> 
#> 
#> [[3]]$event_time$effects
#> [[3]]$event_time$effects$intercept
#> [1] 0
#> 
#> [[3]]$event_time$effects$treatment
#> [1] -0.4
#> 
#> [[3]]$event_time$effects$covariates
#> NULL
#> 
#> 
#> 
#> [[3]]$censoring
#> [[3]]$censoring$mode
#> [1] "target_overall"
#> 
#> [[3]]$censoring$target
#> [1] 0.25
#> 
#> [[3]]$censoring$admin_time
#> [1] Inf
#> 
#> 
#> [[3]]$seed
#> NULL
#> 
#> 
#> [[4]]
#> [[4]]$n
#> [1] 80
#> 
#> [[4]]$covariates
#> [[4]]$covariates$defs
#> [[4]]$covariates$defs[[1]]
#> [[4]]$covariates$defs[[1]]$name
#> [1] "x"
#> 
#> [[4]]$covariates$defs[[1]]$type
#> [1] "continuous"
#> 
#> [[4]]$covariates$defs[[1]]$dist
#> [1] "normal"
#> 
#> [[4]]$covariates$defs[[1]]$params
#> [[4]]$covariates$defs[[1]]$params$mean
#> [1] 0
#> 
#> [[4]]$covariates$defs[[1]]$params$sd
#> [1] 1
#> 
#> 
#> 
#> 
#> 
#> [[4]]$treatment
#> [[4]]$treatment$assignment
#> [1] "randomization"
#> 
#> [[4]]$treatment$allocation
#> [1] "1:1"
#> 
#> 
#> [[4]]$event_time
#> [[4]]$event_time$model
#> [1] "aft_lognormal"
#> 
#> [[4]]$event_time$baseline
#> [[4]]$event_time$baseline$mu
#> [1] 2.7
#> 
#> [[4]]$event_time$baseline$sigma
#> [1] 0.6
#> 
#> 
#> [[4]]$event_time$effects
#> [[4]]$event_time$effects$intercept
#> [1] 0
#> 
#> [[4]]$event_time$effects$treatment
#> [1] -0.4
#> 
#> [[4]]$event_time$effects$covariates
#> NULL
#> 
#> 
#> 
#> [[4]]$censoring
#> [[4]]$censoring$mode
#> [1] "target_overall"
#> 
#> [[4]]$censoring$target
#> [1] 0.25
#> 
#> [[4]]$censoring$admin_time
#> [1] Inf
#> 
#> 
#> [[4]]$seed
#> NULL
#> 
#> 
```

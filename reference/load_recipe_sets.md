# Load datasets from a recipe-sets manifest

Reads a `manifest.rds` created by
[`generate_recipe_sets()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/generate_recipe_sets.md),
loads one dataset per row (preferring `rds` → `rdata` → `csv` → `txt`),
restores attribute `"achieved_censoring"`, and returns a named list of
`list(data = <data.frame>, meta = <list>)`.

## Usage

``` r
load_recipe_sets(manifest_path)
```

## Arguments

- manifest_path:

  Path to `manifest.rds`.

## Value

A named list where each element is `list(data=..., meta=...)`.

## Examples

``` r
if (FALSE) { # \dontrun{
sets <- load_recipe_sets("checks/manifest.rds")
names(sets)
str(sets[[1]]$meta)
} # }
```

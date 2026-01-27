# Generate simulated datasets across scenario combinations (TXT/CSV/RDS/RData)

Builds a grid of scenarios from a base recipe and a named list of
variations, simulates one or more replicates per scenario, and writes
the datasets to files in a target folder as `.txt` (tab), `.csv`,
`.rds`, and/or `.RData`. A manifest is written as `manifest.rds`. No
YAML and **no L/tau**.

## Usage

``` r
generate_recipe_sets(
  base_recipe,
  vary = list(),
  out_dir,
  formats = c("txt", "csv", "rds", "rdata"),
  n_reps = 1L,
  seed_base = NULL,
  filename_template = "sc{scenario_id}_r{rep}"
)
```

## Arguments

- base_recipe:

  A recipe list (use
  [`validate_recipe()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/validate_recipe.md)
  if needed).

- vary:

  Named list; keys are dotted paths inside the recipe (e.g., `"n"`,
  `"censoring.target"`, `"event_time.effects.treatment"`).

- out_dir:

  Directory to write datasets and `manifest.rds` (created if missing).

- formats:

  Character vector subset of `c("txt","csv","rds","rdata")`.

- n_reps:

  Integer; number of replicates per scenario.

- seed_base:

  Optional integer; per-rep seed computed as
  `seed_base + scenario_id*1000 + rep`.

- filename_template:

  Base filename (no extension) with placeholders: `"{scenario_id}"`,
  `"{rep}"` and any dotted path used in `vary`.

## Value

Invisibly returns the manifest `data.frame` and writes `manifest.rds`.

## Examples

``` r
covs <- list(list(name="x", type="continuous", dist="normal", params=list(mean=0, sd=1)))
rec <- recipe_quick_aft(60, "aft_lognormal",
         baseline=list(mu=2.7, sigma=0.6), treat_effect=-0.2,
         covariates=covs, target_censoring=0.2)
dir.create("checks", showWarnings = FALSE)
man <- generate_recipe_sets(rec, vary=list(n=c(60,80)), out_dir="checks",
         formats=c("csv","rds"), n_reps=1, seed_base=123)
```

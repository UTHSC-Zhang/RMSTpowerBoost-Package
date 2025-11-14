# Rebuild manifest for an existing output directory (no re-simulation)

Reads datasets already written in `out_dir` and reconstructs a
`manifest.rds` with rich metadata (model, baseline, effects, etc.).

## Usage

``` r
rebuild_manifest(
  base_recipe,
  vary,
  out_dir,
  filename_template = "sc{scenario_id}_r{rep}"
)
```

## Arguments

- base_recipe:

  A validated recipe list (use
  [`validate_recipe()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/validate_recipe.md)
  if needed).

- vary:

  Named list used originally in
  [`generate_recipe_sets()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/generate_recipe_sets.md)
  for the grid.

- out_dir:

  Directory that already contains the datasets.

- filename_template:

  The same template you used when writing files (default
  `"sc{scenario_id}_r{rep}"`). It may also include tokens for dotted
  paths from `vary`.

## Value

The rebuilt manifest (also writes `manifest.rds` in `out_dir`).

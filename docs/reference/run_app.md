# Launch the RMSTpowerBoost Shiny Application

Launches the Shiny application bundled with the package. App-specific
dependencies are installed only when they are needed.

## Usage

``` r
run_app(install_missing = TRUE, repos = getOption("repos"))
```

## Arguments

- install_missing:

  Logical; if `TRUE`, prompt to install missing app dependencies.

- repos:

  CRAN mirror(s) passed to
  [`utils::install.packages()`](https://rdrr.io/r/utils/install.packages.html)
  when installing missing app dependencies.

## Value

Invisible return value from
[`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html).

## Examples

``` r
if (FALSE) { # interactive()
RMSTpowerBoost::run_app()
}
```

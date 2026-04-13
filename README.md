# **RMSTpowerBoost: Power and Sample Size Calculations for RMST-Based Trials** [![codecov](https://codecov.io/gh/UTHSC-Zhang/RMSTpowerBoost-Package/graph/badge.svg?token=MK7AIUQHBC)](https://codecov.io/gh/UTHSC-Zhang/RMSTpowerBoost-Package)

## Overview

`RMSTpowerBoost` provides power and sample size tools for study designs that use restricted mean survival time (RMST). The package supports direct RMST modeling with analytical and bootstrap-based procedures for settings that include linear covariate adjustment, stratification, semiparametric additive effects, and covariate-dependent censoring.

The package includes both an R interface and a Shiny application for interactive use.

## Key Features

- Linear IPCW-based RMST power and sample size calculations.
- Additive and multiplicative stratified models for multi-center or highly stratified studies.
- Bootstrap-based semiparametric GAM procedures for nonlinear covariate effects.
- Analytical and simulation-based methods for covariate-dependent censoring under a single censoring mechanism.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("UTHSC-Zhang/RMSTpowerBoost-Package")
```

## Codecov Sunburst

![Codecov Sunburst](https://codecov.io/gh/arnabaich96/RMSTpowerBoost-Package/graphs/sunburst.svg?token=5L9OGLTMU6)

## Shiny App

Interactive web application:

[Launch the Shiny App](https://arnab96.shinyapps.io/uthsc-app/)

## Local App Sync + Coverage Gate

Run this from the package repo root to sync app artifacts from `../RMSTpowerBoost-App`, validate the synced app, run tests, generate coverage artifacts, and enforce the `>90%` overall/per-file gate:

```r
Rscript tools/sync_app_from_repo.R ../RMSTpowerBoost-App && \
Rscript tools/validate_synced_app.R && \
Rscript -e "testthat::test_local('.', reporter='summary')" && \
Rscript tools/generate_coverage_artifacts.R && \
Rscript tools/check_coverage_thresholds.R
```

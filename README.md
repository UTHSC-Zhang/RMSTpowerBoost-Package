# **RMSTpowerBoost: Power and Sample Size Calculations for RMST-Based Trials** 

## Overview

`RMSTpowerBoost` provides power and sample size tools for study designs that use restricted mean survival time (RMST) as a summary metric of time-to-event outcomes. The package supports covariate adjustment with analytical and simulation-based procedures for settings that include nonproportional hazards, stratification or multi-center effects, and dependent censoring.

The package includes both an R interface and a Shiny application for interactive use.

## Key Features

- Linear IPCW-based RMST power and sample size calculations.
- Additive and multiplicative stratified models for multi-center or highly stratified studies.
- Bootstrap-based semiparametric GAM procedures for nonlinear covariate effects.
- Analytical and simulation-based methods for covariate-dependent censoring under a single censoring mechanism.

## Unified Interface

The main package surface is organized around three top-level helpers:

- `rmst.sim()` generates pilot or reference survival data across the supported AFT and PH model families.
- `rmst.power()` computes a power curve from a `Surv(time, status) ~ ...` formula.
- `rmst.ss()` searches for the minimum sample size that reaches a target power.

```r
library(RMSTpowerBoost)
library(survival)

set.seed(42)
pilot <- rmst.sim(
  n = 150,
  model = "aft_lognormal",
  baseline = list(mu = 2.2, sigma = 0.5),
  treat_effect = -0.3,
  covariates = list(covar_continuous("age"), covar_binary("female")),
  L = 12,
  seed = 42
)

power_fit <- rmst.power(
  Surv(time, status) ~ age + female,
  data = pilot,
  arm = "arm",
  sample_sizes = c(75, 100, 125),
  L = 12
)
```

`rmst.power()` and `rmst.ss()` route automatically based on a few key arguments:

- `type = "analytical"` or `"boot"` selects the analytical or bootstrap path.
- `strata` and `strata_type` activate additive or multiplicative stratified models.
- `dep_cens = TRUE` selects the dependent-censoring workflow.
- Smooth terms such as `s(age)` in the formula select the GAM bootstrap workflow.

`seed` in `rmst.sim()` is stored in the recipe metadata for provenance, but reproducible simulation still comes from calling `set.seed()` before `rmst.sim()`.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("UTHSC-Zhang/RMSTpowerBoost-Package")
```

## Shiny App

Interactive web application:

[Launch the Shiny App](https://arnab96.shinyapps.io/uthsc-app/)

# **RMSTpowerBoost: Power and Sample Size Calculations for RMST-Based Trials** [![codecov](https://codecov.io/gh/UTHSC-Zhang/RMSTpowerBoost-Package/graph/badge.svg?token=MK7AIUQHBC)](https://app.codecov.io/gh/UTHSC-Zhang/RMSTpowerBoost-Package)

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

## Shiny App

Interactive web application:

[Launch the Shiny App](https://arnab96.shinyapps.io/uthsc-app/)



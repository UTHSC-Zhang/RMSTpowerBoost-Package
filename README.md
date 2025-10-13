
# RMSTpowerBoost: Power and Sample Size Calculations for RMST-Based Trials

---

## Overview

The analysis of time-to-event data is moving beyond the proportional hazards assumption, and the Restricted Mean Survival Time (RMST) has emerged as a clinically and causally interpretable alternative to the hazard ratio. However, tools for designing studies based on modern, direct modeling approaches for the RMST have been lacking.

`RMSTpowerBoost` fills this critical gap by implementing a variety of advanced statistical methods for study design, allowing researchers to accurately plan trials under complex scenarios. This software suite consists of two primary components: an interactive web application for ease of use, and a comprehensive R package for flexibility and reproducibility.

### Key Features

* **Multiple Model Types**: Handles various data structures and assumptions:
    * **Non-Stratified Models (Linear & GAM)**: For single-center studies with linear or non-linear covariate effects.
    * **Dependent Censoring Models**: Appropriate for settings with competing risks.
    * **Stratified Models (Additive & Multiplicative)**: For multi-center trials with different treatment effect assumptions.
* **Dual Calculation Methods**: For many models, the package offers both:
    * A fast **analytical** approach based on asymptotic variance formulas.
    * A robust **bootstrap** (simulation-based) approach for enhanced accuracy.
* **Interactive Application**: Includes a user-friendly Shiny application for point-and-click analysis, making these advanced methods accessible to a broad audience.

## Installation

You can install the development version of `RMSTpowerBoost` from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("UTHSC-Zhang/RMSTpowerBoost-package")
```

## Interactive Shiny Application

For users who prefer a graphical user interface, an interactive Shiny application is available. The app provides a point-and-click interface to all the package's functionalities.

You can access the application directly in your web browser by following this link: [**Launch Web App**](https://arnab96.shinyapps.io/uthsc-app/)




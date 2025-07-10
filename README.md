

# RMSTSS: Power and Sample Size for RMST-Based Trials

[](https://www.google.com/search?q=https://CRAN.R-project.org/package%3DRMSTSS)
[](https://www.google.com/search?q=https://github.com/arnabaich/RMSTSS/actions/workflows/R-CMD-check.yaml)
`RMSTSS` provides a comprehensive suite of tools for calculating statistical power and sample size for clinical trials that use the Restricted Mean Survival Time (RMST) as the primary endpoint.

## Overview

The analysis of time-to-event data is moving beyond the proportional hazards assumption, and the RMST has emerged as a robust and clinically intuitive alternative to the hazard ratio. However, tools for designing studies based on modern, direct-modeling approaches for the RMST have been lacking.

`RMSTSS` fills this critical gap by implementing a variety of advanced statistical methods for study design, allowing researchers to accurately plan trials under complex scenarios. More detailed usage of this package is described in the vignette, which can be accessed using the following link [**Click Here**](https://uthsc-zhang.github.io/RMSTSS-Package/articles/RMSTSS.html).

### Key Features

  * **Multiple Model Types**: Handles various data structures and assumptions.
      * **Linear IPCW Models**: Standard direct regression modeling.
      * **Additive & Multiplicative Stratified Models**: Efficiently handles stratification by variables with many levels (e.g., clinical centers).
      * **Semiparametric GAM Models**: Allows for flexible, non-linear covariate effects.
      * **Dependent Censoring Models**: Appropriate for settings with competing risks.
  * **Dual Calculation Methods**: For many models, the package offers both:
      * A fast **analytical** approach based on asymptotic variance formulas.
      * A robust **bootstrap** (simulation-based) approach for enhanced accuracy and fewer distributional assumptions.
  * **Interactive Application**: Includes a user-friendly Shiny application for point-and-click analysis.

## Installation

You can install the development version of `RMSTSS` from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("https://github.com/UTHSC-Zhang/RMSTSS-package.git")
```

## Quick Example

Here is a basic example of calculating statistical power for a range of sample sizes using a linear RMST model. The calculation is based on pilot data from the `veteran` dataset, adjusting for the Karnofsky performance score (`karno`).

```r
library(RMSTSS)
library(survival)

# Prepare the veteran dataset for analysis
vet <- veteran
vet$arm <- ifelse(vet$trt == 1, 0, 1)

# Calculate power for several sample sizes per arm
power_results <- linear.power.analytical(
  pilot_data = vet,
  time_var = "time",
  status_var = "status",
  arm_var = "arm",
  linear_terms = "karno",
  sample_sizes = c(100, 150, 200, 250),
  tau = 270 # Truncation time of 9 months
)

# View the results data frame
print(power_results$results_data)

# View the power curve plot
print(power_results$results_plot)
```

## Interactive Shiny Application

For users who prefer a graphical user interface, an interactive Shiny application is available. The app provides a point-and-click interface to all the package's functionalities, allowing you to upload data, map variables, select models, and run analyses without writing any code.

You can access the application directly in your web browser by following this link: [**Click Here**](https://arnab96.shinyapps.io/uthsc-app/)

## Citation

If you use `RMSTSS` in your research, please cite it as follows:

```
Aich, A. (2025). RMSTSS: Sample Size and Power Calculations for RMST-based Clinical Trials. R package version 0.1.0. https://github.com/arnabaich/RMSTSS
```

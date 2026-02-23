

# **RMSTpowerBoost: Power and Sample Size Calculations for RMST-Based Trials** [![codecov](https://codecov.io/gh/arnabaich96/RMSTpowerBoost-Package/graph/badge.svg?token=5L9OGLTMU6)](https://codecov.io/gh/arnabaich96/RMSTpowerBoost-Package)

---

## 🔍 Overview

As clinical trials and observational studies increasingly move beyond the restrictive **proportional hazards (PH) assumption**, researchers need flexible alternatives for analyzing time-to-event data. One such approach—**Restricted Mean Survival Time (RMST)**—has gained traction due to its clear clinical interpretation and robustness in non-PH settings.

Yet, while RMST-based analysis is growing in popularity, tools for **study design and power calculation** within this framework have lagged behind.

**`RMSTpowerBoost` bridges that gap**—offering a powerful suite of methods to design, simulate, and plan RMST-based trials, even under complex censoring mechanisms or nonlinear treatment effects. Whether you prefer point-and-click simplicity or script-based control, `RMSTpowerBoost` has you covered.

---

## ✨ Key Features

### ✅ Broad Model Support

Design studies under a range of realistic conditions:

* **Non-Stratified Models**

  * *Linear*: Simple, interpretable modeling.
  * *GAM (Generalized Additive Models)*: Flexible for nonlinear covariate effects.

* **Dependent Censoring Models**

  * Account for censoring that's related to baseline covariates.

* **Stratified Models**

  * *Additive and Multiplicative*: Ideal for multi-center trials with center-specific treatment effects.

### 🧠 Dual Estimation Methods

For many models, you can choose between:

* **Analytical Approach**: Based on asymptotic variance formulas. Fast and efficient for large-scale simulations.
* **Bootstrap-Based Method**: Simulation-driven, providing greater robustness, especially in smaller or complex datasets.

### 🖱️ Interactive Web App

No R experience? No problem. A full-featured **Shiny web application** is available, making advanced RMST-based design accessible through an intuitive graphical interface.

---

## 📦 Installation

Install the latest development version directly from GitHub:

```r
# If not already installed
install.packages("remotes")

# Install RMSTpowerBoost
remotes::install_github("UTHSC-Zhang/RMSTpowerBoost-package")
```

---

## Codecov Sunburst
![Codecov Sunburst](https://codecov.io/gh/arnabaich96/RMSTpowerBoost-Package/graphs/sunburst.svg?token=5L9OGLTMU6)
---

## 🌐 Try the Shiny App

Explore the tool without writing a single line of code:

👉 [**Launch the Interactive Web App**](https://arnab96.shinyapps.io/uthsc-app/)

Perfect for teaching, preliminary analyses, or rapid prototyping.

---

If you need to design a study using RMST, whether for a randomized trial, observational cohort, or simulation study, `RMSTpowerBoost` delivers both flexibility and precision—helping you move confidently beyond the hazard ratio.

---

## Local App Sync + Coverage Gate

Run this from the package repo root to sync app artifacts from `../RMSTpowerBoost-App`, validate the synced app, run tests, generate coverage artifacts, and enforce the `>90%` overall/per-file gate:

```r
Rscript tools/sync_app_from_repo.R ../RMSTpowerBoost-App && \
Rscript tools/validate_synced_app.R && \
Rscript -e "testthat::test_local('.', reporter='summary')" && \
Rscript tools/generate_coverage_artifacts.R && \
Rscript tools/check_coverage_thresholds.R
```

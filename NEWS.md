# RMSTpowerBoost 1.0.3

## Statistical corrections

* **IPCW event indicator**: all IPCW-based estimators (`linear.*`, `additive.*.analytical`, `MS.*.analytical`, `DC.*.analytical`) now use the truncated-outcome indicator \(\Delta^Y = 1\) if the event occurs before `L` *or follow-up reaches `L`*, instead of the raw event indicator. Subjects censored after `L` are complete cases for the RMST at `L` (Zhang & Schaubel, 2024, Biometrical Journal; Tian et al., 2014); dropping them attenuated treatment-effect estimates and distorted power.
* **Censoring model time scale**: censoring distributions (Kaplan-Meier and stratified Cox) are now fit on the original follow-up time rather than the `L`-truncated time, which had recorded post-`L` censorings as censoring events at exactly `L`.
* **Dependent-censoring model**: `DC.power.analytical()` / `DC.ss.analytical()` previously weighted *all* subjects by \(1/\hat G(Y)\), including censored subjects with their censored time as the outcome, yielding a severely attenuated estimator. Weights are now \(\Delta^Y/\hat G(Y)\) and the regression uses complete cases only.
* **Stratified baseline hazard lookup**: `additive.*.analytical()` and `MS.*.analytical()` matched `basehaz()` strata against `"var=level"` labels; current versions of the survival package return bare `"level"` labels, so no stratum ever matched and every IPCW weight silently collapsed to 1 (unweighted complete-case analysis). Both label formats are now accepted and an informative error is raised if a stratum cannot be located.
* **Multiplicative model variance scaling**: the asymptotic variance in `MS.*.analytical()` is now scaled per enrolled pilot subject rather than per observed event, which had overstated power under censoring.

These corrections were validated against the estimator of Zhang & Schaubel (2024): on simulated data with a known true RMST difference, the linear, additive, and dependent-censoring estimators are now consistent (see `tests/validation-testR-alignment.R`). Simulation calibration confirms the sandwich variances (empirical SD vs. analytic SE within 2% at n = 4000; sandwich z-test type-I error 0.058 at nominal 0.05).

* **Reported treatment-effect standard errors**: the `model_output$treatment_effect` tables of `linear.*.analytical()` and `DC.*.analytical()` displayed the asymptotic n = 1 standard error instead of the pilot-scale standard error (about 63x too large for a pilot of 4000), making the displayed confidence intervals meaningless. They are now scaled by `1/sqrt(n_pilot)`.
* **Stratified bootstrap hypothesis**: `GAM.power.boot()`/`GAM.ss.boot()` (stratified) and `MS.power.boot()`/`MS.ss.boot()` previously fit stratum-by-arm interactions and took the minimum p-value across strata without multiplicity adjustment, inflating the effective type-I error (about 1 - (1-alpha)^J for J strata) and overstating power. They now fit stratum-specific intercepts with a common treatment effect, matching the single-effect hypothesis of the analytical counterparts.
* **Documented**: bootstrap linear power uses model-based weighted-lm p-values (a conservative test relative to the sandwich test used by the analytic method); the multiplicative bootstrap drops non-positive pseudo-observations before taking logs.
* **Robustness**: degenerate pilots with no usable IPCW weights now fail informatively instead of with "object not found"; `.meta$n_col` now reports the correct column name for additive-stratified and GAM results; routing `strata_type = "additive"` with `type = "boot"` now messages that the GAM pseudo-observation bootstrap is used.

## Shiny app fixes

* The bundled app's Analytical branch passed an unsupported `point_cb` argument to the package's power/sample-size functions, causing an "unused argument" error on every analytical run. The argument was removed and the computed points are now streamed to the live plot after each call returns.
* The app's "Repeated" method previously simulated the power of a **log-rank test** (ignoring the truncation time `L`, the selected model, and covariates). It now resamples the pilot data and tests the treatment effect on the L-truncated RMST with an IPCW-weighted linear model (Delta-Y complete-case weights), matching the package's linear IPCW methodology; strata enter as fixed effects with a common treatment effect.

## Other changes

* Refined the unified interface around `rmst.power()`, `rmst.ss()`, and `rmst.sim()`.
* Clarified simulation seed handling so documentation and examples match the implementation: `seed` is retained as recipe metadata, while reproducible simulation still relies on `set.seed()`.

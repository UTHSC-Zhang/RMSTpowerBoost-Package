# Development version

* The four linear IPCW power and sample-size functions now use conditional Cox
  censoring survival with treatment and `linear_terms` as predictors and Breslow
  ties, replacing pooled Kaplan-Meier weights. The model is fit on original
  follow-up and evaluated at each subject's truncated outcome. Bootstrap samples
  refit the model. No-censoring data retain unit weights. Analytical sandwich and
  bootstrap WLS inference are unchanged; the analytical sandwich does not
  incorporate Cox fitting uncertainty.
* Weight capping has been removed from the four linear IPCW functions entirely,
  with no opt-in fallback. The weights are used exactly as fitted, so that the
  estimator and its simplified sandwich variance reproduce Equations (11)-(12) of
  Zhang & Schaubel (2024), and so that the analytical and bootstrap engines weight
  identically. The former 99th-percentile cap fired unconditionally rather than
  only on unstable weights, and compressed the upper tail that IPCW exists to
  up-weight. `model_output$censoring_weights` now reports only `raw_summary`;
  the `cap_value` and `capped_fraction` entries are gone, as are the internal
  `.linear_cox_weights(cap =)` and `.estimate_linear_params(cap_weights =)`
  arguments. `.linear_cox_weights()` now returns a plain numeric vector. Point
  estimates and standard errors shift slightly; in the reference pilot the
  treatment-effect SE moved by well under 1%.

* `linear.power.analytical()` and `linear.ss.analytical()` gain a `test` argument
  selecting which hypothesis test the reported power refers to. The default,
  `test = "wls"`, is the weighted-least-squares t-test that `linear.power.boot()`
  simulates and that external reference implementations use: the estimate's
  sampling distribution is taken from the sandwich variance and the rejection
  threshold from the model-based WLS standard error, which is not consistent here
  because IPCW weights are sampling weights rather than inverse-variance weights.
  On a 4000-subject pilot the WLS standard error ran about 20% above the measured
  bootstrap standard deviation of the estimate, and the analytic and simulated
  power now agree within Monte Carlo error (0.125 vs 0.110 at 500/arm, 0.247 vs
  0.263 at 1000/arm, 0.491 vs 0.515 at 2000/arm, B = 400). `test = "sandwich"`
  restores the previous behavior, the power of the sandwich Wald test of Zhang &
  Schaubel (2024), which is the more efficient test but not the one the bootstrap
  engine runs. `model_output$treatment_effect$std_error` and
  `variance_components$se_effect_n1` follow the selected test; both standard errors
  are always returned in `variance_components`.

* Weight capping has been removed from the two stratified analytical engines,
  `additive.power.analytical()` / `additive.ss.analytical()` and
  `MS.power.analytical()` / `MS.ss.analytical()`, for the same reasons it was
  removed from the linear engine: the 99th-percentile cap fired unconditionally
  rather than only on unstable weights, compressed exactly the upper tail that
  IPCW exists to up-weight, and broke the correspondence with Equations (11)-(12)
  of Zhang & Schaubel (2024), which derive the estimator and its simplified
  sandwich under the fitted weights. `.ipcw_stratified_cox_weights()` now returns
  `Delta_Y * W_hat` directly, so incomplete subjects weigh exactly zero in both
  engines, and it fails loudly on invalid censoring survival probabilities
  instead of letting the cap absorb an overflow. No-censoring data retain unit
  weights without fitting a Cox model. `model_output$censoring_weights` now
  reports only `raw_summary`; the `cap_value` and `capped_fraction` entries are
  gone, as is the internal `.cap_weights()` helper. Point estimates, standard
  errors, and power shift accordingly.

  For the multiplicative engine, `censoring_weights$raw_summary` now includes the
  zeros for subjects whose truncated outcome is unobserved, where it previously
  summarised `1/G` across everyone. The additive engine also now reports a
  non-convergent censoring model as an error, as the multiplicative engine
  already did, rather than swallowing the warning from a bare `coxph()` call.

* `MS.power.analytical()` and `MS.ss.analytical()` gain the same `test` argument
  as their linear counterparts, but it defaults the other way, to
  `test = "sandwich"`. That preserves the existing numbers and keeps the engine
  aligned with the Wang et al. (2019) inference. The linear default is `"wls"`
  because `linear.power.boot()` simulates exactly that test; no bootstrap engine
  in this package simulates an IPCW-weighted t-test for a stratified model, since
  `MS.power.boot()` fits an unweighted model to jackknife pseudo-observations and
  the additive model routes to the GAM pseudo-observation bootstrap. Under
  `test = "wls"` the rejection threshold comes from the model-based WLS standard
  error while the sampling distribution stays on the sandwich; the coefficient
  table, `treatment_effect$std_error` and `variance_components$se_effect_n1`
  follow the selected test, and both standard errors are always returned in
  `variance_components`. The additive engine gets no `test` argument: it solves
  its estimating equation in closed form by stratum-centering, so there is no
  weighted-least-squares fit and no second variance to choose between.

* `rmst.power()` and `rmst.ss()` gain a `test` argument, which was previously
  reachable only by calling the engine functions directly. It defaults to `NULL`,
  leaving each engine on its own default, and emits a message when supplied for
  an engine that does not accept it. `summary()` now prints the selected test and
  both standard errors when the engine reports them.

* The dependent-censoring engine is unchanged: `DC.power.analytical()` and
  `DC.ss.analytical()` still cap at the 99th percentile and still floor the
  censoring survival at `1e-6`, and they gain no `test` argument.

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

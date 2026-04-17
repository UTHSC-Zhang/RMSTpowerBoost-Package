# Tests for R/rmst-interface.R
# Covers: .parse_rmst_formula, .resolve_strata, .rmst_route, .rmst_model_label,
#         rmst.power, rmst.ss, and all S3 methods (print/summary/plot).

# ── Shared pilot data helpers ─────────────────────────────────────────────────

.iface_pilot <- function(n = 50L, with_strata = FALSE) {
  set.seed(8888L)
  dat <- data.frame(
    time   = stats::rexp(n, rate = 0.15) + 0.1,
    status = stats::rbinom(n, 1L, 0.75),
    arm    = rep(0L:1L, length.out = n),
    age    = stats::rnorm(n, 60, 8),
    stringsAsFactors = FALSE
  )
  if (with_strata)
    dat$region <- factor(rep(c("A", "B"), length.out = n))
  dat
}

.dc_pilot <- function(n = 200L, seed = 42L) {
  set.seed(seed)
  dat <- data.frame(
    arm = rep(0L:1L, each = n %/% 2L),
    age = stats::rnorm(n, 60, 8)
  )
  dat$time <- stats::rexp(n, rate = ifelse(dat$arm == 0L, 0.10, 0.08))
  ev <- sample(0:2, n, replace = TRUE, prob = c(0.6, 0.2, 0.2))
  dat$status     <- as.integer(ev == 0L)
  dat$comp_event <- as.integer(ev == 1L)
  dat
}

pilot   <- .iface_pilot(50L)
pilot_s <- .iface_pilot(60L, with_strata = TRUE)
dc_dat  <- .dc_pilot(200L)

# ── .parse_rmst_formula ───────────────────────────────────────────────────────

test_that(".parse_rmst_formula parses Surv(t, s) ~ x correctly", {
  p <- RMSTpowerBoost:::.parse_rmst_formula(
    survival::Surv(time, status) ~ age
  )
  expect_identical(p$time_var,     "time")
  expect_identical(p$status_var,   "status")
  expect_identical(p$linear_terms, "age")
  expect_null(p$smooth_terms)
})

test_that(".parse_rmst_formula detects s() smooth terms", {
  p <- RMSTpowerBoost:::.parse_rmst_formula(
    survival::Surv(time, status) ~ s(age)
  )
  expect_null(p$linear_terms)
  expect_identical(p$smooth_terms, "age")
})

test_that(".parse_rmst_formula handles mixed smooth + linear terms", {
  p <- RMSTpowerBoost:::.parse_rmst_formula(
    survival::Surv(time, status) ~ age + s(bmi)
  )
  expect_identical(p$linear_terms, "age")
  expect_identical(p$smooth_terms, "bmi")
})

test_that(".parse_rmst_formula handles intercept-only RHS (no covariates)", {
  p <- RMSTpowerBoost:::.parse_rmst_formula(
    survival::Surv(time, status) ~ 1
  )
  expect_null(p$linear_terms)
  expect_null(p$smooth_terms)
})

test_that(".parse_rmst_formula errors when LHS is not Surv()", {
  expect_error(
    RMSTpowerBoost:::.parse_rmst_formula(time ~ age),
    "Surv"
  )
})

# ── .resolve_strata ───────────────────────────────────────────────────────────

test_that(".resolve_strata returns NULL for NULL input", {
  expect_null(RMSTpowerBoost:::.resolve_strata(NULL))
})

test_that(".resolve_strata returns character string unchanged", {
  expect_identical(RMSTpowerBoost:::.resolve_strata("region"), "region")
})

test_that(".resolve_strata extracts column name from one-sided formula", {
  expect_identical(RMSTpowerBoost:::.resolve_strata(~region), "region")
})

test_that(".resolve_strata errors when formula has two terms", {
  expect_error(
    RMSTpowerBoost:::.resolve_strata(~a + b),
    "exactly one column"
  )
})

test_that(".resolve_strata errors on unsupported input type", {
  expect_error(RMSTpowerBoost:::.resolve_strata(42))
  expect_error(RMSTpowerBoost:::.resolve_strata(list("x")))
})

# ── .rmst_route ───────────────────────────────────────────────────────────────

test_that(".rmst_route dep_cens+boot auto-switches to analytical with message", {
  expect_message(
    r <- RMSTpowerBoost:::.rmst_route(TRUE, FALSE, NULL, "additive", "boot"),
    "switching type"
  )
  expect_identical(r$model,       "DC")
  expect_identical(r$method_used, "analytical")
  expect_identical(r$fn_power,    RMSTpowerBoost::DC.power.analytical)
  expect_identical(r$fn_ss,       RMSTpowerBoost::DC.ss.analytical)
})

test_that(".rmst_route dep_cens+analytical returns DC analytical silently", {
  expect_silent(
    r <- RMSTpowerBoost:::.rmst_route(TRUE, FALSE, NULL, "additive", "analytical")
  )
  expect_identical(r$model,       "DC")
  expect_identical(r$method_used, "analytical")
})

test_that(".rmst_route has_smooth+analytical auto-switches to boot with message", {
  expect_message(
    r <- RMSTpowerBoost:::.rmst_route(FALSE, TRUE, NULL, "additive", "analytical"),
    "switching type"
  )
  expect_identical(r$model,       "GAM")
  expect_identical(r$method_used, "boot")
})

test_that(".rmst_route has_smooth+boot returns GAM boot silently", {
  expect_silent(
    r <- RMSTpowerBoost:::.rmst_route(FALSE, TRUE, "region", "additive", "boot")
  )
  expect_identical(r$model,    "GAM")
  expect_identical(r$fn_power, RMSTpowerBoost::GAM.power.boot)
})

test_that(".rmst_route no-strata+analytical returns linear analytical", {
  r <- RMSTpowerBoost:::.rmst_route(FALSE, FALSE, NULL, "additive", "analytical")
  expect_identical(r$model,    "linear")
  expect_identical(r$fn_power, RMSTpowerBoost::linear.power.analytical)
  expect_identical(r$fn_ss,    RMSTpowerBoost::linear.ss.analytical)
})

test_that(".rmst_route no-strata+boot returns linear boot", {
  r <- RMSTpowerBoost:::.rmst_route(FALSE, FALSE, NULL, "additive", "boot")
  expect_identical(r$model,    "linear")
  expect_identical(r$fn_power, RMSTpowerBoost::linear.power.boot)
  expect_identical(r$fn_ss,    RMSTpowerBoost::linear.ss.boot)
})

test_that(".rmst_route strata+additive+analytical returns additive analytical", {
  r <- RMSTpowerBoost:::.rmst_route(FALSE, FALSE, "region", "additive", "analytical")
  expect_identical(r$model,    "additive")
  expect_identical(r$fn_power, RMSTpowerBoost::additive.power.analytical)
  expect_identical(r$fn_ss,    RMSTpowerBoost::additive.ss.analytical)
})

test_that(".rmst_route strata+additive+boot returns GAM boot", {
  r <- RMSTpowerBoost:::.rmst_route(FALSE, FALSE, "region", "additive", "boot")
  expect_identical(r$model,    "GAM")
  expect_identical(r$fn_power, RMSTpowerBoost::GAM.power.boot)
})

test_that(".rmst_route strata+multiplicative+analytical returns MS analytical", {
  r <- RMSTpowerBoost:::.rmst_route(FALSE, FALSE, "region", "multiplicative", "analytical")
  expect_identical(r$model,    "multiplicative")
  expect_identical(r$fn_power, RMSTpowerBoost::MS.power.analytical)
  expect_identical(r$fn_ss,    RMSTpowerBoost::MS.ss.analytical)
})

test_that(".rmst_route strata+multiplicative+boot returns MS boot", {
  r <- RMSTpowerBoost:::.rmst_route(FALSE, FALSE, "region", "multiplicative", "boot")
  expect_identical(r$model,    "multiplicative")
  expect_identical(r$fn_power, RMSTpowerBoost::MS.power.boot)
  expect_identical(r$fn_ss,    RMSTpowerBoost::MS.ss.boot)
})

# ── .rmst_model_label ─────────────────────────────────────────────────────────

test_that(".rmst_model_label produces correct labels for all model/method pairs", {
  mk <- function(m, mt) list(model = m, method = mt)
  expect_true(grepl("Linear IPCW",         RMSTpowerBoost:::.rmst_model_label(mk("linear",         "analytical"))))
  expect_true(grepl("Analytical",          RMSTpowerBoost:::.rmst_model_label(mk("linear",         "analytical"))))
  expect_true(grepl("Bootstrap",           RMSTpowerBoost:::.rmst_model_label(mk("linear",         "boot"))))
  expect_true(grepl("Additive Stratified", RMSTpowerBoost:::.rmst_model_label(mk("additive",       "boot"))))
  expect_true(grepl("Multiplicative",      RMSTpowerBoost:::.rmst_model_label(mk("multiplicative", "analytical"))))
  expect_true(grepl("GAM",                 RMSTpowerBoost:::.rmst_model_label(mk("GAM",            "boot"))))
  expect_true(grepl("Dependent Censoring", RMSTpowerBoost:::.rmst_model_label(mk("DC",             "analytical"))))
  expect_true(grepl("my_model",            RMSTpowerBoost:::.rmst_model_label(mk("my_model",       "boot"))))
})

# ── rmst.power ────────────────────────────────────────────────────────────────

test_that("rmst.power linear analytical returns rmst_power with 5-element structure", {
  r <- rmst.power(survival::Surv(time, status) ~ age,
                  data = pilot, arm = "arm",
                  sample_sizes = c(30L, 50L), L = 5)
  expect_s3_class(r, "rmst_power")
  expect_named(r, c("results_data", "results_plot", "results_summary",
                    "model_output", ".meta"))
  expect_s3_class(r$results_data, "data.frame")
  expect_identical(r$.meta$model,  "linear")
  expect_identical(r$.meta$method, "analytical")
  expect_identical(r$.meta$truncation, 5)
  expect_identical(r$.meta$n_col,  "N_per_Arm")
})

test_that("rmst.power linear boot uses n_sim and returns correct class", {
  r <- rmst.power(survival::Surv(time, status) ~ age,
                  data = pilot, arm = "arm",
                  sample_sizes = c(20L), L = 5,
                  type = "boot", n_sim = 2L)
  expect_s3_class(r, "rmst_power")
  expect_identical(r$.meta$method, "boot")
})

test_that("rmst.power strata as formula sets strata metadata", {
  r <- rmst.power(survival::Surv(time, status) ~ age,
                  data = pilot_s, arm = "arm",
                  sample_sizes = c(30L), L = 5,
                  strata = ~region,
                  strata_type = "multiplicative")
  expect_identical(r$.meta$strata, "region")
})

test_that("rmst.power strata as character string sets strata metadata", {
  r <- rmst.power(survival::Surv(time, status) ~ age,
                  data = pilot_s, arm = "arm",
                  sample_sizes = c(30L), L = 5,
                  strata = "region",
                  strata_type = "multiplicative")
  expect_identical(r$.meta$strata, "region")
})

test_that("rmst.power multiplicative strata sets N_per_Stratum n_col", {
  r <- rmst.power(survival::Surv(time, status) ~ age,
                  data = pilot_s, arm = "arm",
                  sample_sizes = c(30L), L = 5,
                  strata = "region", strata_type = "multiplicative")
  expect_identical(r$.meta$n_col, "N_per_Stratum")
})

test_that("rmst.power DC boot->analytical emits switching message", {
  expect_message(
    r <- rmst.power(survival::Surv(time, status) ~ age,
                    data = dc_dat, arm = "arm",
                    sample_sizes = c(100L), L = 8,
                    dep_cens = TRUE, type = "boot"),
    "switching type"
  )
  expect_identical(r$.meta$model,          "DC")
  expect_identical(r$.meta$method,         "analytical")
  expect_identical(r$.meta$type_requested, "boot")
})

test_that("rmst.power errors on formula without Surv() LHS", {
  expect_error(
    rmst.power(time ~ age, data = pilot, arm = "arm",
               sample_sizes = c(30L), L = 5),
    "Surv"
  )
})

# ── rmst.ss ───────────────────────────────────────────────────────────────────

test_that("rmst.ss linear analytical returns rmst_ss with 5-element structure", {
  r <- rmst.ss(survival::Surv(time, status) ~ age,
               data = pilot, arm = "arm",
               target_power = 0.02, L = 5,
               n_start = 20L, n_step = 10L, max_n = 60L)
  expect_s3_class(r, "rmst_ss")
  expect_named(r, c("results_data", "results_plot", "results_summary",
                    "model_output", ".meta"))
  expect_identical(r$.meta$model,  "linear")
  expect_identical(r$.meta$method, "analytical")
  expect_identical(r$.meta$n_col,  "N_per_Arm")
})

test_that("rmst.ss DC boot->analytical emits switching message", {
  expect_message(
    r <- rmst.ss(survival::Surv(time, status) ~ age,
                 data = dc_dat, arm = "arm",
                 target_power = 0.50, L = 8,
                 dep_cens = TRUE, type = "boot",
                 n_start = 50L, n_step = 50L, max_n = 300L),
    "switching type"
  )
  expect_identical(r$.meta$model, "DC")
})

# ── print.rmst_power ──────────────────────────────────────────────────────────

test_that("print.rmst_power prints without error and returns x invisibly", {
  r   <- rmst.power(survival::Surv(time, status) ~ age,
                    data = pilot, arm = "arm",
                    sample_sizes = c(30L, 50L), L = 5)
  out <- capture.output(ret <- print(r))
  expect_identical(ret, r)
  expect_true(any(grepl("RMST Power Analysis", out)))
  expect_true(any(grepl("Truncation time",     out)))
  expect_true(any(grepl("Power Results",       out)))
})

test_that("print.rmst_power shows strata line when strata is non-NULL", {
  r   <- rmst.power(survival::Surv(time, status) ~ age,
                    data = pilot_s, arm = "arm",
                    sample_sizes = c(30L), L = 5,
                    strata = "region",
                    strata_type = "multiplicative")
  out <- capture.output(print(r))
  expect_true(any(grepl("Strata", out)))
})

test_that("print.rmst_power shows Note when type was auto-switched", {
  r   <- suppressMessages(
    rmst.power(survival::Surv(time, status) ~ age,
               data = dc_dat, arm = "arm",
               sample_sizes = c(100L), L = 8,
               dep_cens = TRUE, type = "boot")
  )
  out <- capture.output(print(r))
  expect_true(any(grepl("Note", out)))
})

# ── print.rmst_ss ─────────────────────────────────────────────────────────────

test_that("print.rmst_ss prints without error and returns x invisibly", {
  r   <- rmst.ss(survival::Surv(time, status) ~ age,
                 data = pilot, arm = "arm",
                 target_power = 0.02, L = 5,
                 n_start = 20L, n_step = 10L, max_n = 60L)
  out <- capture.output(ret <- print(r))
  expect_identical(ret, r)
  expect_true(any(grepl("RMST Sample Size", out)))
  expect_true(any(grepl("Truncation time",  out)))
})

test_that("print.rmst_ss shows Note when type was auto-switched", {
  r   <- suppressMessages(
    rmst.ss(survival::Surv(time, status) ~ age,
            data = dc_dat, arm = "arm",
            target_power = 0.50, L = 8,
            dep_cens = TRUE, type = "boot",
            n_start = 50L, n_step = 50L, max_n = 300L)
  )
  out <- capture.output(print(r))
  expect_true(any(grepl("Note", out)))
})

# ── summary.rmst_power ────────────────────────────────────────────────────────

test_that("summary.rmst_power returns summary.rmst_power with ten named fields", {
  r <- rmst.power(survival::Surv(time, status) ~ age,
                  data = pilot, arm = "arm",
                  sample_sizes = c(30L), L = 5)
  s <- summary(r)
  expect_s3_class(s, "summary.rmst_power")
  expect_true(all(c("model_info", "power_results", "effect_summary",
                    "coefficient_table", "treatment_effect",
                    "arm_specific_rmst", "variance_components",
                    "censoring_weights", "diagnostics",
                    "simulation_draws") %in% names(s)))
  expect_s3_class(s$model_info, "data.frame")
  expect_identical(s$model_info$Model, "linear")
})

test_that("summary.rmst_power strata meta shows non-None Strata field", {
  r <- rmst.power(survival::Surv(time, status) ~ age,
                  data = pilot_s, arm = "arm",
                  sample_sizes = c(30L), L = 5,
                  strata = "region",
                  strata_type = "multiplicative")
  s <- summary(r)
  expect_identical(s$model_info$Strata, "region")
})

# ── summary.rmst_ss ───────────────────────────────────────────────────────────

test_that("summary.rmst_ss returns summary.rmst_ss with required fields", {
  r <- rmst.ss(survival::Surv(time, status) ~ age,
               data = pilot, arm = "arm",
               target_power = 0.02, L = 5,
               n_start = 20L, n_step = 10L, max_n = 60L)
  s <- summary(r)
  expect_s3_class(s, "summary.rmst_ss")
  expect_true(all(c("model_info", "ss_results", "effect_summary",
                    "treatment_effect", "diagnostics") %in% names(s)))
})

# ── print.summary.rmst_power ──────────────────────────────────────────────────

test_that("print.summary.rmst_power (analytical) covers main sections", {
  r   <- rmst.power(survival::Surv(time, status) ~ age,
                    data = pilot, arm = "arm",
                    sample_sizes = c(30L), L = 5)
  s   <- summary(r)
  out <- capture.output(ret <- print(s))
  expect_identical(ret, s)
  expect_true(any(grepl("Model Information", out)))
  expect_true(any(grepl("Power Results",     out)))
})

test_that("print.summary.rmst_power (boot) shows NULL coefficient_table note", {
  r   <- rmst.power(survival::Surv(time, status) ~ age,
                    data = pilot, arm = "arm",
                    sample_sizes = c(20L), L = 5,
                    type = "boot", n_sim = 2L)
  s   <- summary(r)
  out <- capture.output(print(s))
  expect_true(any(grepl("not available", out)))
})

test_that("print.summary.rmst_power mock covers all optional output sections", {
  mock <- structure(
    list(
      model_info = data.frame(
        Model = "linear", Method = "analytical", Truncation_time = 10,
        Alpha = 0.05, Arm_variable = "arm", Strata = "None",
        Dep_censoring = FALSE, stringsAsFactors = FALSE
      ),
      power_results  = data.frame(N_per_Arm = c(50L, 100L), Power = c(0.5, 0.8)),
      effect_summary = data.frame(estimand = "RMST diff", estimate = 1.2),
      coefficient_table = data.frame(term = "arm", estimate = 1.2, p.value = 0.03),
      treatment_effect  = data.frame(estimand = "RMST diff", estimate = 1.2,
                                     std_error = 0.5),
      arm_specific_rmst = data.frame(arm = c(0L, 1L), rmst = c(8.0, 9.2),
                                     std_error = c(NA_real_, NA_real_)),
      variance_components = list(
        se_effect_n1 = 0.5,
        A_hat        = matrix(1, 1L, 1L),
        B_hat        = matrix(1, 1L, 1L),
        V_hat_n      = matrix(1, 1L, 1L)
      ),
      censoring_weights = list(
        raw_summary     = stats::quantile(c(0.5, 0.8, 1.0), c(0.25, 0.5, 0.75)),
        cap_value       = 2.0,
        capped_fraction = 0.05
      ),
      diagnostics = list(n_used = 100L, n_events = 75L, n_sim = 10L,
                         convergence_ok = TRUE, singular_flag = FALSE),
      simulation_draws = data.frame(
        replicate = 1:3, estimate = c(1.1, 1.2, 1.3),
        std_error = NA_real_, p_value = c(0.01, 0.02, 0.03)
      )
    ),
    class = "summary.rmst_power"
  )
  out <- capture.output(ret <- print(mock))
  expect_identical(ret, mock)
  expect_true(any(grepl("Coefficient Table",  out)))
  expect_true(any(grepl("Arm-Specific RMST",  out)))
  expect_true(any(grepl("Variance Components", out)))
  expect_true(any(grepl("Censoring Weights",  out)))
  expect_true(any(grepl("Simulation Draws",   out)))
  expect_true(any(grepl("Note: SE",           out)))  # NA std_error note
})

# ── print.summary.rmst_ss ─────────────────────────────────────────────────────

test_that("print.summary.rmst_ss prints all sections and returns x invisibly", {
  r   <- rmst.ss(survival::Surv(time, status) ~ age,
                 data = pilot, arm = "arm",
                 target_power = 0.02, L = 5,
                 n_start = 20L, n_step = 10L, max_n = 60L)
  s   <- summary(r)
  out <- capture.output(ret <- print(s))
  expect_identical(ret, s)
  expect_true(any(grepl("Sample Size Result", out)))
})

# ── plot.rmst_power / plot.rmst_ss ────────────────────────────────────────────

test_that("plot.rmst_power returns ggplot invisibly without error", {
  r   <- rmst.power(survival::Surv(time, status) ~ age,
                    data = pilot, arm = "arm",
                    sample_sizes = c(30L, 50L), L = 5)
  ret <- expect_no_error(plot(r))
  expect_true(inherits(ret, "ggplot"))
})

test_that("plot.rmst_ss returns ggplot invisibly without error", {
  r   <- rmst.ss(survival::Surv(time, status) ~ age,
                 data = pilot, arm = "arm",
                 target_power = 0.02, L = 5,
                 n_start = 20L, n_step = 10L, max_n = 60L)
  ret <- expect_no_error(plot(r))
  expect_true(inherits(ret, "ggplot"))
})

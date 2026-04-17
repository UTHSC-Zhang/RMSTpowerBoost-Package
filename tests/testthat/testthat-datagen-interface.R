# Tests for R/datagen-interface.R
# Covers: covar_continuous, covar_binary, covar_categorical,
#         recipe_quick_ph, rmst.sim, and their S3 methods.

# ── covar_continuous ──────────────────────────────────────────────────────────

test_that("covar_continuous returns correct list structure", {
  out <- covar_continuous("age", dist = "normal", mean = 50, sd = 10)
  expect_identical(out$name,   "age")
  expect_identical(out$type,   "continuous")
  expect_identical(out$dist,   "normal")
  expect_identical(out$params, list(mean = 50, sd = 10))
})

test_that("covar_continuous default dist is 'normal' with empty params", {
  out <- covar_continuous("x")
  expect_identical(out$dist,   "normal")
  expect_identical(out$params, list())
})

test_that("covar_continuous accepts all valid distributions", {
  for (d in c("normal", "lognormal", "gamma", "weibull", "uniform", "beta", "t")) {
    out <- covar_continuous("x", dist = d)
    expect_identical(out$dist, d)
  }
})

test_that("covar_continuous errors on invalid distribution via match.arg", {
  expect_error(covar_continuous("x", dist = "binomial"))
  expect_error(covar_continuous("x", dist = "poisson"))
})

# ── covar_binary ──────────────────────────────────────────────────────────────

test_that("covar_binary returns correct list structure", {
  out <- covar_binary("female", p = 0.52)
  expect_identical(out$name,      "female")
  expect_identical(out$type,      "categorical")
  expect_identical(out$dist,      "bernoulli")
  expect_equal(out$params$p,      0.52)
})

test_that("covar_binary default p = 0.5", {
  out <- covar_binary("male")
  expect_equal(out$params$p, 0.5)
})

test_that("covar_binary accepts boundary values p=0 and p=1", {
  expect_no_error(covar_binary("x", p = 0))
  expect_no_error(covar_binary("x", p = 1))
})

test_that("covar_binary errors on out-of-range p", {
  expect_error(covar_binary("x", p = -0.1),  "'p' must be")
  expect_error(covar_binary("x", p =  1.01), "'p' must be")
})

test_that("covar_binary errors on non-scalar p", {
  expect_error(covar_binary("x", p = c(0.3, 0.7)), "'p' must be")
})

test_that("covar_binary errors on non-numeric p", {
  expect_error(covar_binary("x", p = "0.5"), "'p' must be")
})

# ── covar_categorical ─────────────────────────────────────────────────────────

test_that("covar_categorical returns correct structure with auto-labels", {
  out <- covar_categorical("stage", probs = c(0.4, 0.35, 0.25))
  expect_identical(out$name,          "stage")
  expect_identical(out$type,          "categorical")
  expect_identical(out$dist,          "categorical")
  expect_equal(out$params$prob,       c(0.4, 0.35, 0.25))
  expect_identical(out$params$labels, c("1", "2", "3"))
})

test_that("covar_categorical uses custom labels when supplied", {
  out <- covar_categorical("grade", probs = c(0.5, 0.3, 0.2),
                           labels = c("I", "II", "III"))
  expect_identical(out$params$labels, c("I", "II", "III"))
})

test_that("covar_categorical errors when probs do not sum to 1", {
  expect_error(covar_categorical("x", probs = c(0.3, 0.3)), "'probs' must sum")
})

test_that("covar_categorical errors on negative probabilities", {
  expect_error(covar_categorical("x", probs = c(0.5, -0.5, 1.0)), "non-negative")
})

test_that("covar_categorical errors when labels length mismatches probs", {
  expect_error(
    covar_categorical("x", probs = c(0.5, 0.5), labels = c("A")),
    "same length"
  )
})

# ── recipe_quick_ph ───────────────────────────────────────────────────────────

test_that("recipe_quick_ph ph_exponential returns valid validated recipe", {
  r <- recipe_quick_ph(50L, "ph_exponential",
                       baseline = list(rate = 0.2),
                       treat_effect = -0.3, seed = 1L)
  expect_true(is.list(r))
  expect_identical(r$event_time$model, "ph_exponential")
  expect_identical(r$n, 50L)
})

test_that("recipe_quick_ph ph_exponential simulates correct nrow", {
  r   <- recipe_quick_ph(50L, "ph_exponential",
                         baseline = list(rate = 0.2),
                         treat_effect = -0.3, seed = 2L)
  dat <- simulate_from_recipe(r, seed = 2L)
  expect_true(is.data.frame(dat))
  expect_equal(nrow(dat), 50L)
})

test_that("recipe_quick_ph ph_weibull includes covariate column in output", {
  r <- recipe_quick_ph(60L, "ph_weibull",
                       baseline = list(shape = 1.5, scale = 10),
                       treat_effect = -0.4,
                       covariates = list(covar_continuous("age", mean = 60, sd = 10)),
                       seed = 3L)
  expect_identical(r$event_time$model, "ph_weibull")
  dat <- simulate_from_recipe(r, seed = 3L)
  expect_true("age" %in% names(dat))
})

test_that("recipe_quick_ph ph_gompertz produces finite event times", {
  r   <- recipe_quick_ph(40L, "ph_gompertz",
                         baseline = list(shape = 0.05, rate = 0.1),
                         treat_effect = -0.2, seed = 4L)
  expect_identical(r$event_time$model, "ph_gompertz")
  dat <- simulate_from_recipe(r, seed = 4L)
  expect_true(all(is.finite(dat$time)))
  expect_equal(nrow(dat), 40L)
})

test_that("recipe_quick_ph cox_pwexp generates data of correct size", {
  r   <- recipe_quick_ph(50L, "cox_pwexp",
                         baseline = list(rates = c(0.1, 0.2), cuts = 5),
                         treat_effect = -0.3, seed = 5L)
  dat <- simulate_from_recipe(r, seed = 5L)
  expect_equal(nrow(dat), 50L)
})

test_that("recipe_quick_ph errors on invalid model name", {
  expect_error(
    recipe_quick_ph(50L, "ph_invalid",
                    baseline = list(rate = 0.2), treat_effect = 0)
  )
})

test_that("recipe_quick_ph passes binary covariate to generated data", {
  r   <- recipe_quick_ph(40L, "ph_exponential",
                         baseline = list(rate = 0.2),
                         treat_effect = 0,
                         covariates = list(covar_binary("sex")),
                         seed = 6L)
  dat <- simulate_from_recipe(r, seed = 6L)
  expect_true("sex" %in% names(dat))
})

test_that("recipe_quick_ph respects non-default allocation", {
  r <- recipe_quick_ph(60L, "ph_exponential",
                       baseline = list(rate = 0.1), treat_effect = 0,
                       allocation = "1:2", seed = 7L)
  expect_identical(r$treatment$allocation, "1:2")
})

# ── rmst.sim ──────────────────────────────────────────────────────────────────

test_that("rmst.sim aft_lognormal returns rmst_simdata with correct attrs", {
  df <- rmst.sim(n = 80L, model = "aft_lognormal",
                 baseline = list(mu = 2.2, sigma = 0.5),
                 treat_effect = -0.3, L = 12, seed = 42L)
  expect_s3_class(df, "rmst_simdata")
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 80L)
  expect_identical(attr(df, "L"), 12)
  expect_false(is.null(attr(df, "recipe")))
  expect_true(all(c("time", "status", "arm") %in% names(df)))
})

test_that("rmst.sim aft_weibull (recipe_quick_aft path) works", {
  df <- rmst.sim(n = 60L, model = "aft_weibull",
                 baseline = list(shape = 1.5, scale = 10),
                 treat_effect = -0.2, seed = 10L)
  expect_s3_class(df, "rmst_simdata")
  expect_equal(nrow(df), 60L)
  expect_identical(attr(df, "recipe")$event_time$model, "aft_weibull")
})

test_that("rmst.sim aft_loglogistic (direct validate_recipe path) works", {
  df <- rmst.sim(n = 60L, model = "aft_loglogistic",
                 baseline = list(scale = 8, shape = 1.2),
                 treat_effect = -0.2, seed = 20L)
  expect_s3_class(df, "rmst_simdata")
  expect_identical(attr(df, "recipe")$event_time$model, "aft_loglogistic")
})

test_that("rmst.sim ph_exponential (recipe_quick_ph path) works", {
  df <- rmst.sim(n = 60L, model = "ph_exponential",
                 baseline = list(rate = 0.15),
                 treat_effect = -0.4, seed = 30L)
  expect_s3_class(df, "rmst_simdata")
  expect_equal(nrow(df), 60L)
})

test_that("rmst.sim ph_weibull works", {
  df <- rmst.sim(n = 60L, model = "ph_weibull",
                 baseline = list(shape = 1.5, scale = 10),
                 treat_effect = -0.3, seed = 31L)
  expect_s3_class(df, "rmst_simdata")
})

test_that("rmst.sim ph_gompertz works", {
  df <- rmst.sim(n = 60L, model = "ph_gompertz",
                 baseline = list(shape = 0.05, rate = 0.1),
                 treat_effect = -0.2, seed = 32L)
  expect_s3_class(df, "rmst_simdata")
})

test_that("rmst.sim cox_pwexp works", {
  df <- rmst.sim(n = 60L, model = "cox_pwexp",
                 baseline = list(rates = c(0.1, 0.2), cuts = 5),
                 treat_effect = -0.3, seed = 33L)
  expect_s3_class(df, "rmst_simdata")
})

test_that("rmst.sim includes covariate columns from all helper types", {
  df <- rmst.sim(n = 80L, model = "aft_lognormal",
                 baseline = list(mu = 2.2, sigma = 0.5),
                 treat_effect = -0.3,
                 covariates = list(
                   covar_continuous("age",   mean = 60, sd = 10),
                   covar_binary("female",    p = 0.5),
                   covar_categorical("stage", probs = c(0.4, 0.35, 0.25))
                 ),
                 seed = 50L)
  expect_true(all(c("age", "female", "stage") %in% names(df)))
})

test_that("rmst.sim seed produces reproducible output", {
  df1 <- rmst.sim(n = 50L, model = "aft_lognormal",
                  baseline = list(mu = 2, sigma = 0.5), seed = 99L)
  df2 <- rmst.sim(n = 50L, model = "aft_lognormal",
                  baseline = list(mu = 2, sigma = 0.5), seed = 99L)
  expect_identical(df1$time,   df2$time)
  expect_identical(df1$status, df2$status)
})

test_that("rmst.sim L=NULL stores NULL attribute", {
  df <- rmst.sim(n = 50L, model = "aft_lognormal",
                 baseline = list(mu = 2, sigma = 0.5), seed = 1L)
  expect_null(attr(df, "L"))
})

test_that("rmst.sim achieved_censoring attribute is numeric in [0, 1]", {
  df <- rmst.sim(n = 80L, model = "aft_lognormal",
                 baseline = list(mu = 2.2, sigma = 0.5),
                 target_censoring = 0.25, seed = 77L)
  ac <- attr(df, "achieved_censoring")
  expect_true(is.numeric(ac) && length(ac) == 1L)
  expect_true(ac >= 0 && ac <= 1)
})

test_that("rmst.sim errors on unrecognised model name", {
  expect_error(
    rmst.sim(n = 50L, model = "wrong_model",
             baseline = list(mu = 2, sigma = 0.5))
  )
})

# ── print.rmst_simdata ────────────────────────────────────────────────────────

test_that("print.rmst_simdata prints without error and returns x invisibly", {
  df  <- rmst.sim(n = 60L, model = "aft_lognormal",
                  baseline = list(mu = 2.2, sigma = 0.5),
                  treat_effect = -0.3, L = 10, seed = 1L)
  out <- capture.output(ret <- print(df))
  expect_identical(ret, df)
  expect_true(any(grepl("Simulated RMST Dataset", out)))
  expect_true(any(grepl("Truncation time",        out)))
  expect_true(any(grepl("Arm 0",                  out)))
})

test_that("print.rmst_simdata skips Arm breakdown and truncation when absent", {
  df      <- rmst.sim(n = 40L, model = "aft_lognormal",
                      baseline = list(mu = 2, sigma = 0.5), seed = 2L)
  df_bare <- structure(
    as.data.frame(df),
    class             = c("rmst_simdata", "data.frame"),
    recipe            = attr(df, "recipe"),
    L                 = NULL,
    achieved_censoring = attr(df, "achieved_censoring")
  )
  df_bare$arm <- NULL
  out <- capture.output(print(df_bare))
  expect_true( any(grepl("N \\(total\\)", out)))
  expect_false(any(grepl("Arm 0",          out)))
  expect_false(any(grepl("Truncation time", out)))
})

# ── summary.rmst_simdata ──────────────────────────────────────────────────────

test_that("summary.rmst_simdata returns summary.rmst_simdata with six named fields", {
  df <- rmst.sim(n = 80L, model = "aft_lognormal",
                 baseline = list(mu = 2.2, sigma = 0.5),
                 treat_effect = -0.3,
                 covariates = list(covar_continuous("age")),
                 seed = 5L)
  s <- summary(df)
  expect_s3_class(s, "summary.rmst_simdata")
  expect_true(is.list(s))
  expect_true(all(c("header", "effects", "treatment",
                    "censoring", "covariates", "baseline") %in% names(s)))
  expect_true(is.data.frame(s$header))
})

test_that("summary.rmst_simdata works when there are no covariates", {
  df <- rmst.sim(n = 50L, model = "ph_weibull",
                 baseline = list(shape = 1.5, scale = 10),
                 treat_effect = -0.3, seed = 7L)
  s  <- summary(df)
  expect_s3_class(s, "summary.rmst_simdata")
  expect_true(is.null(s$covariates) || is.data.frame(s$covariates))
})

test_that("print.summary.rmst_simdata prints all sections and returns x invisibly", {
  df  <- rmst.sim(n = 60L, model = "aft_lognormal",
                  baseline = list(mu = 2.2, sigma = 0.5),
                  treat_effect = -0.3, seed = 6L)
  s   <- summary(df)
  out <- capture.output(ret <- print(s))
  expect_identical(ret, s)
  expect_true(any(grepl("Simulated Dataset Summary", out)))
})

# ── .rmst_sim_model_label helper ──────────────────────────────────────────────

test_that(".rmst_sim_model_label returns friendly strings for all seven models", {
  models <- c("aft_lognormal", "aft_weibull", "aft_loglogistic",
              "ph_exponential", "ph_weibull", "ph_gompertz", "cox_pwexp")
  labels <- vapply(models, RMSTpowerBoost:::.rmst_sim_model_label, character(1L))
  expect_false(any(is.na(labels)))
  expect_identical(labels[["aft_lognormal"]], "AFT Log-Normal")
  expect_identical(labels[["ph_weibull"]],    "PH Weibull")
  expect_identical(labels[["cox_pwexp"]],     "PH Piecewise-Exponential (Cox)")
})

test_that(".rmst_sim_model_label falls back to the raw model string for unknowns", {
  expect_identical(
    RMSTpowerBoost:::.rmst_sim_model_label("unknown_xyz"),
    "unknown_xyz"
  )
})

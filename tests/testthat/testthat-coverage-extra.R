test_that("run_app uses override env var to trigger missing app error", {
  old <- Sys.getenv("RMSTPOWERBOOST_APP_DIR", unset = NA_character_)
  on.exit(Sys.setenv(RMSTPOWERBOOST_APP_DIR = old), add = TRUE)
  Sys.setenv(RMSTPOWERBOOST_APP_DIR = "")
  expect_error(run_app(), "Could not find the app directory")
})

test_that("recipe internals cover model aliases, transforms, and assignments", {
  nm <- RMSTpowerBoost:::.normalize_model_name("ph_exp", baseline = list(rate = 0.2))
  expect_identical(nm$model, "cox_pwexp")
  expect_true(!is.null(nm$baseline$rates))

  nm2 <- RMSTpowerBoost:::.normalize_model_name("ph_pwexp", baseline = list(rates = c(0.1, 0.2), cuts = 3))
  expect_identical(nm2$model, "cox_pwexp")

  x <- c(1, 2, 3)
  expect_equal(RMSTpowerBoost:::.apply_transforms(x, c("center(2)", "scale(2)")), c(-0.5, 0, 0.5))

  cat_def <- list(dist = "categorical", params = list(prob = c(0.5, 0.5), labels = c("A", "B")))
  ord_def <- list(dist = "ordinal", params = list(prob = c(0.2, 0.8), labels = c("L", "H")))
  cat_vals <- RMSTpowerBoost:::.rdraw(cat_def, 10)
  ord_vals <- RMSTpowerBoost:::.rdraw(ord_def, 10)
  expect_true(is.factor(cat_vals))
  expect_true(is.ordered(ord_vals))

  X <- data.frame(stratum = factor(rep(c("S1", "S2"), each = 5)), x = rnorm(10))
  tr <- list(assignment = "stratified", allocation = "1:1", stratify_by = "stratum")
  arm_s <- RMSTpowerBoost:::.assign_treatment(10, X, tr)
  expect_length(arm_s, 10)

  tr2 <- list(assignment = "logistic_ps", ps_model = list(formula = "~ x", beta = c(0, 0.5)))
  arm_l <- RMSTpowerBoost:::.assign_treatment(10, X, tr2)
  expect_length(arm_l, 10)
})

test_that("frailty, lp, and censoring branches are exercised", {
  X <- data.frame(x = c(0, 1, 2), g = factor(c("a", "b", "a")))
  effects <- list(intercept = 0.2, treatment = -0.3, formula = "~ x", beta = c(1, 0.5))
  lp <- RMSTpowerBoost:::.build_lp(effects, X, arm = c(0, 1, 0))
  expect_length(lp, 3)

  eta_g <- RMSTpowerBoost:::.apply_frailty(lp, model = "ph_weibull",
                                           frailty = list(type = "gamma", var = 0.2, group = "g"),
                                           groups = X)
  eta_l <- RMSTpowerBoost:::.apply_frailty(lp, model = "ph_weibull",
                                           frailty = list(type = "lognormal", var = 0.2, group = "g"),
                                           groups = X)
  expect_length(eta_g, 3)
  expect_length(eta_l, 3)

  T_event <- c(5, 6, 7)
  expect_equal(RMSTpowerBoost:::.solve_rate_for_target(T_event, target = 0.5, admin_time = 0), 0)
})

test_that("explicit censoring path and pwexp simulator work", {
  set.seed(123)
  covs <- list(
    list(name = "x", type = "continuous", dist = "normal", params = list(mean = 0, sd = 1))
  )
  rec <- list(
    n = 30,
    covariates = list(defs = covs),
    treatment = list(assignment = "randomization", allocation = "1:1"),
    event_time = list(
      model = "ph_exponential",
      baseline = list(rate = 0.2),
      effects = list(intercept = 0, treatment = -0.2, covariates = list(x = 0.1))
    ),
    censoring = list(
      mode = "explicit",
      administrative = list(time = 2),
      random = list(dist = "exponential", params = list(rate = 0.1)),
      dependent = list(formula = "~ x", beta = c(0, 0.5), base = 0.05)
    )
  )
  dat <- simulate_from_recipe(rec, seed = 99)
  expect_true(is.data.frame(dat))
  expect_true(is.finite(attr(dat, "achieved_censoring")))

  t_pw <- RMSTpowerBoost:::.sim_pwexp(rates = c(0.2, 0.05), cuts = 3, lp = rep(0, 5))
  expect_length(t_pw, 5)
  expect_true(all(is.finite(t_pw)))
})

test_that("generate_recipe_sets handles vary grid, templates, and formats", {
  covs <- list(
    list(name = "g", type = "categorical", dist = "categorical",
         params = list(prob = c(0.4, 0.6), labels = c("A", "B")))
  )
  rec <- recipe_quick_aft(
    n = 20,
    model = "aft_lognormal",
    baseline = list(mu = 2.4, sigma = 0.6),
    treat_effect = -0.1,
    covariates = covs,
    target_censoring = 0.1,
    allocation = "1:1"
  )
  out_dir <- file.path(tempdir(), "rmst-gensets-extra")
  man <- generate_recipe_sets(
    rec,
    vary = list(n = c(20, 25), "event_time.effects.treatment" = c(-0.1, -0.2)),
    out_dir = out_dir,
    formats = c("rdata"),
    n_reps = 1,
    seed_base = 10,
    filename_template = "sc{scenario_id}_n{n}_trt{event_time.effects.treatment}_r{rep}"
  )
  expect_true(nrow(man) == 4)
  expect_true(all(file.exists(man$file_rdata)))

  grid <- recipe_grid(rec, list(n = c(10, 11)))
  expect_true(length(grid) == 2)
  expect_true(all(vapply(grid, function(x) x$n, integer(1)) %in% c(10, 11)))
})

test_that("simulation runners cover app helpers with smooth terms", {
  set.seed(321)
  pilot <- data.frame(
    time = stats::rexp(40, rate = 0.2),
    status = stats::rbinom(40, 1, 0.9),
    arm = rep(0:1, length.out = 40),
    age = stats::rnorm(40, 60, 6),
    region = factor(rep(c("A", "B"), each = 20))
  )
  runner <- .get_gam_simulation_runner(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm", strata_var = "region",
    linear_terms = NULL, smooth_terms = "age", L = 3, alpha = 0.1, n_sim = 2, parallel.cores = 1
  )
  sim <- runner(5)
  expect_true(is.list(sim))
  expect_true(all(c("power", "estimates", "std_errors") %in% names(sim)))

  runner_lin <- .get_linear_ipcw_simulation_runner(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    linear_terms = "age", L = 3, alpha = 0.1, n_sim = 2, parallel.cores = 1
  )
  sim_lin <- runner_lin(5)
  expect_true(is.list(sim_lin))
  expect_true(all(c("power", "estimates", "std_errors") %in% names(sim_lin)))
})

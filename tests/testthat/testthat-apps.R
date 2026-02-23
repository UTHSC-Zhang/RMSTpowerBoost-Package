# Tests for app-facing canonical functions

make_pilot_data <- function(n = 60, stratified = FALSE) {
  df <- data.frame(
    time = stats::rexp(n, rate = 0.12),
    status = stats::rbinom(n, 1, 0.9),
    arm = rep(0:1, length.out = n),
    age = stats::rnorm(n, 60, 8)
  )
  if (stratified) {
    if (n %% 4 != 0) {
      n <- n + (4 - (n %% 4))
      df <- data.frame(
        time = stats::rexp(n, rate = 0.12),
        status = stats::rbinom(n, 1, 0.9),
        arm = rep(0:1, length.out = n),
        age = stats::rnorm(n, 60, 8)
      )
    }
    df$region <- factor(rep(c("A", "B"), each = n / 2))
    df$arm <- rep(rep(0:1, each = n / 4), times = 2)
  }
  df
}

test_that("analytical functions return structured outputs", {
  set.seed(101)
  pilot <- make_pilot_data(n = 80, stratified = TRUE)

  lin_power <- linear.power.analytical(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    sample_sizes = c(10, 20), linear_terms = "age", L = 5
  )
  expect_named(lin_power, c("results_data", "results_plot", "results_summary"))
  expect_s3_class(lin_power$results_plot, "ggplot")

  lin_ss <- linear.ss.analytical(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    target_power = 0, linear_terms = "age", L = 5,
    n_start = 10, n_step = 10, max_n_per_arm = 10
  )
  expect_true(is.data.frame(lin_ss$results_data))

  dc_power <- DC.power.analytical(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    sample_sizes = c(10), linear_terms = "age", L = 5
  )
  expect_true(is.data.frame(dc_power$results_data))

  dc_ss <- DC.ss.analytical(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    target_power = 0, linear_terms = "age", L = 5,
    n_start = 10, n_step = 10, max_n_per_arm = 10
  )
  expect_true(is.data.frame(dc_ss$results_data))

  add_power <- additive.power.analytical(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    strata_var = "region", sample_sizes = c(10), linear_terms = "age", L = 5
  )
  expect_true(is.data.frame(add_power$results_data))

  add_ss <- additive.ss.analytical(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    strata_var = "region", target_power = 0, linear_terms = "age", L = 5,
    n_start = 10, n_step = 10, max_n_per_arm = 10
  )
  expect_true(is.data.frame(add_ss$results_data))

  ms_power <- MS.power.analytical(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    strata_var = "region", sample_sizes = c(10), linear_terms = "age", L = 5
  )
  expect_true(is.data.frame(ms_power$results_data))

  ms_ss <- MS.ss.analytical(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    strata_var = "region", target_power = 0, linear_terms = "age", L = 5,
    n_start = 10, n_step = 10, max_n_per_arm = 10
  )
  expect_true(is.data.frame(ms_ss$results_data))
})

test_that("bootstrap functions return expected structures", {
  set.seed(202)
  pilot_strat <- make_pilot_data(n = 60, stratified = TRUE)
  pilot <- make_pilot_data(n = 60, stratified = FALSE)

  add_power <- GAM.power.boot(
    pilot_data = pilot_strat, time_var = "time", status_var = "status", arm_var = "arm",
    strata_var = "region", sample_sizes = c(8), linear_terms = "age",
    L = 5, n_sim = 2, alpha = 0.1, parallel.cores = 1
  )
  expect_named(add_power, c("results_data", "results_plot", "results_summary"))

  add_ss <- GAM.ss.boot(
    pilot_data = pilot_strat, time_var = "time", status_var = "status", arm_var = "arm",
    strata_var = "region", target_power = 0, linear_terms = "age",
    L = 5, n_sim = 2, alpha = 0.1, parallel.cores = 1,
    patience = 1, n_start = 8, n_step = 5, max_n_per_arm = 8
  )
  expect_true(is.data.frame(add_ss$results_data))

  lin_power <- linear.power.boot(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    sample_sizes = c(8), linear_terms = "age", L = 5, n_sim = 2, alpha = 0.1
  )
  expect_true(is.data.frame(lin_power$results_data))

  lin_ss <- linear.ss.boot(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    target_power = 0, linear_terms = "age", L = 5, n_sim = 2, alpha = 0.1,
    patience = 1, n_start = 8, n_step = 5, max_n_per_arm = 8
  )
  expect_true(is.data.frame(lin_ss$results_data))

  ms_power <- MS.power.boot(
    pilot_data = pilot_strat, time_var = "time", status_var = "status", arm_var = "arm",
    strata_var = "region", sample_sizes = c(8), linear_terms = "age",
    L = 5, n_sim = 2, alpha = 0.1, parallel.cores = 1
  )
  expect_true(is.data.frame(ms_power$results_data))

  ms_ss <- MS.ss.boot(
    pilot_data = pilot_strat, time_var = "time", status_var = "status", arm_var = "arm",
    strata_var = "region", target_power = 0, linear_terms = "age",
    L = 5, n_sim = 2, alpha = 0.1, parallel.cores = 1,
    patience = 1, n_start = 8, n_step = 5, max_n_per_arm = 8
  )
  expect_true(is.data.frame(ms_ss$results_data))
})

test_that("run_app calls shiny::runApp with the app directory", {
  testthat::skip_if_not_installed("shiny")

  shiny_ns <- asNamespace("shiny")
  orig_run_app <- get("runApp", envir = shiny_ns)
  called <- FALSE
  assignInNamespace("runApp", function(appDir, ...) {
    called <<- TRUE
    appDir
  }, ns = "shiny")
  on.exit(assignInNamespace("runApp", orig_run_app, ns = "shiny"), add = TRUE)

  res <- testthat::with_mocked_bindings(
    run_app(),
    .ensure_app_dependencies = function(...) invisible(TRUE),
    .package = "RMSTpowerBoost"
  )
  expect_true(called)
  expect_true(dir.exists(res))
})

test_that("ensure_app_dependencies cancels when user declines install", {
  expect_error(
    testthat::with_mocked_bindings(
      RMSTpowerBoost:::.ensure_app_dependencies(install_missing = TRUE),
      .get_missing_app_dependencies = function() c("shiny", "mice"),
      .prompt_install_app_dependencies = function(missing) FALSE,
      .package = "RMSTpowerBoost"
    ),
    "App launch canceled"
  )
})

test_that("ensure_app_dependencies installs only missing packages", {
  calls <- 0L
  installed <- NULL
  testthat::with_mocked_bindings(
    {
      expect_silent(RMSTpowerBoost:::.ensure_app_dependencies(install_missing = TRUE))
      expect_equal(calls, 2L)
      expect_identical(installed, c("shiny", "mice"))
    },
    .get_missing_app_dependencies = function() {
      calls <<- calls + 1L
      if (calls == 1L) c("shiny", "mice") else character(0)
    },
    .prompt_install_app_dependencies = function(missing) TRUE,
    .install_app_dependencies = function(missing, repos) {
      installed <<- missing
      invisible(TRUE)
    },
    .package = "RMSTpowerBoost"
  )
})

test_that("ensure_app_dependencies skips install when none missing", {
  testthat::with_mocked_bindings(
    {
      expect_silent(RMSTpowerBoost:::.ensure_app_dependencies(install_missing = TRUE))
    },
    .get_missing_app_dependencies = function() character(0),
    .install_app_dependencies = function(missing, repos) stop("should not install"),
    .package = "RMSTpowerBoost"
  )
})

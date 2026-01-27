# Tests for app-facing helpers and wrappers

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

test_that("app analytical helpers return structured outputs", {
  set.seed(101)
  pilot <- make_pilot_data(n = 80, stratified = TRUE)

  lin_power <- linear.power.analytical.app(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    sample_sizes = c(10, 20), linear_terms = "age", L = 5
  )
  expect_named(lin_power, c("results_data", "results_plot", "results_summary"))
  expect_s3_class(lin_power$results_plot, "ggplot")

  lin_ss <- linear.ss.analytical.app(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    target_power = 0, linear_terms = "age", L = 5,
    n_start = 10, n_step = 10, max_n_per_arm = 10
  )
  expect_true(is.data.frame(lin_ss$results_data))

  dc_power <- DC.power.analytical.app(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    dep_cens_status_var = "status", sample_sizes = c(10),
    linear_terms = "age", L = 5
  )
  expect_true(is.data.frame(dc_power$results_data))

  dc_ss <- DC.ss.analytical.app(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    dep_cens_status_var = "status", target_power = 0,
    linear_terms = "age", L = 5,
    n_start = 10, n_step = 10, max_n_per_arm = 10
  )
  expect_true(is.data.frame(dc_ss$results_data))

  add_power <- additive.power.analytical.app(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    strata_var = "region", sample_sizes = c(10), linear_terms = "age", L = 5
  )
  expect_true(is.data.frame(add_power$results_data))

  add_ss <- additive.ss.analytical.app(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    strata_var = "region", target_power = 0, linear_terms = "age", L = 5,
    n_start = 10, n_step = 10, max_n_per_arm = 10
  )
  expect_true(is.data.frame(add_ss$results_data))

  ms_power <- MS.power.analytical.app(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    strata_var = "region", sample_sizes = c(10), linear_terms = "age", L = 5
  )
  expect_true(is.data.frame(ms_power$results_data))

  ms_ss <- MS.ss.analytical.app(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    strata_var = "region", target_power = 0, linear_terms = "age", L = 5,
    n_start = 10, n_step = 10, max_n_per_arm = 10
  )
  expect_true(is.data.frame(ms_ss$results_data))
})

test_that("app bootstrap wrappers return expected structures", {
  set.seed(202)
  pilot_strat <- make_pilot_data(n = 60, stratified = TRUE)
  pilot <- make_pilot_data(n = 60, stratified = FALSE)

  add_power <- additive.power.boot.app(
    pilot_data = pilot_strat, time_var = "time", status_var = "status", arm_var = "arm",
    strata_var = "region", sample_sizes = c(8), linear_terms = "age",
    L = 5, n_sim = 2, alpha = 0.1, parallel.cores = 1
  )
  expect_named(add_power, c("results_data", "results_plot", "results_summary"))

  add_ss <- additive.ss.boot.app(
    pilot_data = pilot_strat, time_var = "time", status_var = "status", arm_var = "arm",
    strata_var = "region", target_power = 0, linear_terms = "age",
    L = 5, n_sim = 2, alpha = 0.1, parallel.cores = 1,
    patience = 1, n_start = 8, n_step = 5, max_n_per_arm = 8
  )
  expect_true(is.data.frame(add_ss$results_data))

  lin_power <- linear.power.boot.app(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    sample_sizes = c(8), linear_terms = "age", L = 5, n_sim = 2, alpha = 0.1, parallel.cores = 1
  )
  expect_true(is.data.frame(lin_power$results_data))

  lin_ss <- linear.ss.boot.app(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm",
    target_power = 0, linear_terms = "age", L = 5, n_sim = 2, alpha = 0.1,
    patience = 1, n_start = 8, n_step = 5, max_n_per_arm = 8, parallel.cores = 1
  )
  expect_true(is.data.frame(lin_ss$results_data))

  ms_power <- MS.power.boot.app(
    pilot_data = pilot_strat, time_var = "time", status_var = "status", arm_var = "arm",
    strata_var = "region", sample_sizes = c(8), linear_terms = "age",
    L = 5, n_sim = 2, alpha = 0.1, parallel.cores = 1
  )
  expect_true(is.data.frame(ms_power$results_data))

  ms_ss <- MS.ss.boot.app(
    pilot_data = pilot_strat, time_var = "time", status_var = "status", arm_var = "arm",
    strata_var = "region", target_power = 0, linear_terms = "age",
    L = 5, n_sim = 2, alpha = 0.1, parallel.cores = 1,
    patience = 1, n_start = 8, n_step = 5, max_n_per_arm = 8
  )
  expect_true(is.data.frame(ms_ss$results_data))
})

test_that("survival diagnostics cover log-rank and plotting branches", {
  set.seed(303)
  pilot <- make_pilot_data(n = 40, stratified = TRUE)

  res_basic <- .run_survival_diagnostics(
    pilot_data = pilot, time_var = "time", status_var = "status", arm_var = "arm", alpha = 0.1
  )
  expect_s3_class(res_basic$km_plot, "ggplot")
  expect_true(any(res_basic$logrank_summary$Statistic == "P-Value"))

  pilot_one_arm <- pilot
  pilot_one_arm$arm <- 0
  res_strata <- .run_survival_diagnostics(
    pilot_data = pilot_one_arm, time_var = "time", status_var = "status",
    arm_var = "arm", strata_var = "region", alpha = 0.1
  )
  expect_true(any(grepl("Not Applicable", res_strata$logrank_summary$Value)))
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

  res <- run_app()
  expect_true(called)
  expect_true(dir.exists(res))
})

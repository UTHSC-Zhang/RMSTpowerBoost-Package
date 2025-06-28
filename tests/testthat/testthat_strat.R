# --- Unit Tests for Multiplicative Stratified RMST Functions ---
#
# This script uses the 'testthat' package to verify the functionality
# of the power and sample size calculation functions for the stratified
# multiplicative RMST model.

library(testthat)
source(here::here("R/multiplicative_stratified.R"))
testthat::test_that("design_rmst_strat_power: input validation works", {
   # --- Test Setup ---
   set.seed(123)
   pilot_data_strat_test <- data.frame(
      time = stats::rexp(120, rate = 0.15),
      status = stats::rbinom(120, 1, 0.6),
      arm = rep(0:1, each = 60),
      region = factor(rep(c("A", "B", "C"), each = 40)),
      age = stats::rnorm(120, 55, 8)
   )
   test_time_var <- "time"; test_status_var <- "status"; test_arm_var <- "arm"
   test_strata_var <- "region"; test_tau <- 10; test_n_sim <- 10

   # --- Test Execution ---
   # This test checks for the specific error message from inside the function,
   # which is possible because `sample_sizes` defaults to NULL.
   testthat::expect_error(
      design_rmst_strat_power(
         pilot_data = pilot_data_strat_test,
         time_var = test_time_var, status_var = test_status_var, arm_var = test_arm_var, strata_var = test_strata_var,
         tau = test_tau, n_sim = test_n_sim
      ),
      "You must provide the 'sample_sizes' argument."
   )
})


testthat::test_that("design_rmst_strat_ss: input validation works", {
   # --- Test Setup ---
   set.seed(123)
   pilot_data_strat_test <- data.frame(
      time = stats::rexp(120, rate = 0.15),
      status = stats::rbinom(120, 1, 0.6),
      arm = rep(0:1, each = 60),
      region = factor(rep(c("A", "B", "C"), each = 40)),
      age = stats::rnorm(120, 55, 8)
   )
   test_time_var <- "time"; test_status_var <- "status"; test_arm_var <- "arm"
   test_strata_var <- "region"; test_tau <- 10; test_n_sim <- 10

   # --- Test Execution ---
   # Check for error when 'target_power' is missing (defaults to NULL).
   testthat::expect_error(
      design_rmst_strat_ss(
         pilot_data = pilot_data_strat_test,
         time_var = test_time_var, status_var = test_status_var, arm_var = test_arm_var, strata_var = test_strata_var,
         tau = test_tau, n_sim = test_n_sim
      ),
      "You must provide a single numeric value for 'target_power'."
   )

   # Check for the same error when the input type is wrong.
   testthat::expect_error(
      design_rmst_strat_ss(
         pilot_data = pilot_data_strat_test,
         time_var = test_time_var, status_var = test_status_var, arm_var = test_arm_var, strata_var = test_strata_var,
         target_power = c(0.8, 0.9),
         tau = test_tau, n_sim = test_n_sim
      ),
      "You must provide a single numeric value for 'target_power'."
   )
})


testthat::test_that("design_rmst_strat_power: output structure is correct", {
   # --- Test Setup ---
   set.seed(123)
   pilot_data_strat_test <- data.frame(
      time = stats::rexp(120, rate = 0.15),
      status = stats::rbinom(120, 1, 0.6),
      arm = rep(0:1, each = 60),
      region = factor(rep(c("A", "B", "C"), each = 40)),
      age = stats::rnorm(120, 55, 8)
   )
   test_time_var <- "time"; test_status_var <- "status"; test_arm_var <- "arm"
   test_strata_var <- "region"; test_tau <- 10; test_n_sim <- 10

   # --- Test Execution ---
   results <- design_rmst_strat_power(
      pilot_data = pilot_data_strat_test,
      time_var = test_time_var, status_var = test_status_var, arm_var = test_arm_var, strata_var = test_strata_var,
      sample_sizes = c(50),
      tau = test_tau, n_sim = test_n_sim
   )

   testthat::expect_type(results, "list")
   testthat::expect_named(results, c("results_data", "results_plot", "results_summary"))
   testthat::expect_s3_class(results$results_data, "data.frame")
   testthat::expect_s3_class(results$results_plot, "ggplot")
   testthat::expect_true(is.data.frame(results$results_summary) || is.null(results$results_summary))
   testthat::expect_named(results$results_data, c("N_per_Stratum", "Power"))
})


testthat::test_that("design_rmst_strat_power: output structure is correct", {
   # --- Test Setup ---
   set.seed(123)
   pilot_data_strat_test <- data.frame(
      time = stats::rexp(120, rate = 0.15),
      status = stats::rbinom(120, 1, 0.6),
      arm = rep(0:1, each = 60),
      region = factor(rep(c("A", "B", "C"), each = 40)),
      age = stats::rnorm(120, 55, 8)
   )
   test_time_var <- "time"; test_status_var <- "status"; test_arm_var <- "arm"
   test_strata_var <- "region"; test_tau <- 10; test_n_sim <- 10

   # --- Test Execution ---
   results <- design_rmst_strat_power(
      pilot_data = pilot_data_strat_test,
      time_var = test_time_var, status_var = test_status_var, arm_var = test_arm_var, strata_var = test_strata_var,
      sample_sizes = c(50),
      tau = test_tau, n_sim = test_n_sim
   )

   testthat::expect_type(results, "list")
   testthat::expect_named(results, c("results_data", "results_plot", "results_summary"))
   testthat::expect_s3_class(results$results_data, "data.frame")
   testthat::expect_s3_class(results$results_plot, "ggplot")
   testthat::expect_true(is.data.frame(results$results_summary) || is.null(results$results_summary))
   testthat::expect_named(results$results_data, c("N_per_Stratum", "Power"))
})

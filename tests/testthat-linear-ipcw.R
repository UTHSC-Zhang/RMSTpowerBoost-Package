# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html
# This script contains unit tests for the two independent IPCW functions.
# It assumes that the functions `design_rmst_ipcw_power` and `design_rmst_ipcw_ss`
# have been loaded into the environment.

library(testthat)
source(here::here("R/linear_ipcw.R"))
# --- 1. Tests for design_rmst_ipcw_power ---

test_that("Power function runs and returns correct structure", {

   # Setup: Create a reliable synthetic pilot dataset for testing
   set.seed(42)
   pilot_data_power <- base::data.frame(
      time = stats::rexp(100, rate = 0.05),
      status = stats::rbinom(100, 1, 0.8),
      arm = rep(0:1, each = 50),
      age = stats::rnorm(100, 60, 10)
   )

   # Action: Run the function with a standard call
   test_results <- design_rmst_ipcw_power(
      pilot_data = pilot_data_power,
      time_var = "time",
      status_var = "status",
      arm_var = "arm",
      sample_sizes = c(80, 100),
      linear_terms = "age",
      tau = 20,
      n_sim = 4 # Use very low n_sim for fast testing
   )

   # Assertions: Check the output's integrity
   testthat::expect_type(test_results, "list")
   testthat::expect_named(test_results, c("results_data", "results_plot", "results_summary"))

   # Check the results data frame
   testthat::expect_s3_class(test_results$results_data, "data.frame")
   testthat::expect_equal(nrow(test_results$results_data), 2)
   testthat::expect_named(test_results$results_data, c("N_per_Arm", "Power"))
   testthat::expect_true(all(test_results$results_data$Power >= 0 | is.na(test_results$results_data$Power)))
})

test_that("Power function handles automatic covariate detection", {

   set.seed(42)
   pilot_data_power <- base::data.frame(
      time = stats::rexp(100, rate = 0.05),
      status = stats::rbinom(100, 1, 0.8),
      arm = rep(0:1, each = 50),
      age = stats::rnorm(100, 60, 10),
      score = stats::rnorm(100, 10, 2)
   )

   # Action: We expect a message confirming auto-detection
   testthat::expect_message(
      results <- design_rmst_ipcw_power(
         pilot_data = pilot_data_power,
         time_var = "time",
         status_var = "status",
         arm_var = "arm",
         sample_sizes = 80,
         linear_terms = NULL, # Let the function find covariates
         tau = 20,
         n_sim = 4
      ),
      "Using unspecified columns as linear terms"
   )

   # Assertions: Check that output is still valid
   testthat::expect_s3_class(results$results_data, "data.frame")
})


# --- 2. Tests for design_rmst_ipcw_ss ---

test_that("Sample size function finds a plausible N when an effect exists", {

   # Setup: Create data WITH a treatment effect
   # (arm 1 has longer survival times: rate = 0.025 vs 0.05)
   set.seed(123)
   pilot_data_with_effect <- data.frame(
      time = c(stats::rexp(50, rate = 0.05), stats::rexp(50, rate = 0.025)),
      status = stats::rbinom(100, 1, 0.9),
      arm = rep(0:1, each = 50)
   )

   # Action: This search should now succeed without warnings
   testthat::expect_no_warning({
      results <- design_rmst_ipcw_ss(
         pilot_data = pilot_data_with_effect,
         time_var = "time",
         status_var = "status",
         arm_var = "arm",
         target_power = 0.50, # A reasonable target
         tau = 20,
         n_sim = 10,
         n_start = 50,
         n_step = 50,
         patience = 3
      )
   })

   # Assertions
   testthat::expect_s3_class(results$results_data, "data.frame")
   testthat::expect_false(is.na(results$results_data$Required_N_per_Arm))
   testthat::expect_true(results$results_data$Required_N_per_Arm > 0)
})

test_that("Sample size function handles stagnation correctly", {

   set.seed(99)
   pilot_data_ss <- base::data.frame(
      time = stats::rexp(100, rate = 0.01), # Low event rate makes power hard to find
      status = stats::rbinom(100, 1, 0.2),
      arm = rep(0:1, each = 50)
   )

   # Action: Run with an unachievable target power to trigger the patience rule.
   testthat::expect_warning(
      results <- design_rmst_ipcw_ss(
         pilot_data = pilot_data_ss,
         time_var = "time",
         status_var = "status",
         arm_var = "arm",
         target_power = 0.99, # Unachievable power
         tau = 20,
         n_sim = 4,
         n_start = 50,
         n_step = 25,
         patience = 2
      ),
      "Search terminated due to stagnation"
   )

   # Assertions: Check that it still returns the best N it found
   testthat::expect_false(is.na(results$results_data$Required_N_per_Arm))
})

test_that("Sample size function validates its parameters", {

   set.seed(99)
   pilot_data_ss <- base::data.frame(time=1, status=1, arm=1)

   # Assertion: Check that it stops if target_power is a vector
   testthat::expect_error(
      design_rmst_ipcw_ss(
         pilot_data = pilot_data_ss,
         time_var = "time",
         status_var = "status",
         arm_var = "arm",
         target_power = c(0.8, 0.9), # Invalid input
         tau = 1
      ),
      "You must provide a single numeric value for 'target_power'."
   )
})

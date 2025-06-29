
# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html
# This script contains unit tests for the two independent IPCW functions.




# --- 1. Preamble ---
library(testthat)
source(here::here("R/dependent_censoring.R"))


# --- Test Suite for Dependent Censoring Analytic Functions ---

test_that("Power function runs and returns correct structure", {
   # 1. SETUP: Simulate a pilot dataset for this test
   set.seed(42)
   n_pilot <- 500
   pilot_data <- data.frame(
      age = stats::rnorm(n_pilot, 60, 8),
      arm = rep(0:1, each = n_pilot / 2)
   )
   hazard_rates <- ifelse(pilot_data$arm == 0, 0.1, 0.07)
   pilot_data$time <- stats::rexp(n_pilot, rate = hazard_rates)
   event_probs <- stats::runif(n_pilot)
   pilot_data$primary_event <- ifelse(event_probs < 0.5, 1, 0)
   pilot_data$competing_event <- ifelse(event_probs >= 0.5 & event_probs < 0.75, 1, 0)
   pilot_data$primary_event[pilot_data$competing_event == 1] <- 0

   # 2. ACTION: Call the function with the correct name
   results <- design_rmst_dc_power(
      pilot_data = pilot_data,
      time_var = "time",
      status_var = "primary_event",
      dep_cens_status_var = "competing_event",
      arm_var = "arm",
      linear_terms = c("age"),
      sample_sizes = c(500, 1000),
      tau = 15
   )

   # 3. EXPECTATIONS: Make assertions about the output
   expect_named(results, c("results_data", "results_plot"))
   expect_s3_class(results$results_data, "data.frame")
   expect_s3_class(results$results_plot, "ggplot")
   expect_equal(nrow(results$results_data), 2)
   expect_true(all(results$results_data$Power >= 0 & results$results_data$Power <= 1))
})


test_that("Sample size function runs and returns correct structure", {
   # 1. SETUP
   set.seed(123)
   n_pilot <- 500
   pilot_data <- data.frame(
      age = stats::rnorm(n_pilot, 60, 8),
      arm = rep(0:1, each = n_pilot / 2)
   )
   hazard_rates <- ifelse(pilot_data$arm == 0, 0.1, 0.07)
   pilot_data$time <- stats::rexp(n_pilot, rate = hazard_rates)
   event_probs <- stats::runif(n_pilot)
   pilot_data$primary_event <- ifelse(event_probs < 0.5, 1, 0)
   pilot_data$competing_event <- ifelse(event_probs >= 0.5 & event_probs < 0.75, 1, 0)
   pilot_data$primary_event[pilot_data$competing_event == 1] <- 0

   # 2. ACTION: Call the function with the correct name
   results <- design_rmst_dc_ss(
      pilot_data = pilot_data,
      time_var = "time",
      status_var = "primary_event",
      dep_cens_status_var = "competing_event",
      arm_var = "arm",
      linear_terms = c("age"),
      target_power = 0.80,
      tau = 15
   )

   # 3. EXPECTATIONS
   expect_named(results, c("results_data", "results_plot", "results_summary"))
   expect_s3_class(results$results_data, "data.frame")
   expect_gt(results$results_data$Required_N_per_Arm, 0)
})


test_that("Weaker effect correctly requires a larger sample size", {
   # 1. SETUP
   set.seed(101)
   n_pilot <- 500
   pilot_data_baseline <- data.frame(arm = rep(0:1, each = n_pilot/2))
   hazard_baseline <- ifelse(pilot_data_baseline$arm == 0, 0.1, 0.06) # Strong effect
   pilot_data_baseline$time <- stats::rexp(n_pilot, rate = hazard_baseline)
   pilot_data_baseline$primary_event <- rbinom(n_pilot, 1, 0.6)
   pilot_data_baseline$competing_event <- rbinom(n_pilot, 1, 0.2)
   pilot_data_baseline$primary_event[pilot_data_baseline$competing_event==1] <- 0

   pilot_data_weaker <- pilot_data_baseline
   hazard_weaker <- ifelse(pilot_data_weaker$arm == 0, 0.1, 0.08) # Weaker effect
   pilot_data_weaker$time <- stats::rexp(n_pilot, rate = hazard_weaker)

   # 2. ACTION: Call the function with the correct name
   ss_baseline <- design_rmst_dc_ss(
      pilot_data = pilot_data_baseline, time_var = "time", status_var = "primary_event",
      dep_cens_status_var = "competing_event", arm_var = "arm",
      target_power = 0.80, tau = 20, linear_terms = NULL
   )

   ss_weaker <- design_rmst_dc_ss(
      pilot_data = pilot_data_weaker, time_var = "time", status_var = "primary_event",
      dep_cens_status_var = "competing_event", arm_var = "arm",
      target_power = 0.80, tau = 20, linear_terms = NULL
   )

   # 3. EXPECTATION
   expect_gt(ss_weaker$results_data$Required_N_per_Arm, ss_baseline$results_data$Required_N_per_Arm)
})


test_that("Functions work correctly with no covariates (linear_terms = NULL)", {
   # 1. SETUP
   set.seed(202)
   n_pilot <- 500
   pilot_data <- data.frame(arm = rep(0:1, each = n_pilot / 2))
   hazard_rates <- ifelse(pilot_data$arm == 0, 0.1, 0.07)
   pilot_data$time <- stats::rexp(n_pilot, rate = hazard_rates)
   event_probs <- stats::runif(n_pilot)
   pilot_data$primary_event <- ifelse(event_probs < 0.5, 1, 0)
   pilot_data$competing_event <- ifelse(event_probs >= 0.5 & event_probs < 0.75, 1, 0)
   pilot_data$primary_event[pilot_data$competing_event == 1] <- 0

   # 2. ACTION & EXPECTATIONS
   # Expect the functions to run without error when no covariates are provided
   expect_no_error(
      design_rmst_dc_power(
         pilot_data = pilot_data, time_var = "time", status_var = "primary_event",
         dep_cens_status_var = "competing_event", arm_var = "arm",
         linear_terms = NULL, sample_sizes = c(500, 1000), tau = 15
      )
   )

   expect_no_error(
      design_rmst_dc_ss(
         pilot_data = pilot_data, time_var = "time", status_var = "primary_event",
         dep_cens_status_var = "competing_event", arm_var = "arm",
         linear_terms = NULL, target_power = 0.80, tau = 15
      )
   )
})



















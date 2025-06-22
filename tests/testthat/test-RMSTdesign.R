#
# ==============================================================================
#
# Unit Tests for the Rmst_design Function
#
# This script uses the 'testthat' package to validate the functionality of
# the Rmst_design function across its various methods and parameters.
# It assumes that the Rmst_design function and required packages are loaded.
#
# To run:
# 1. Make sure 'testthat', 'survival', 'ggplot2', and 'mgcv' are installed.
# 2. Source the file containing the Rmst_design function.
# 3. Source this file, or run it using testthat::test_file().
#
# ==============================================================================

# --- Test Suite Setup ---

# Create a consistent pilot dataset for all tests
# This code is run once to prepare the data needed for testing.
source(here::here("R/RMSTdesign.R"))
data(pbc, package = "survival")
pilot_data <- pbc[1:312, ]
# Clean the data and create necessary variables for all method types
pilot_data <- pilot_data[stats::complete.cases(pilot_data$bili, pilot_data$age, pilot_data$stage), ]
pilot_data$stratum <- as.factor(pilot_data$stage)
# Create a status variable for the main event (death = status 2)
pilot_data$event_status <- ifelse(pilot_data$status == 2, 1, 0)
# Create a status variable for a competing event (transplant = status 1)
pilot_data$competing_event_status <- ifelse(pilot_data$status == 1, 1, 0)


# Define common parameters for tests to keep them fast
# NOTE: n_sim is kept very low. This is NOT for statistical accuracy,
#       but to ensure the code runs through all paths without error.
common_args <- list(
   pilot_data = pilot_data,
   time_var = "time",
   status_var = "status",
   arm_var = "transplant_status",
   event_indicator = 2,
   treatment_indicator = "d-penicil",
   tau = 3652.5,
   n_sim = 2, # VERY LOW for speed
   n_start = 50,
   n_step = 25,
   max_n_per_arm = 75,
   alpha = 0.05
)


# ==============================================================================
# Define Test Cases using testthat::test_that
# ==============================================================================

testthat::test_that("Input validation works correctly", {

   # Error if neither sample_sizes nor target_powers is provided
   testthat::expect_error(
      Rmst_design(pilot_data = common_args$pilot_data, time_var = "time", status_var = "status",
                  arm_var = "transplant_status", event_indicator = 2, treatment_indicator = "d-penicil", tau = 3652.5),
      "Must provide 'sample_sizes' or 'target_powers'."
   )

   # Error if both are provided
   testthat::expect_error(
      Rmst_design(sample_sizes = c(100), target_powers = c(0.8),
                  pilot_data = common_args$pilot_data, time_var = "time", status_var = "status",
                  arm_var = "transplant_status", event_indicator = 2, treatment_indicator = "d-penicil", tau = 3652.5),
      "Cannot provide both 'sample_sizes' and 'target_powers'."
   )

   # Error on invalid method name
   testthat::expect_error(
      Rmst_design(method = "invalid_method_name", sample_sizes = c(100),
                  pilot_data = common_args$pilot_data, time_var = "time", status_var = "status",
                  arm_var = "transplant_status", event_indicator = 2, treatment_indicator = "d-penicil", tau = 3652.5),
      "Invalid 'method' specified."
   )
})


testthat::test_that("Stratification validation works", {

   # Error if a non-stratified method gets a strata_var
   testthat::expect_error(
      Rmst_design(method = "linear_ipcw", sample_sizes = c(100), strata_var = "stratum",
                  pilot_data = common_args$pilot_data, time_var = "time", status_var = "status",
                  arm_var = "transplant_status", event_indicator = 2, treatment_indicator = "d-penicil", tau = 3652.5),
      "does not support stratification"
   )

   # Error if a stratified method does NOT get a strata_var
   testthat::expect_error(
      Rmst_design(method = "multiplicative_stratified", sample_sizes = c(100), strata_var = NULL,
                  pilot_data = common_args$pilot_data, time_var = "time", status_var = "status",
                  arm_var = "transplant_status", event_indicator = 2, treatment_indicator = "d-penicil", tau = 3652.5),
      "requires a 'strata_var'"
   )
})


testthat::test_that("Method 'additive_gam' runs in all modes", {

   # Non-stratified, power mode
   res_ns_p <- Rmst_design(
      method = "additive_gam", sample_sizes = c(100), strata_var = NULL,
      linear_terms = "age", smooth_terms = "bili",
      pilot_data = common_args$pilot_data, time_var = "time", status_var = "status",
      arm_var = "transplant_status", event_indicator = 2, treatment_indicator = "d-penicil",
      tau = common_args$tau, n_sim = common_args$n_sim
   )
   testthat::expect_true(is.list(res_ns_p))
   testthat::expect_equal(names(res_ns_p), c("results_data", "results_plot"))
   testthat::expect_s3_class(res_ns_p$results_data, "data.frame")
   testthat::expect_s3_class(res_ns_p$results_plot, "ggplot")
   testthat::expect_equal(names(res_ns_p$results_data), c("N_per_Arm", "Power"))

   # Stratified, sample size mode
   res_s_ss <- Rmst_design(
      method = "additive_gam", target_powers = c(0.8), strata_var = "stratum",
      linear_terms = "age", smooth_terms = "bili",
      pilot_data = common_args$pilot_data, time_var = "time", status_var = "status",
      arm_var = "transplant_status", event_indicator = 2, treatment_indicator = "d-penicil",
      tau = common_args$tau, n_sim = common_args$n_sim, n_start = 50, n_step = 25, max_n_per_arm = 75
   )
   testthat::expect_equal(names(res_s_ss$results_data), c("Target_Power", "Required_N_per_Arm"))
})


testthat::test_that("Method 'linear_ipcw' runs correctly", {

   res_p <- Rmst_design(
      method = "linear_ipcw", sample_sizes = c(100), strata_var = NULL,
      linear_terms = c("age", "bili"),
      pilot_data = common_args$pilot_data, time_var = "time", status_var = "status",
      arm_var = "transplant_status", event_indicator = 2, treatment_indicator = "d-penicil",
      tau = common_args$tau, n_sim = common_args$n_sim
   )
   testthat::expect_true(is.list(res_p))
   testthat::expect_equal(names(res_p), c("results_data", "results_plot"))
})


testthat::test_that("Method 'multiplicative_stratified' runs correctly", {

   res_ss <- Rmst_design(
      method = "multiplicative_stratified", target_powers = c(0.8), strata_var = "stratum",
      linear_terms = c("age", "bili"),
      pilot_data = common_args$pilot_data, time_var = "time", status_var = "status",
      arm_var = "transplant_status", event_indicator = 2, treatment_indicator = "d-penicil",
      tau = common_args$tau, n_sim = common_args$n_sim, n_start = 50, n_step = 25, max_n_per_arm = 75
   )
   testthat::expect_true(is.list(res_ss))
   testthat::expect_equal(names(res_ss$results_data), c("Target_Power", "Required_N_per_Arm"))
})


testthat::test_that("Method 'linear_dependent_censoring' runs correctly", {

   # Error if dependent censoring variable is missing
   testthat::expect_error(
      Rmst_design(method = "linear_dependent_censoring", sample_sizes = c(100),
                  pilot_data = common_args$pilot_data, time_var = "time", status_var = "status",
                  arm_var = "transplant_status", event_indicator = 2, treatment_indicator = "d-penicil", tau = 3652.5),
      "requires 'dep_cens_status_var'"
   )

   # Runs correctly with all parameters
   res_p <- Rmst_design(
      method = "linear_dependent_censoring", sample_sizes = c(150),
      linear_terms = c("age", "bili"),
      pilot_data = pilot_data, # Using the test-specific data
      time_var = "time",
      status_var = "status",
      arm_var = "transplant_status",
      dep_cens_status_var = "status",
      event_indicator = 2,
      dep_cens_event_indicator = 1, # status=1 is competing event
      treatment_indicator = "d-penicil",
      tau = common_args$tau,
      n_sim = common_args$n_sim
   )
   testthat::expect_true(is.list(res_p))
   testthat::expect_equal(names(res_p$results_data), c("N_per_Arm", "Power"))
})

cat("\nAll tests passed successfully!\n")

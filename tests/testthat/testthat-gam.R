
# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html
# This script contains unit tests for the two independent IPCW functions.



# --- Data for non-stratified tests ---
pilot_data_gam_simple <- base::data.frame(
   time = stats::rexp(100, rate = 0.08),
   status = stats::rbinom(100, 1, 0.7),
   arm = rep(0:1, each = 50),
   age = stats::rnorm(100, 60, 10),
   biomarker = stats::rnorm(100, 15, 4)
)

# --- Data for stratified tests ---
pilot_data_gam_strat <- base::data.frame(
   time = stats::rexp(120, rate = 0.08),
   status = stats::rbinom(120, 1, 0.7),
   arm = rep(0:1, each = 60),
   age = stats::rnorm(120, 60, 10),
   biomarker = stats::rnorm(120, 15, 4),
   region = base::factor(rep(c("NA", "EU", "AS"), each = 40))
)


# --- 1. Tests for design_rmst_gam_power ---

testthat::test_that("GAM Power function runs (non-stratified)", {
   # Action & Assertion: Check for successful completion
   testthat::expect_no_error({
      results <- design_rmst_gam_power(
         pilot_data = pilot_data_gam_simple,
         time_var = "time",
         status_var = "status",
         arm_var = "arm",
         sample_sizes = c(100),
         linear_terms = "age",
         smooth_terms = "biomarker",
         tau = 15,
         n_sim = 4, # Low n_sim for fast testing
         parallel.cores = 1
      )
   })

   # Assertions: Check output structure
   testthat::expect_type(results, "list")
   testthat::expect_named(results, c("results_data", "results_plot", "results_summary"))
   testthat::expect_s3_class(results$results_data, "data.frame")
   testthat::expect_named(results$results_data, c("N_per_Group", "Power"))
})

testthat::test_that("GAM Power function runs (stratified & parallel)", {
   # Action & Assertion: Check stratified case with parallel processing
   testthat::expect_no_error({
      results <- design_rmst_gam_power(
         pilot_data = pilot_data_gam_strat,
         time_var = "time",
         status_var = "status",
         arm_var = "arm",
         strata_var = "region",
         sample_sizes = c(50), # N per stratum
         linear_terms = "age",
         smooth_terms = "biomarker",
         tau = 15,
         n_sim = 4,
         parallel.cores = 2 # Test parallel execution
      )
   })

   # Assertions
   testthat::expect_type(results, "list")
   testthat::expect_s3_class(results$results_data, "data.frame")
})

testthat::test_that("GAM Power function issues message for auto-detection", {
   testthat::expect_message(
      design_rmst_gam_power(
         pilot_data = pilot_data_gam_simple,
         time_var = "time",
         status_var = "status",
         arm_var = "arm",
         sample_sizes = 100,
         linear_terms = NULL, # Trigger auto-detection
         smooth_terms = NULL,
         tau = 15,
         n_sim = 2
      ),
      "Using unspecified columns as linear terms"
   )
})

testthat::test_that("GAM Power function errors if sample_sizes is missing", {
   # CORRECTED: The test now expects the actual error message from R.
   testthat::expect_error(
      design_rmst_gam_power(
         pilot_data = pilot_data_gam_simple,
         time_var = "time",
         status_var = "status",
         arm_var = "arm",
         tau = 15,
         n_sim = 2
      ),
      "argument \"sample_sizes\" is missing, with no default"
   )
})


# --- 2. Tests for design_rmst_gam_ss ---

testthat::test_that("GAM Sample Size function finds a plausible N when an effect exists", {

   # Setup: Create data WITH a treatment effect for a stable test
   set.seed(321)
   pilot_data_ss_effect <- data.frame(
      time = c(stats::rexp(50, rate = 0.1), stats::rexp(50, rate = 0.04)),
      status = stats::rbinom(100, 1, 0.9),
      arm = rep(0:1, each = 50),
      age = stats::rnorm(100, 60, 10)
   )

   # Action & Assertion: Should succeed without warnings
   testthat::expect_no_warning({
      results <- design_rmst_gam_ss(
         pilot_data = pilot_data_ss_effect,
         time_var = "time",
         status_var = "status",
         arm_var = "arm",
         target_power = 0.60, # A reasonable target
         linear_terms = "age",
         tau = 15,
         n_sim = 10,
         n_start = 50,
         n_step = 50,
         patience = 2,
         parallel.cores = 1
      )
   })

   # Assertions
   testthat::expect_s3_class(results$results_data, "data.frame")
   testthat::expect_false(is.na(results$results_data$Required_N_per_Group))
})

testthat::test_that("GAM Sample Size function handles stagnation correctly", {

   # Setup: Use data with no effect to ensure stagnation
   pilot_data_ss_no_effect <- data.frame(
      time = stats::rexp(100, rate = 0.05),
      status = stats::rbinom(100, 1, 0.5),
      arm = rep(0:1, each = 50)
   )

   # Action & Assertion: Expect a stagnation warning
   testthat::expect_warning(
      results <- design_rmst_gam_ss(
         pilot_data = pilot_data_ss_no_effect,
         time_var = "time",
         status_var = "status",
         arm_var = "arm",
         target_power = 0.95, # Unachievable power
         tau = 20,
         n_sim = 4,
         n_start = 50,
         n_step = 25,
         patience = 2
      ),
      "Search terminated due to stagnation"
   )

   # Assertions: Still returns the best N found
   testthat::expect_false(is.na(results$results_data$Required_N_per_Group))
})

testthat::test_that("GAM Sample Size function validates target_power input", {
   testthat::expect_error(
      design_rmst_gam_ss(
         pilot_data = pilot_data_gam_simple,
         time_var = "time",
         status_var = "status",
         arm_var = "arm",
         target_power = c(0.8, 0.9), # Should be a single value
         tau = 15,
         n_sim = 2
      ),
      "You must provide a single numeric value for 'target_power'"
   )
})

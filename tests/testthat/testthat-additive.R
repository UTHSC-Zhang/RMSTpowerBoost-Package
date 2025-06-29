# Test script for Additive RMST model functions (Boot & Analytic)
library(testthat)
devtools::load_all(".")


# --- Test Data Setup ---
# Unstratified data with an additive effect
pilot_data_add_simple <- data.frame(
   time = stats::rexp(100, rate = 0.08),
   status = stats::rbinom(100, 1, 0.8),
   arm = rep(0:1, each = 50),
   age = stats::rnorm(100, 60, 10)
)
pilot_data_add_simple$time[pilot_data_add_simple$arm == 1] <-
   pilot_data_add_simple$time[pilot_data_add_simple$arm == 1] + 2.5 # Additive effect

# Stratified data with an additive effect
pilot_data_add_strat <- data.frame(
   time = stats::rexp(120, rate = 0.08),
   status = stats::rbinom(120, 1, 0.8),
   arm = rep(0:1, each = 60),
   age = stats::rnorm(120, 60, 10),
   region = factor(rep(c("NA", "EU", "AS"), each = 40))
)
pilot_data_add_strat$time[pilot_data_add_strat$arm == 1] <-
   pilot_data_add_strat$time[pilot_data_add_strat$arm == 1] + 2.5 # Additive effect


# --- Test Suite for Additive Bootstrap Functions ---

test_that("GAM.power.boot works for non-stratified case", {
   set.seed(123)
   results <- GAM.power.boot(
      pilot_data = pilot_data_add_simple,
      time_var = "time", status_var = "status", arm_var = "arm",
      sample_sizes = c(100), tau = 15, n_sim = 10
   )
   # CORRECTED: Use expect_type for a base list
   expect_type(results, "list")
   expect_named(results, c("results_data", "results_plot", "results_summary"))
   expect_s3_class(results$results_data, "data.frame")
   expect_named(results$results_data, c("N_per_Group", "Power"))
})

test_that("GAM.power.boot works for stratified case", {
   set.seed(456)
   results <- GAM.power.boot(
      pilot_data = pilot_data_add_strat,
      time_var = "time", status_var = "status", arm_var = "arm", strata_var = "region",
      sample_sizes = c(50), tau = 15, n_sim = 10
   )
   # CORRECTED: Use expect_type for a base list
   expect_type(results, "list")
   expect_s3_class(results$results_data, "data.frame")
})

test_that("GAM.ss.boot finds a plausible N", {
   set.seed(789)
   results <- GAM.ss.boot(
      pilot_data = pilot_data_add_simple,
      time_var = "time", status_var = "status", arm_var = "arm",
      target_power = 0.6, tau = 15, n_sim = 10,
      n_start = 50, n_step = 50, patience = 2
   )
   expect_s3_class(results$results_data, "data.frame")
   expect_true(results$results_data$Required_N_per_Group > 0)
})


# --- Test Suite for Additive Analytical Functions ---

test_that("additive.power.analytical works for stratified case", {
   set.seed(111)
   results <- additive.power.analytical(
      pilot_data = pilot_data_add_strat,
      time_var = "time", status_var = "status", arm_var = "arm", strata_var = "region",
      sample_sizes = c(100, 150), tau = 15
   )
   expect_type(results, "list")
   expect_named(results$results_data, c("N_per_Stratum", "Power"))
})

test_that("additive.ss.analytical finds a plausible N", {
   set.seed(222)
   results <- additive.ss.analytical(
      pilot_data = pilot_data_add_strat,
      time_var = "time", status_var = "status", arm_var = "arm", strata_var = "region",
      target_power = 0.8, tau = 15
   )
   expect_s3_class(results$results_data, "data.frame")
   expect_true(results$results_data$Required_N_per_Stratum > 0)
})

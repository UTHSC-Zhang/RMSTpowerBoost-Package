# Test script for Dependent Censoring functions
library(testthat)

# --- Test Data Setup ---
# A pilot dataset with a primary event and a competing event
generate_dc_data <- function(seed, n = 200, effect = 1.2) {
   set.seed(seed)
   pilot_data <- data.frame(
      age = stats::rnorm(n, 60, 8),
      arm = rep(0:1, each = n / 2)
   )
   hazard_rates <- ifelse(pilot_data$arm == 0, 0.1, 0.1 / effect)
   pilot_data$time <- stats::rexp(n, rate = hazard_rates)

   # Assign events: 60% primary, 20% competing, 20% independent censoring
   event_type <- sample(0:2, n, replace = TRUE, prob = c(0.6, 0.2, 0.2))
   pilot_data$status <- ifelse(event_type == 0, 1, 0)
   pilot_data$comp_event <- ifelse(event_type == 1, 1, 0)
   # Shorten time for non-primary events
   pilot_data$time[event_type != 0] <- pilot_data$time[event_type != 0] * 0.7

   return(pilot_data)
}


# --- Test Suite for DC.power.analytical ---

test_that("DC.power.analytical returns correct structure", {
   pilot_data <- generate_dc_data(seed = 123)

   results <- DC.power.analytical(
      pilot_data = pilot_data,
      time_var = "time",
      status_var = "status",
      dep_cens_status_var = "comp_event",
      arm_var = "arm",
      sample_sizes = c(500, 1000),
      tau = 15
   )

   # CORRECTED: Use expect_type for a base list
   expect_type(results, "list")
   expect_named(results, c("results_data", "results_plot"))
   expect_s3_class(results$results_data, "data.frame")
   expect_s3_class(results$results_plot, "ggplot")
   expect_true(all(results$results_data$Power >= 0 & results$results_data$Power <= 1))
})


# --- Test Suite for DC.ss.analytical ---

test_that("DC.ss.analytical returns correct structure", {
   pilot_data <- generate_dc_data(seed = 456)
   results <- DC.ss.analytical(
      pilot_data = pilot_data,
      time_var = "time",
      status_var = "status",
      dep_cens_status_var = "comp_event",
      arm_var = "arm",
      target_power = 0.80,
      tau = 15
   )
   # CORRECTED: Use expect_type for a base list
   expect_type(results, "list")
   expect_named(results, c("results_data", "results_plot", "results_summary"))
   expect_gt(results$results_data$Required_N_per_Arm, 0)
})

test_that("DC.ss.analytical requires larger N for weaker effect", {
   data_strong_effect <- generate_dc_data(seed = 101, effect = 1.7)
   data_weak_effect <- generate_dc_data(seed = 101, effect = 1.2)

   ss_strong <- DC.ss.analytical(
      pilot_data = data_strong_effect, time_var = "time", status_var = "status",
      dep_cens_status_var = "comp_event", arm_var = "arm",
      target_power = 0.80, tau = 20
   )

   ss_weak <- DC.ss.analytical(
      pilot_data = data_weak_effect, time_var = "time", status_var = "status",
      dep_cens_status_var = "comp_event", arm_var = "arm",
      target_power = 0.80, tau = 20
   )

   expect_gt(ss_weak$results_data$Required_N_per_Arm, ss_strong$results_data$Required_N_per_Arm)
})

test_that("DC functions work correctly with no covariates", {
   pilot_data <- generate_dc_data(seed = 202)

   expect_no_error(
      DC.power.analytical(
         pilot_data = pilot_data, time_var = "time", status_var = "status",
         dep_cens_status_var = "comp_event", arm_var = "arm",
         linear_terms = NULL, sample_sizes = c(500), tau = 15
      )
   )

   expect_no_error(
      DC.ss.analytical(
         pilot_data = pilot_data, time_var = "time", status_var = "status",
         dep_cens_status_var = "comp_event", arm_var = "arm",
         linear_terms = NULL, target_power = 0.80, tau = 15
      )
   )
})

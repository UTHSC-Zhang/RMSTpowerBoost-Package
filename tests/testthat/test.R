library(testthat)
library(RMSTSS)

set.seed(123)

test_that("estimate_rmst returns correct RMST for known data", {
   time <- c(1, 2, 3, 4, 5)
   status <- c(1, 1, 0, 1, 0)
   tau <- 4
   result <- estimate_rmst(time, status, tau)
   expect_true(is.numeric(result))
   expect_gt(result, 0)
})

test_that("simulate_stratified_survival generates correct structure", {
   sim_data <- simulate_stratified_survival(
      n_per_stratum = c(10, 20),
      hazard_rates = c(0.1, 0.2),
      censoring_rate = 10,
      strata_names = c("A", "B")
   )
   expect_s3_class(sim_data, "data.frame")
   expect_true(all(c("time", "status", "stratum") %in% names(sim_data)))
   expect_equal(length(unique(sim_data$stratum)), 2)
   expect_equal(sum(sim_data$stratum == "A"), 10)
   expect_equal(sum(sim_data$stratum == "B"), 20)
})

test_that("rmst_sample_size returns numeric vector with positive integers", {
   sim_data <- simulate_stratified_survival(
      n_per_stratum = c(10, 10),
      hazard_rates = c(0.05, 0.1),
      censoring_rate = 5,
      strata_names = c("S1", "S2")
   )

   result <- rmst_sample_size(
      data = sim_data,
      tau = 3,
      effect_size = 0.5,
      alpha = 0.05,
      power = 0.8
   )
   expect_type(result, "double")
   expect_true(all(result > 0))
   expect_named(result)
})

test_that("RMSt.power computes power correctly", {
   power <- RMSt.power(delta = 1, sigma = 1.5, n_total = 100, alpha = 0.05)
   expect_true(is.numeric(power))
   expect_gt(power, 0)
   expect_lt(power, 1)

   power_high <- RMSt.power(delta = 3, sigma = 1, n_total = 100, alpha = 0.05)
   expect_gt(power_high, power)
})

test_that("RMSt.sample.size returns valid sample size", {
   n <- RMSt.sample.size(delta = 0.5, sigma2 = 0.25, alpha = 0.05, power = 0.8)
   expect_type(n, "double")
   expect_gt(n, 0)
   expect_equal(n, ceiling(n))  # should be an integer
})

test_that("pooled_variance computes correct value", {
   var_vec <- c(1, 2, 3)
   n_vec <- c(10, 20, 30)
   pooled <- pooled_variance(var_vec, n_vec)
   expect_true(is.numeric(pooled))
   expect_gt(pooled, 0)
})

test_that("estimate_rmst returns correct structure", {
   n <- 100
   time <- rexp(n)
   event <- rbinom(n, 1, 0.7)
   group <- rbinom(n, 1, 0.5)
   tau <- 5

   result <- estimate_rmst(time, event, group, tau)
   expect_true(is.list(result))
   expect_named(result, c("rmst_trt", "rmst_ctrl", "rmst_diff", "rmst_var"))
})

test_that("rmst_ss_calc handles edge cases", {
   expect_error(rmst_ss_calc(delta = 0, sigma2 = 0.2, alpha = 0.05, power = 0.8),
                "Delta must be non-zero")
   expect_error(rmst_ss_calc(delta = 1, sigma2 = -1, alpha = 0.05, power = 0.8),
                "Variance must be positive")
})

test_that("rmst_power_calc returns power between 0 and 1", {
   power <- rmst_power_calc(delta = 2, sigma2 = 1, alpha = 0.05, n1 = 50, n2 = 50)
   expect_true(is.numeric(power))
   expect_true(power > 0 && power < 1)
})

test_that("rmst_power_sim returns expected structure", {
   result <- rmst_power_sim(n1 = 50, n2 = 50, delta = 1,
                            lambda1 = 0.05, lambda0 = 0.1,
                            tau = 10, nsim = 10, alpha = 0.05)
   expect_named(result, c("estimated_power", "power_se"))
   expect_true(result$estimated_power >= 0 && result$estimated_power <= 1)
})

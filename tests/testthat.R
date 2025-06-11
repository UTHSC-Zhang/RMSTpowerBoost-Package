library(testthat)
library(RMSTSS)

set.seed(123)

# ---- estimate_rmst ----
test_that("estimate_rmst computes RMST correctly", {
   time <- c(1, 2, 3, 4, 5)
   status <- c(1, 1, 0, 1, 0)
   tau <- 4
   result <- estimate_rmst(time, status, tau)
   expect_true(is.numeric(result))
   expect_gt(result, 0)
})

test_that("estimate_rmst handles tau beyond max time", {
   time <- c(1, 2, 3)
   status <- c(1, 0, 1)
   tau <- 10
   result <- estimate_rmst(time, status, tau)
   expect_true(is.numeric(result))
})

test_that("estimate_rmst errors with mismatched lengths", {
   expect_error(
      estimate_rmst(time = c(1, 2), status = c(1, 0, 1), tau = 3),
      regexp = "length"
   )
})


# ---- simulate_stratified_survival ----
test_that("simulate_stratified_survival output is valid", {
   sim_data <- simulate_stratified_survival(
      n_per_stratum = c(5, 10),
      hazard_rates = c(0.1, 0.2),
      censoring_rate = 5
   )
   expect_s3_class(sim_data, "data.frame")
   expect_true(all(c("time", "status", "stratum") %in% names(sim_data)))
})

test_that("simulate_stratified_survival errors with mismatched inputs", {
   expect_error(
      simulate_stratified_survival(n_per_stratum = c(5, 10), hazard_rates = c(0.1), censoring_rate = 5),
      regexp = "Lengths of n_per_stratum and hazard_rates must be equal"
   )
})


# ---- rmst_sample_size ----
test_that("rmst_sample_size returns sample sizes per stratum", {
   sim_data <- simulate_stratified_survival(
      n_per_stratum = c(10, 10),
      hazard_rates = c(0.05, 0.1),
      censoring_rate = 3,
      strata_names = c("S1", "S2")
   )

   result <- rmst_sample_size(
      data = sim_data,
      tau = 2,
      effect_size = 0.5,
      alpha = 0.05,
      power = 0.8
   )
   expect_type(result, "double")
   expect_true(all(result > 0))
   expect_named(result, c("S1", "S2"))
})

test_that("rmst_sample_size errors without stratum levels", {
   bad_data <- data.frame(time = c(1, 2), status = c(1, 0), stratum = c("A", "B"))
   expect_error(
      rmst_sample_size(data = bad_data, tau = 3, effect_size = 0.5),
      regexp = "must have a factor 'stratum'"
   )
})

test_that("rmst_sample_size errors with zero effect size", {
   sim_data <- simulate_stratified_survival(
      n_per_stratum = c(10, 10),
      hazard_rates = c(0.1, 0.2),
      censoring_rate = 4
   )
   expect_error(
      rmst_sample_size(sim_data, tau = 3, effect_size = 0),
      regexp = "zero"
   )
})


# ---- RMSt.power ----
test_that("RMSt.power returns power within (0, 1)", {
   power <- RMSt.power(delta = 1, sigma = 2, n_total = 100)
   expect_true(is.numeric(power))
   expect_gt(power, 0)
   expect_lt(power, 1)
})

test_that("RMSt.power handles odd total sample size", {
   power <- suppressWarnings(RMSt.power(delta = 1, sigma = 2, n_total = 99))
   expect_true(is.numeric(power))
})

test_that("RMSt.power returns near 1 for large effect size", {
   power <- RMSt.power(delta = 5, sigma = 1, n_total = 100)
   expect_gt(power, 0.99)
})


# ---- RMSt.sample.size ----
test_that("RMSt.sample.size computes valid sample size", {
   n <- RMSt.sample.size(delta = 0.5, sigma2 = 0.25)
   expect_type(n, "double")
   expect_gt(n, 0)
   expect_equal(n, ceiling(n))
})

test_that("RMSt.sample.size errors if delta is 0", {
   expect_error(RMSt.sample.size(delta = 0, sigma2 = 0.1), regexp = "Effect size delta cannot be zero")
})


# ---- pooled_variance ----
test_that("pooled_variance computes correctly", {
   var_vec <- c(1, 2)
   n_vec <- c(10, 20)
   result <- pooled_variance(var_vec, n_vec)
   expect_true(is.numeric(result))
   expect_gt(result, 0)
})

test_that("pooled_variance errors on mismatched inputs", {
   expect_error(pooled_variance(c(1, 2), c(10)), regexp = "Length of var_vec and n_vec")
})

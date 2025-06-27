

# Helper function to check if a value is between 0 and 1
expect_between <- function(object, lower, upper) {
   act <- testthat::quasi_label(rlang::enquo(object), arg = "object")
   testthat::expect(
      object >= lower && object <= upper,
      sprintf("%s is not between %s and %s. Actual value: %s", act$lab, lower, upper, object)
   )
}



test_that("RMST calculation is correct for a simple case", {
   # Single event at time 5, tau = 10
   # Survival is 1 up to time 5, and 0 after.
   # Area should be 1 * 5 = 5.
   time <- c(5)
   status <- c(1)
   tau <- 10
   rmst <- estimate_rmst(time, status, tau)
   expect_equal(rmst, 5)
})

test_that("RMST calculation handles censoring correctly", {
   # One event at 5, one censored at 8. tau = 10
   # S(t) is 1 until t=5, then 0.5 until t=10.
   # Area = (5-0)*1 + (10-5)*0.5 = 5 + 2.5 = 7.5
   time <- c(5, 8)
   status <- c(1, 0)
   tau <- 10
   rmst <- estimate_rmst(time, status, tau)
   expect_equal(rmst, 7.5)
})

test_that("RMST calculation respects truncation time tau", {
   # Event at 5, censored at 8. tau = 6
   # S(t) is 1 until t=5, then 0.5.
   # Area = (5-0)*1 + (6-5)*0.5 = 5 + 0.5 = 5.5
   time <- c(5, 8)
   status <- c(1, 0)
   tau <- 6
   rmst <- estimate_rmst(time, status, tau)
   expect_equal(rmst, 5.5)
})

test_that("RMST is equal to tau when there are no events", {
   time <- c(10, 12, 15)
   status <- c(0, 0, 0)
   tau <- 8
   rmst <- estimate_rmst(time, status, tau)
   expect_equal(rmst, tau)
})




test_that("Simulation function produces correct data structure", {
   sim_data <- simulate_stratified_survival(
      n_per_stratum = c(10, 20),
      hazard_rates = c(0.1, 0.5),
      censoring_rate = 10
   )
   expect_s3_class(sim_data, "data.frame")
   expect_equal(nrow(sim_data), 30)
   expect_equal(colnames(sim_data), c("time", "status", "stratum"))
   expect_s3_class(sim_data$stratum, "factor")
   expect_equal(levels(sim_data$stratum), c("Stratum1", "Stratum2"))
})

test_that("Simulation function respects custom strata names", {
   sim_data <- simulate_stratified_survival(
      n_per_stratum = c(5, 5),
      hazard_rates = c(0.1, 0.5),
      censoring_rate = 10,
      strata_names = c("A", "B")
   )
   expect_equal(levels(sim_data$stratum), c("A", "B"))
})

test_that("Simulation function validates input lengths", {
   expect_error(
      simulate_stratified_survival(
         n_per_stratum = c(10, 20),
         hazard_rates = c(0.1), # Mismatched length
         censoring_rate = 10
      ),
      "Lengths of n_per_stratum and hazard_rates must be equal."
   )
})




test_that("Pooled variance calculation is correct", {
   # Simple case with two strata
   var_vec <- c(10, 20)
   n_vec <- c(100, 150)
   # Expected: ((99*10) + (149*20)) / (100+150-2)
   expected_pooled_var <- ((99 * 10) + (149 * 20)) / (250 - 2)
   expect_equal(pooled_variance(var_vec, n_vec), expected_pooled_var)
})

test_that("Pooled variance handles equal variance case", {
   expect_equal(pooled_variance(c(25, 25), c(10, 10)), 25)
})

test_that("Pooled variance validates input lengths", {
   expect_error(
      pooled_variance(c(10, 20), c(100)),
      "Length of var_vec and n_vec must be equal."
   )
})




test_that("Sample size function returns correct structure", {
   set.seed(123) # for reproducibility of bootstrap
   sim_data <- simulate_stratified_survival(
      n_per_stratum = c(50, 50),
      hazard_rates = c(0.1, 0.2),
      censoring_rate = 15,
      strata_names = c("Control", "Treatment")
   )

   ss_req <- rmst_sample_size(
      data = sim_data,
      tau = 10,
      effect_size = 0.5
   )

   expect_true(is.numeric(ss_req))
   expect_equal(length(ss_req), 2)
   expect_equal(names(ss_req), c("Control", "Treatment"))
   expect_true(all(ss_req > 0))
})

test_that("Sample size function requires a factor stratum column", {
   sim_data <- data.frame(time = 1:10, status = 1, stratum = "A")
   expect_error(
      rmst_sample_size(data = sim_data, tau = 5, effect_size = 1),
      "Data must have a factor 'stratum' column with levels."
   )
})




test_that("RMSt.power calculations are logical", {
   # Power should be > alpha if delta is non-zero
   expect_gt(RMSt.power(delta = 1.2, sigma = 2, n_total = 100, alpha = 0.05), 0.05)
   # Power approaches 1 for large effect size or sample size
   expect_equal(round(RMSt.power(delta = 10, sigma = 2, n_total = 100), 2), 1.00)
   expect_equal(round(RMSt.power(delta = 1, sigma = 2, n_total = 1000), 2), 1.00)
   # Power decreases with more variance
   power1 <- RMSt.power(delta = 1, sigma = 2, n_total = 100)
   power2 <- RMSt.power(delta = 1, sigma = 4, n_total = 100)
   expect_lt(power2, power1)
   # Power should be between 0 and 1
   expect_between(power1, 0, 1)
})

test_that("RMSt.power warns for odd sample size", {
   expect_warning(
      RMSt.power(delta = 1, sigma = 2, n_total = 101),
      "Total sample size n_total should be even. Rounding down."
   )
})

test_that("RMSt.sample.size calculations are logical", {
   # Required N increases with variance and power
   n1 <- RMSt.sample.size(delta = 0.5, sigma2 = 0.25, power = 0.8)
   n2 <- RMSt.sample.size(delta = 0.5, sigma2 = 0.50, power = 0.8) # more variance
   n3 <- RMSt.sample.size(delta = 0.5, sigma2 = 0.25, power = 0.9) # more power
   expect_gt(n2, n1)
   expect_gt(n3, n1)

   # Required N decreases with larger effect size
   n4 <- RMSt.sample.size(delta = 1.0, sigma2 = 0.25, power = 0.8)
   expect_lt(n4, n1)

   # Check a known value
   # For alpha=0.05, power=0.8, z_alpha=1.96, z_beta=0.84
   # n = (2 * sigma2 * (1.96 + 0.84)^2) / delta^2
   delta <- 0.5
   sigma2 <- 0.25
   z_alpha <- qnorm(1 - 0.05 / 2)
   z_beta <- qnorm(0.8)
   expected_n <- (2 * sigma2 * (z_alpha + z_beta)^2) / delta^2
   expect_equal(RMSt.sample.size(delta, sigma2), ceiling(expected_n))
})

test_that("RMSt.sample.size rejects delta of zero", {
   expect_error(
      RMSt.sample.size(delta = 0, sigma2 = 0.25),
      "Effect size delta cannot be zero."
   )
})



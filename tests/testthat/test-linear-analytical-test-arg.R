# The analytic engine must describe the same test the bootstrap engine runs.
analytic_pilot <- function() {
  set.seed(4417)
  d <- data.frame(arm = rep(0:1, each = 250), z = rnorm(500, 4, 1))
  event <- rexp(500, exp(-2 + .15 * d$z - .3 * d$arm))
  censor <- rexp(500, exp(-3 + .4 * d$z + .5 * d$arm))
  d$time <- pmin(event, censor)
  d$status <- as.integer(event <= censor)
  d
}

test_that("test = 'wls' uses the WLS threshold with the sandwich sampling variance", {
  d <- analytic_pilot()
  e <- RMSTpowerBoost:::.estimate_linear_params(d, "time", "status", "arm", "z", 2)
  k <- e$arm_coeff_name
  # The model-based SE is the one summary(lm) reports, on the n = 1 scale.
  expect_equal(e$se_beta_model_n1 / sqrt(e$n_pilot),
               unname(summary(e$fit_lm)$coefficients[k, 2]), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(e$se_beta_n1, e$se_beta_model_n1)))

  z <- qnorm(.975)
  for (n_arm in c(200, 400)) {
    N <- 2 * n_arm
    s_sw <- e$se_beta_n1 / sqrt(N)
    s_ml <- e$se_beta_model_n1 / sqrt(N)
    b <- abs(e$beta_effect)
    expect_equal(
      unname(RMSTpowerBoost:::.rmst_wald_power(e$beta_effect, e$se_beta_n1, N, z,
                                               e$se_beta_model_n1)),
      unname(pnorm((b - z * s_ml) / s_sw) + pnorm((-b - z * s_ml) / s_sw)),
      tolerance = 1e-12)
  }
})

test_that("the default is wls, sandwich is opt-in, and other engines are unaffected", {
  d <- analytic_pilot()
  sizes <- c(200, 400)
  default <- linear.power.analytical(d, "time", "status", "arm", sizes, "z", 2)
  wls <- linear.power.analytical(d, "time", "status", "arm", sizes, "z", 2, test = "wls")
  sand <- linear.power.analytical(d, "time", "status", "arm", sizes, "z", 2, test = "sandwich")
  expect_equal(default$results_data$Power, wls$results_data$Power)
  # The WLS test pays for an inflated standard error, so it is the less powerful one.
  expect_true(all(wls$results_data$Power < sand$results_data$Power))
  # Reported SEs follow the selected test.
  e <- RMSTpowerBoost:::.estimate_linear_params(d, "time", "status", "arm", "z", 2)
  expect_equal(wls$model_output$treatment_effect$std_error,
               unname(summary(e$fit_lm)$coefficients[e$arm_coeff_name, 2]), tolerance = 1e-10)
  expect_equal(sand$model_output$treatment_effect$std_error,
               unname(e$se_beta_n1 / sqrt(e$n_pilot)), tolerance = 1e-10)
  expect_equal(wls$model_output$variance_components$se_effect_sandwich_n1,
               sand$model_output$variance_components$se_effect_sandwich_n1)
  # A larger required N under the less powerful test.
  ss_w <- linear.ss.analytical(d, "time", "status", "arm", .8, "z", 2,
                               n_start = 100, n_step = 100, max_n_per_arm = 20000)
  ss_s <- linear.ss.analytical(d, "time", "status", "arm", .8, "z", 2,
                               n_start = 100, n_step = 100, max_n_per_arm = 20000,
                               test = "sandwich")
  expect_gte(ss_w$results_data$Required_N_per_Arm, ss_s$results_data$Required_N_per_Arm)
  # Omitting se_test_n1 must reproduce the original one-tailed formula exactly,
  # which is what the additive, multiplicative and DC engines still call.
  expect_identical(RMSTpowerBoost:::.rmst_wald_power(.5, 3, 400, qnorm(.975)),
                   pnorm(.5 / (3 / sqrt(400)) - qnorm(.975)))
})

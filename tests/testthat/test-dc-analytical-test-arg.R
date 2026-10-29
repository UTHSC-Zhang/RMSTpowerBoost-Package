# The dependent-censoring engine exposes the same WLS-vs-sandwich choice as the
# linear one, with the same default: the weighted-least-squares t-test.
dc_test_pilot <- function(n = 500) {
  set.seed(4417)
  d <- data.frame(arm = rep(0:1, each = n / 2), z = stats::rnorm(n, 4, 1))
  event  <- stats::rexp(n, exp(-2 + .15 * d$z - .3 * d$arm))
  censor <- stats::rexp(n, exp(-3 + .4 * d$z))
  d$time <- pmin(event, censor)
  d$status <- as.integer(event <= censor)
  d
}

test_that("the DC model-based SE is the one summary(lm) reports", {
  d <- dc_test_pilot()
  e <- RMSTpowerBoost:::.estimate_dc_params(d, "time", "status", "arm", "z", 2)
  expect_equal(unname(e$se_beta_model_n1 / sqrt(e$n_pilot)),
               unname(summary(e$fit_wls)$coefficients["arm", 2]),
               tolerance = 1e-10)
  # The sandwich and the model-based variance must not coincide: IPCW weights
  # are sampling weights, not inverse-variance weights.
  expect_false(isTRUE(all.equal(e$se_beta_n1, e$se_beta_model_n1)))

  z <- stats::qnorm(.975)
  for (n_arm in c(200, 400)) {
    N <- 2 * n_arm
    s_sw <- e$se_beta_n1 / sqrt(N)
    s_ml <- e$se_beta_model_n1 / sqrt(N)
    b <- abs(e$beta_effect)
    expect_equal(
      unname(RMSTpowerBoost:::.rmst_wald_power(e$beta_effect, e$se_beta_n1, N, z,
                                               e$se_beta_model_n1)),
      unname(stats::pnorm((b - z * s_ml) / s_sw) +
             stats::pnorm((-b - z * s_ml) / s_sw)),
      tolerance = 1e-12)
  }
})

test_that("DC defaults to wls, sandwich is opt-in and more powerful", {
  d <- dc_test_pilot()
  sizes <- c(200, 400)
  default <- DC.power.analytical(d, "time", "status", "arm", sizes, "z", 2)
  wls <- DC.power.analytical(d, "time", "status", "arm", sizes, "z", 2, test = "wls")
  sand <- DC.power.analytical(d, "time", "status", "arm", sizes, "z", 2,
                              test = "sandwich")
  expect_equal(default$results_data$Power, wls$results_data$Power)
  # The WLS test pays for an inflated standard error, so it is the less powerful one.
  expect_true(all(wls$results_data$Power < sand$results_data$Power))

  e <- RMSTpowerBoost:::.estimate_dc_params(d, "time", "status", "arm", "z", 2)
  # Reported SEs follow the selected test.
  expect_equal(wls$model_output$treatment_effect$std_error,
               unname(summary(e$fit_wls)$coefficients["arm", 2]), tolerance = 1e-10)
  expect_equal(sand$model_output$treatment_effect$std_error,
               unname(e$se_beta_n1 / sqrt(e$n_pilot)), tolerance = 1e-10)
  # Both are always available regardless of the choice.
  expect_equal(wls$model_output$variance_components$se_effect_sandwich_n1,
               sand$model_output$variance_components$se_effect_sandwich_n1)
  expect_equal(wls$model_output$variance_components$se_effect_wls_n1,
               sand$model_output$variance_components$se_effect_wls_n1)
  expect_identical(wls$model_output$variance_components$test, "wls")
  expect_identical(sand$model_output$variance_components$test, "sandwich")

  # The coefficient table follows the test too: the WLS t-table under "wls",
  # the robust sandwich table under "sandwich".
  expect_equal(wls$model_output$coefficient_table$std_error,
               unname(summary(e$fit_wls)$coefficients[, 2]), tolerance = 1e-10)
  expect_equal(sand$model_output$coefficient_table$std_error,
               unname(sqrt(diag(e$V_hat_n) / e$n_pilot)), tolerance = 1e-10)
})

test_that("the DC sample-size search honours the selected test", {
  d <- dc_test_pilot()
  args <- list(d, "time", "status", "arm", .8, "z", 2,
               n_start = 100, n_step = 100, max_n_per_arm = 20000)
  ss_default <- suppressWarnings(do.call(DC.ss.analytical, args))
  ss_w <- suppressWarnings(do.call(DC.ss.analytical, c(args, list(test = "wls"))))
  ss_s <- suppressWarnings(do.call(DC.ss.analytical, c(args, list(test = "sandwich"))))
  expect_identical(ss_default$results_data$Required_N_per_Arm,
                   ss_w$results_data$Required_N_per_Arm)
  # A larger required N under the less powerful test.
  expect_gte(ss_w$results_data$Required_N_per_Arm,
             ss_s$results_data$Required_N_per_Arm)
})

test_that("rmst.power and rmst.ss pass test through to the DC engine", {
  d <- dc_test_pilot()
  sand <- DC.power.analytical(d, "time", "status", "arm", c(200, 400), "z", 2,
                              test = "sandwich")
  r <- expect_no_message(rmst.power(survival::Surv(time, status) ~ z, data = d,
                                    arm = "arm", sample_sizes = c(200, 400), L = 2,
                                    dep_cens = TRUE, test = "sandwich"))
  expect_equal(r$results_data$Power, sand$results_data$Power)
  expect_identical(r$model_output$variance_components$test, "sandwich")
  # Left NULL, the engine keeps its own default.
  r0 <- rmst.power(survival::Surv(time, status) ~ z, data = d, arm = "arm",
                   sample_sizes = c(200, 400), L = 2, dep_cens = TRUE)
  expect_identical(r0$model_output$variance_components$test, "wls")

  s <- suppressWarnings(rmst.ss(survival::Surv(time, status) ~ z, data = d,
                                arm = "arm", target_power = .8, L = 2,
                                dep_cens = TRUE, test = "sandwich",
                                n_start = 100, n_step = 100, max_n = 20000))
  expect_identical(s$model_output$variance_components$test, "sandwich")
})

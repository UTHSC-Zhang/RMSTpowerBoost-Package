# The multiplicative engine exposes the same WLS-vs-sandwich choice as the
# linear one, but defaults the other way: no bootstrap engine in the package
# simulates an IPCW-weighted t-test for a stratified model.
ms_pilot <- function(n = 600) {
  set.seed(4417)
  d <- data.frame(arm = rep(0:1, each = n / 2), z = stats::rnorm(n, 4, 1),
                  region = factor(rep(c("A", "B", "C"), length.out = n)))
  event  <- stats::rexp(n, exp(-2 + .15 * d$z - .3 * d$arm))
  censor <- stats::rexp(n, exp(-3 + .4 * d$z + .5 * d$arm))
  d$time <- pmin(event, censor)
  d$status <- as.integer(event <= censor)
  d
}

test_that("the MS model-based SE is the one summary(lm) reports", {
  d <- ms_pilot()
  e <- RMSTpowerBoost:::.estimate_ms_params(d, "time", "status", "arm",
                                            "region", "z", 2)
  expect_equal(unname(e$se_beta_model_n1 / sqrt(e$n_pilot)),
               unname(summary(e$fit_log_lm)$coefficients["arm", 2]),
               tolerance = 1e-10)
  # The sandwich and the model-based variance must not coincide: IPCW weights
  # are sampling weights, not inverse-variance weights.
  expect_false(isTRUE(all.equal(e$se_beta_n1, e$se_beta_model_n1)))

  z <- stats::qnorm(.975)
  for (n_stratum in c(200, 400)) {
    N <- 3 * n_stratum
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

test_that("MS defaults to sandwich, wls is opt-in and less powerful", {
  d <- ms_pilot()
  sizes <- c(100, 200)
  default <- MS.power.analytical(d, "time", "status", "arm", "region", sizes, "z", 2)
  sand <- MS.power.analytical(d, "time", "status", "arm", "region", sizes, "z", 2,
                              test = "sandwich")
  wls <- MS.power.analytical(d, "time", "status", "arm", "region", sizes, "z", 2,
                             test = "wls")
  # Unlike the linear engine, the default here preserves the sandwich numbers.
  expect_equal(default$results_data$Power, sand$results_data$Power)
  # The WLS test pays for an inflated standard error, so it is the less powerful one.
  expect_true(all(wls$results_data$Power < sand$results_data$Power))

  e <- RMSTpowerBoost:::.estimate_ms_params(d, "time", "status", "arm",
                                            "region", "z", 2)
  # Reported SEs follow the selected test.
  expect_equal(wls$model_output$treatment_effect$std_error[1],
               unname(summary(e$fit_log_lm)$coefficients["arm", 2]), tolerance = 1e-10)
  expect_equal(sand$model_output$treatment_effect$std_error[1],
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
               unname(summary(e$fit_log_lm)$coefficients[, 2]), tolerance = 1e-10)
  expect_equal(sand$model_output$coefficient_table$std_error,
               unname(sqrt(diag(e$V_hat_n) / e$n_pilot)), tolerance = 1e-10)
})

test_that("the MS sample-size search honours the selected test", {
  d <- ms_pilot()
  args <- list(d, "time", "status", "arm", "region", .8, "z", 2,
               n_start = 100, n_step = 100, max_n_per_arm = 20000)
  ss_default <- suppressWarnings(do.call(MS.ss.analytical, args))
  ss_s <- suppressWarnings(do.call(MS.ss.analytical, c(args, list(test = "sandwich"))))
  ss_w <- suppressWarnings(do.call(MS.ss.analytical, c(args, list(test = "wls"))))
  expect_identical(ss_default$results_data$Required_N_per_Stratum,
                   ss_s$results_data$Required_N_per_Stratum)
  # A larger required N under the less powerful test.
  expect_gte(ss_w$results_data$Required_N_per_Stratum,
             ss_s$results_data$Required_N_per_Stratum)
})

test_that("the additive engine stays sandwich-only", {
  d <- ms_pilot()
  # No lm() fit means no WLS t-test and nothing to select between.
  expect_false("test" %in% names(formals(additive.power.analytical)))
  expect_false("test" %in% names(formals(additive.ss.analytical)))
  out <- additive.power.analytical(d, "time", "status", "arm", "region", 100, "z", 2)
  expect_null(out$model_output$variance_components$test)
  expect_true(all(is.finite(out$results_data$Power)))
})

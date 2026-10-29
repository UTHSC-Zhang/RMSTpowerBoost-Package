# Compare the dependent-censoring engine to an independently assembled
# covariate-only Cox/Breslow IPCW fit, and pin down that its weights are no
# longer capped or floored.
dc_pilot <- function(n = 500) {
  set.seed(7241)
  d <- data.frame(arm = rep(0:1, each = n / 2), z = stats::rnorm(n, 4, 1))
  event  <- stats::rexp(n, exp(-2 + .15 * d$z - .3 * d$arm))
  # Censoring depends on z only, which is the DC engine's censoring model.
  censor <- stats::rexp(n, exp(-3 + .4 * d$z))
  d$time <- pmin(event, censor)
  d$status <- as.integer(event <= censor)
  d
}

# predict(type = "survival") keeps the baseline and linear-predictor centering
# consistent on its own, so it is an independent route to the same weights.
dc_reference <- function(d, L = 2) {
  cm <- survival::coxph(survival::Surv(time, status == 0) ~ z,
                        data = d, ties = "breslow")
  nd <- d
  nd$time <- pmin(d$time, L)
  g <- as.numeric(stats::predict(cm, newdata = nd, type = "survival"))
  # Eqs. (11)-(12) of Zhang & Schaubel (2024) use the fitted weights as-is.
  w <- as.numeric(d$status == 1 | d$time >= L) / g
  d$Y <- pmin(d$time, L)
  list(w = w, fit = stats::lm(Y ~ arm + z, data = d[w > 0, ], weights = w[w > 0]),
       cox = cm)
}

test_that("DC Cox IPCW matches independent conditional Breslow weights", {
  d <- dc_pilot()
  e <- RMSTpowerBoost:::.estimate_dc_params(d, "time", "status", "arm", "z", 2)
  r <- dc_reference(d)
  expect_equal(e$df$weights, r$w, tolerance = 1e-10)
  expect_equal(unname(e$beta_effect), unname(coef(r$fit)["arm"]), tolerance = 1e-10)

  # The previous construction paired the uncentered baseline hazard with the
  # default, centered, linear predictor. With mean(z) near 4 that misstates G
  # badly, so the corrected weights must be far from it.
  bh <- survival::basehaz(r$cox, centered = FALSE)
  H0 <- stats::stepfun(bh$time, c(0, bh$hazard))(pmin(d$time, 2))
  g_centered <- exp(-H0 * exp(stats::predict(r$cox, newdata = d, type = "lp")))
  w_centered <- as.numeric(d$status == 1 | d$time >= 2) / g_centered
  expect_gt(max(abs(w_centered - r$w)), .1)

  # The weights do not depend on where the covariate is centred.
  shifted <- d
  shifted$z <- shifted$z + 100
  es <- RMSTpowerBoost:::.estimate_dc_params(shifted, "time", "status", "arm", "z", 2)
  expect_equal(es$df$weights, e$df$weights, tolerance = 1e-8)

  # Conditional censoring must not silently revert to pooled KM.
  km <- survival::survfit(survival::Surv(time, status == 0) ~ 1, data = d)
  g <- stats::stepfun(km$time, c(1, km$surv))(pmin(d$time, 2))
  kmw <- as.numeric(d$status == 1 | d$time >= 2) / g
  expect_gt(max(abs(r$w - kmw)), .1)
})

test_that("the intercept-only DC censoring model matches Nelson-Aalen", {
  d <- dc_pilot()
  e <- RMSTpowerBoost:::.estimate_dc_params(d, "time", "status", "arm", NULL, 2)
  # With no covariates the Breslow baseline hazard is the Nelson-Aalen estimate.
  na <- survival::survfit(survival::Surv(time, status == 0) ~ 1, data = d)
  H <- stats::stepfun(na$time, c(0, na$cumhaz))(pmin(d$time, 2))
  wref <- as.numeric(d$status == 1 | d$time >= 2) * exp(H)
  expect_equal(e$df$weights, wref, tolerance = 1e-10)
})

test_that("DC IPCW weights are never capped or truncated", {
  d <- dc_pilot()
  e <- RMSTpowerBoost:::.estimate_dc_params(d, "time", "status", "arm", "z", 2)
  w <- RMSTpowerBoost:::.dc_cox_weights(e$df, "time", "status", "z")
  # A plain numeric vector: there is no cap to report alongside it.
  expect_type(w, "double")
  expect_null(names(w))
  expect_equal(e$df$weights, w, tolerance = 1e-12)
  # The upper tail survives intact, so no winsorizing happened.
  pos <- w[w > 0]
  expect_gt(max(pos), unname(stats::quantile(pos, .99)))
  expect_gt(sum(pos > unname(stats::quantile(pos, .99))), 0)
  expect_null(e$weight_cap)
  # Nothing cap-shaped is reported to the user any more.
  cw <- DC.power.analytical(d, "time", "status", "arm", 100, "z", 2)$model_output$censoring_weights
  expect_named(cw, "raw_summary")
})

test_that("DC weights carry Delta_Y and handle no censoring", {
  d <- dc_pilot()
  e <- RMSTpowerBoost:::.estimate_dc_params(d, "time", "status", "arm", "z", 2)
  expect_true(all(e$df$weights[!e$df$is_complete] == 0))
  expect_true(all(e$df$weights[e$df$is_complete] > 0))

  # A subject censored before L is dropped; one censored at or past L is kept.
  d2 <- d
  d2$time[1:3] <- c(1, 2, 3)
  d2$status[1:3] <- 0
  e2 <- RMSTpowerBoost:::.estimate_dc_params(d2, "time", "status", "arm", "z", 2)
  expect_equal(e2$df$weights[1], 0)
  expect_true(all(e2$df$weights[2:3] > 0))

  # With no censoring, G is identically one and no Cox fit is needed.
  d3 <- d
  d3$status <- 1L
  for (lt in list("z", NULL)) {
    expect_equal(RMSTpowerBoost:::.estimate_dc_params(
      d3, "time", "status", "arm", lt, 2)$df$weights, rep(1, nrow(d3)))
  }
})

test_that("the DC engine reports censoring-model failure", {
  # Every censored subject has z = 1 and every event z = 0, so the censoring
  # Cox likelihood is monotone and its coefficient runs off to infinity.
  d <- data.frame(time = 1:20, status = rep(0:1, each = 10),
                  z = rep(1:0, each = 10), arm = rep(0:1, 10))
  # Previously the bare coxph() call swallowed the non-convergence warning.
  expect_error(
    RMSTpowerBoost:::.estimate_dc_params(d, "time", "status", "arm", "z", 15),
    "censoring model")
})

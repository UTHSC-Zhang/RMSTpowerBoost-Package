# Compare the stratified engines to an independently assembled stratified
# Cox/Breslow IPCW fit, and pin down that the weights are no longer capped.
strat_pilot <- function(n = 600) {
  set.seed(7241)
  d <- data.frame(arm = rep(0:1, each = n / 2), z = stats::rnorm(n, 4, 1),
                  region = factor(rep(c("A", "B", "C"), length.out = n)))
  event  <- stats::rexp(n, exp(-2 + .15 * d$z - .3 * d$arm))
  censor <- stats::rexp(n, exp(-3 + .4 * d$z + .5 * d$arm))
  d$time <- pmin(event, censor)
  d$status <- as.integer(event <= censor)
  d
}

# Uncentered baseline hazard paired with an uncentered linear predictor, by
# stratum, evaluated at the truncated outcome. basehaz() labels strata either
# "level" or "var=level" depending on the survival version; accept both.
strat_reference <- function(d, L = 2, covs = c("arm", "z"), sv = "region") {
  f <- stats::as.formula(paste0("survival::Surv(time, status == 0) ~ ",
                                paste(covs, collapse = " + "),
                                " + survival::strata(", sv, ")"))
  cm <- survival::coxph(f, data = d, ties = "breslow")
  bh <- survival::basehaz(cm, centered = FALSE)
  y  <- pmin(d$time, L)
  H  <- numeric(nrow(d))
  lab <- as.character(bh$strata)
  for (st in unique(as.character(d[[sv]]))) {
    i <- as.character(d[[sv]]) == st
    j <- lab == st | lab == paste0(sv, "=", st)
    H[i] <- stats::stepfun(bh$time[j], c(0, bh$hazard[j]))(y[i])
  }
  lp <- stats::predict(cm, newdata = d, type = "lp", reference = "zero")
  # Eqs. (11)-(12) of Zhang & Schaubel (2024) use the fitted weights as-is.
  g <- exp(-H * exp(lp))
  list(w = unname(as.numeric(d$status == 1 | d$time >= L) / g), fit = cm)
}

test_that("stratified Cox IPCW matches independent stratified Breslow weights", {
  d <- strat_pilot()
  wref <- strat_reference(d)$w
  ea <- RMSTpowerBoost:::.estimate_additive_params(d, "time", "status", "arm",
                                                   "region", "z", 2)
  em <- RMSTpowerBoost:::.estimate_ms_params(d, "time", "status", "arm",
                                             "region", "z", 2)
  expect_equal(ea$df$weights, wref, tolerance = 1e-10)
  expect_equal(em$df$weights, wref, tolerance = 1e-10)
  # Conditional, stratified censoring must not silently revert to pooled KM.
  km <- survival::survfit(survival::Surv(time, status == 0) ~ 1, data = d)
  g <- stats::stepfun(km$time, c(1, km$surv))(pmin(d$time, 2))
  kmw <- as.numeric(d$status == 1 | d$time >= 2) / g
  expect_gt(max(abs(wref - kmw)), .1)
})

test_that("stratified IPCW weights are never capped or truncated", {
  d <- strat_pilot()
  ea <- RMSTpowerBoost:::.estimate_additive_params(d, "time", "status", "arm",
                                                   "region", "z", 2)
  em <- RMSTpowerBoost:::.estimate_ms_params(d, "time", "status", "arm",
                                             "region", "z", 2)
  # The upper tail survives intact, so no winsorizing happened.
  for (w in list(ea$df$weights, em$df$weights)) {
    pos <- w[w > 0]
    expect_gt(max(pos), unname(stats::quantile(pos, .99)))
    expect_gt(sum(pos > unname(stats::quantile(pos, .99))), 0)
  }
  # There is no cap to report alongside the weights any more.
  expect_null(ea$weight_cap)
  expect_null(em$weight_cap)
  # Nothing cap-shaped is reported to the user any more.
  pa <- additive.power.analytical(d, "time", "status", "arm", "region", 100, "z", 2)
  pm <- MS.power.analytical(d, "time", "status", "arm", "region", 100, "z", 2)
  expect_named(pa$model_output$censoring_weights, "raw_summary")
  expect_named(pm$model_output$censoring_weights, "raw_summary")
})

test_that("stratified weights carry Delta_Y and handle no censoring", {
  d <- strat_pilot()
  ea <- RMSTpowerBoost:::.estimate_additive_params(d, "time", "status", "arm",
                                                   "region", "z", 2)
  # Both engines now receive Delta_Y * W_hat, so incomplete subjects weigh zero
  # and every observed truncated outcome keeps a strictly positive weight.
  expect_true(all(ea$df$weights[!ea$df$is_complete] == 0))
  expect_true(all(ea$df$weights[ea$df$is_complete] > 0))
  em <- RMSTpowerBoost:::.estimate_ms_params(d, "time", "status", "arm",
                                             "region", "z", 2)
  expect_true(all(em$df$weights[!em$df$is_complete] == 0))

  # A subject censored before L is dropped; one censored at or past L is kept.
  d2 <- d
  d2$time[1:3] <- c(1, 2, 3)
  d2$status[1:3] <- 0
  e2 <- RMSTpowerBoost:::.estimate_additive_params(d2, "time", "status", "arm",
                                                   "region", "z", 2)
  expect_equal(e2$df$weights[1], 0)
  expect_true(all(e2$df$weights[2:3] > 0))

  # With no censoring, G is identically one and no Cox fit is needed.
  d3 <- d
  d3$status <- 1L
  expect_equal(unname(RMSTpowerBoost:::.estimate_additive_params(
    d3, "time", "status", "arm", "region", "z", 2)$df$weights), rep(1, nrow(d3)))
  expect_equal(unname(RMSTpowerBoost:::.estimate_ms_params(
    d3, "time", "status", "arm", "region", "z", 2)$df$weights), rep(1, nrow(d3)))
})

test_that("the additive engine now reports censoring-model failure like MS does", {
  d <- strat_pilot()[1:6, ]
  # Previously the bare coxph() call swallowed the non-convergence warning.
  expect_error(
    RMSTpowerBoost:::.estimate_additive_params(d, "time", "status", "arm",
                                               "region", "z", 2))
})

# Compare the package to an independently assembled Cox/Breslow IPCW fit.
cox_pilot <- function() {
  set.seed(7241)
  d <- data.frame(arm = rep(0:1, each = 200), z = rnorm(400, 4, 1))
  event <- rexp(400, exp(-2 + .15 * d$z - .3 * d$arm))
  censor <- rexp(400, exp(-3 + .4 * d$z + .5 * d$arm))
  d$time <- pmin(event, censor)
  d$status <- as.integer(event <= censor)
  d
}

cox_reference <- function(d, L = 2) {
  cm <- survival::coxph(survival::Surv(time, status == 0) ~ arm + z,
                        data = d, ties = "breslow")
  baseline <- survival::basehaz(cm, centered = FALSE)
  hazard <- stats::stepfun(baseline$time, c(0, baseline$hazard))(pmin(d$time, L))
  lp <- stats::predict(cm, newdata = d, type = "lp", reference = "zero")
  # Eqs. (11)-(12) of Zhang & Schaubel (2024) use the fitted weights as-is.
  w <- as.numeric(d$status == 1 | d$time >= L) * exp(hazard * exp(lp))
  d$Y <- pmin(d$time, L)
  fit <- stats::lm(Y ~ factor(arm) + z, data = d, weights = w)
  list(w = w, fit = fit)
}

test_that("linear Cox IPCW matches independent conditional Breslow weights", {
  d <- cox_pilot()
  e <- RMSTpowerBoost:::.estimate_linear_params(d, "time", "status", "arm", "z", 2)
  r <- cox_reference(d)
  expect_equal(e$df$weights, unname(r$w), tolerance = 1e-10)
  expect_equal(unname(e$beta_effect), unname(coef(r$fit)[2]), tolerance = 1e-10)
  # Conditional censoring must not silently revert to pooled KM.
  km <- survival::survfit(survival::Surv(time, status == 0) ~ 1, data = d)
  g <- stats::stepfun(km$time, c(1, km$surv))(pmin(d$time, 2))
  kmw <- as.numeric(d$status == 1 | d$time >= 2) / g
  expect_gt(max(abs(r$w - kmw)), .1)
  shifted <- d
  shifted$z <- shifted$z + 100
  es <- RMSTpowerBoost:::.estimate_linear_params(shifted, "time", "status", "arm", "z", 2)
  expect_equal(es$df$weights, e$df$weights, tolerance = 1e-8)
  dfactor <- d
  dfactor$arm <- factor(dfactor$arm)
  ef <- RMSTpowerBoost:::.estimate_linear_params(dfactor, "time", "status", "arm", "z", 2)
  expect_equal(ef$df$weights, e$df$weights, tolerance = 1e-10)
})

test_that("linear IPCW weights are never capped or truncated", {
  d <- cox_pilot()
  w <- RMSTpowerBoost:::.linear_cox_weights(d, "time", "status", "arm", "z", 2)
  # A plain numeric vector: there is no cap to report alongside it.
  expect_type(w, "double")
  expect_null(names(w))
  expect_length(w, nrow(d))
  # The upper tail survives intact, so no winsorizing happened.
  pos <- w[w > 0]
  expect_gt(max(pos), unname(stats::quantile(pos, .99)))
  expect_gt(sum(pos > unname(stats::quantile(pos, .99))), 0)
  e <- RMSTpowerBoost:::.estimate_linear_params(d, "time", "status", "arm", "z", 2)
  expect_equal(e$df$weights, w, tolerance = 1e-12)
  expect_null(e$weight_cap)
  # Nothing cap-shaped is reported to the user any more.
  cw <- linear.power.analytical(d, "time", "status", "arm", 100, "z", 2)$model_output$censoring_weights
  expect_named(cw, "raw_summary")
})

test_that("linear Cox weights retain horizon observations and handle no censoring", {
  d <- cox_pilot()
  d$time[1:3] <- c(1, 2, 3)
  d$status[1:3] <- 0
  w <- RMSTpowerBoost:::.linear_cox_weights(d, "time", "status", "arm", "z", 2)
  expect_equal(w[1], 0)
  expect_true(all(w[2:3] > 0))
  d$status <- 1L
  expect_equal(RMSTpowerBoost:::.linear_cox_weights(d, "time", "status", "arm", "z", 2),
               rep(1, nrow(d)))
})

test_that("power and sample-size bootstrap refit the same conditional Cox model", {
  d <- cox_pilot()
  set.seed(884)
  sampled <- do.call(rbind, lapply(split(d, d$arm), function(g)
    g[sample(seq_len(nrow(g)), 150, replace = TRUE), ]))
  r <- summary(cox_reference(sampled)$fit)$coefficients[2, ]
  set.seed(884)
  p <- linear.power.boot(d, "time", "status", "arm", 150, "z", 2, n_sim = 1)
  draws <- p$model_output$simulation_draws
  expect_equal(draws$estimate, unname(r[1]), tolerance = 1e-10)
  expect_equal(draws$std_error, unname(r[2]), tolerance = 1e-10)
  expect_equal(draws$p_value, unname(r[4]), tolerance = 1e-10)
  set.seed(884)
  s <- suppressWarnings(linear.ss.boot(d, "time", "status", "arm", .8, "z", 2,
                                      n_sim = 1, n_start = 150, max_n_per_arm = 150))
  expect_equal(s$results_plot$data$Power, p$results_data$Power)
})

# The GAM bootstrap keeps its jackknife pseudo-observations; these tests pin
# down that they are computed exactly and that each replicate tests exactly the
# treatment coefficient, with numeric stratum codes treated as strata.
gam_pilot <- function(n = 160) {
  set.seed(5)
  d <- data.frame(A = rep(0:1, each = n / 2), Age = stats::rnorm(n, 60, 8),
                  Z = stats::rnorm(n), S = rep(1:3, length.out = n))
  event <- stats::rexp(n, .4 * exp(-.35 * d$A + .06 * (d$Age - 60) + .5 * sin(2 * d$Z)))
  censor <- stats::rexp(n, .15)
  d$time <- pmin(event, censor)
  d$status <- as.integer(event <= censor)
  d
}

# Kaplan-Meier written out from its product-limit definition, without survfit,
# and its area up to L summed rectangle by rectangle.
km_rmst_reference <- function(time, status, L) {
  surv <- 1; area <- 0; last <- 0
  for (t in sort(unique(time[status == 1]))) {
    if (t >= L) break
    area <- area + surv * (t - last)
    surv <- surv * (1 - sum(time == t & status == 1) / sum(time >= t))
    last <- t
  }
  area + surv * (L - last)
}
pseudo_reference <- function(time, status, L) {
  n <- length(time)
  loo <- vapply(seq_len(n), function(i) km_rmst_reference(time[-i], status[-i], L), 0)
  n * km_rmst_reference(time, status, L) - (n - 1) * loo
}

test_that("pseudo-observations integrate the Kaplan-Meier curve exactly", {
  set.seed(11)
  time <- round(stats::rexp(120, .5), 1) + .1
  status <- stats::rbinom(120, 1, .7)
  # Tied times, an event and a censoring exactly at L = 2, and follow-up past L.
  time[1:4] <- c(2, 2, 5, 1.5)
  status[1:4] <- c(1, 0, 1, 1)
  for (L in c(2, 3.3)) {
    expect_equal(RMSTpowerBoost:::.rmst_jackknife_pseudo_obs(time, status, L),
                 pseudo_reference(time, status, L), tolerance = 1e-10)
  }
  # With all follow-up beyond L, every pseudo-observation is exactly L.
  expect_equal(RMSTpowerBoost:::.rmst_jackknife_pseudo_obs(c(5, 6, 7), c(1, 0, 1), 2),
               rep(2, 3))
})

test_that("GAM tests exactly the treatment coefficient", {
  d <- gam_pilot()
  set.seed(31)
  draws <- GAM.power.boot(d, "time", "status", "A", sample_sizes = 60,
                          linear_terms = "Age", L = 3, n_sim = 3)$model_output$simulation_draws
  # Rebuild the same resamples on the same future.seed streams, then fit
  # independently and read the treatment row by name.
  set.seed(31)
  ref <- do.call(rbind, future.apply::future_lapply(1:3, function(j) {
    boot <- do.call(rbind, lapply(split(d, d$A), function(g)
      g[sample(1:nrow(g), size = 60, replace = TRUE), ]))
    boot$pseudo_obs <- NA_real_
    for (a in 0:1) {
      i <- boot$A == a
      boot$pseudo_obs[i] <- pseudo_reference(boot$time[i], boot$status[i], 3)
    }
    boot$A <- factor(boot$A, levels = c(0, 1))
    pt <- mgcv::summary.gam(mgcv::gam(pseudo_obs ~ A + Age, data = boot))$p.table
    c(estimate = pt["A1", "Estimate"], p = pt["A1", "Pr(>|t|)"],
      p_age = pt["Age", "Pr(>|t|)"])
  }, future.seed = TRUE))
  expect_equal(draws$estimate, unname(ref[, "estimate"]), tolerance = 1e-8)
  expect_equal(draws$p_value, unname(ref[, "p"]), tolerance = 1e-8)
  # The Age row also begins with "A"; whenever it is the more significant of
  # the two, a prefix match with the minimum p-value would have reported it.
  expect_true(any(ref[, "p_age"] < ref[, "p"]))
})

test_that("numeric stratum codes give stratum intercepts, as factors do", {
  d <- gam_pilot()
  d$S_factor <- factor(d$S)
  run <- function(sv) {
    set.seed(32)
    GAM.power.boot(d, "time", "status", "A", strata_var = sv, sample_sizes = 40,
                   linear_terms = "Z", L = 3, n_sim = 3)$model_output$simulation_draws
  }
  numeric_codes <- run("S")
  factor_codes <- run("S_factor")
  expect_equal(numeric_codes$estimate, factor_codes$estimate, tolerance = 1e-10)
  expect_equal(numeric_codes$p_value, factor_codes$p_value, tolerance = 1e-10)
})

test_that("GAM.ss.boot runs the same replicates as GAM.power.boot", {
  d <- gam_pilot()
  set.seed(33)
  p <- GAM.power.boot(d, "time", "status", "A", sample_sizes = 50,
                      linear_terms = "Age", L = 3, n_sim = 10)
  set.seed(33)
  s <- suppressWarnings(GAM.ss.boot(d, "time", "status", "A", target_power = .999,
                                    linear_terms = "Age", L = 3, n_sim = 10,
                                    n_start = 50, max_n_per_arm = 50))
  expect_equal(s$results_plot$data$Power, p$results_data$Power)
})

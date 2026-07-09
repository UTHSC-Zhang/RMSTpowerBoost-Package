# Validation script: compare package linear IPCW to Miscellaneous/test.R reference
# Run manually (excluded from R CMD check via .Rbuildignore); requires truncnorm.
suppressPackageStartupMessages({
  library(truncnorm)
  library(survival)
  if (dir.exists("RMSTpowerBoost-Package")) {
    devtools::load_all("RMSTpowerBoost-Package")
  } else if (file.exists("DESCRIPTION")) {
    devtools::load_all(".")
  } else {
    library(RMSTpowerBoost)
  }
})

generateData <- function(n, tau, beta0, beta1, beta2, lambda) {
  Z <- truncnorm::rtruncnorm(n, a = -3, b = 3, mean = 0, sd = 1)
  A <- rbinom(n, 1, 0.5)
  e <- rnorm(n, mean = 0, sd = 0.5)
  D <- exp(beta0 + beta1 * Z + beta2 * A + e)
  C <- rexp(n, rate = lambda)
  X <- pmin(D, C, tau)
  deltaX <- ifelse(D < pmin(C, tau), 1, 0)
  cbind.data.frame(Z, A, D, C, X, deltaX)
}

analyze_ref <- function(dat, L) {
  dat$Y <- pmin(dat$X, L)
  dat$deltaY <- ifelse(dat$deltaX == 1, 1, ifelse(L < dat$X, 1, 0))
  modC <- coxph(Surv(X, 1 - deltaX) ~ Z + A, data = dat, ties = "breslow")
  dat.new <- dat
  dat.new$X <- dat.new$Y
  dat$W <- dat$deltaY / predict(modC, newdata = dat.new, type = "survival")
  mod <- lm(Y ~ Z + A, data = dat, weights = W)
  summary(mod)$coefficients[3, ]
}

calcPower_ref <- function(n.sim, n, tau, L, beta0, beta1, beta2, lambda) {
  res <- rep(NA, n.sim)
  for (i in seq_len(n.sim)) {
    dat <- generateData(n, tau, beta0, beta1, beta2, lambda)
    p <- analyze_ref(dat, L = L)[4]
    res[i] <- ifelse(p < 0.05, 1, 0)
  }
  mean(res)
}

set.seed(123)
tau <- 7.5
beta0 <- 0.25
beta1 <- 0.05
beta2 <- 0.05
lambda <- 1.2
L <- 1

dat <- generateData(n = 4000, tau = tau, beta0 = beta0, beta1 = beta1, beta2 = beta2, lambda = lambda)
ref_coef <- analyze_ref(dat, L = L)

pkg <- linear.power.analytical(
  pilot_data = dat,
  time_var = "X",
  status_var = "deltaX",
  arm_var = "A",
  sample_sizes = 2000,
  linear_terms = "Z",
  L = L,
  alpha = 0.05
)

pkg_coef <- pkg$model_output$coefficient_table
pkg_arm <- pkg_coef[pkg_coef$term == "factor(A)1", ]

cat("Reference A estimate:", ref_coef[1], "SE:", ref_coef[2], "p:", ref_coef[4], "\n")
cat("Package A estimate:", pkg_arm$estimate, "SE:", pkg_arm$std_error, "p:", pkg_arm$p_value, "\n")
cat("Package analytic power at n=2000/arm:", pkg$results_data$Power[1], "\n")

ref_power <- calcPower_ref(
  n.sim = 200,
  n = 4000,
  tau = tau,
  L = L,
  beta0 = beta0,
  beta1 = beta1,
  beta2 = beta2,
  lambda = lambda
)
cat("Reference simulation power (200 reps):", ref_power, "\n")

boot <- linear.power.boot(
  pilot_data = dat,
  time_var = "X",
  status_var = "deltaX",
  arm_var = "A",
  sample_sizes = 2000,
  linear_terms = "Z",
  L = L,
  alpha = 0.05,
  n_sim = 200
)
cat("Package bootstrap power at n=2000/arm:", boot$results_data$Power[1], "\n")

# Alignment checks
stopifnot(abs(ref_coef[1] - pkg_arm$estimate) < 0.01)
stopifnot(abs(ref_power - 0.24) < 0.1)
cat("Validation passed: pilot effect matches test.R; simulation power near 0.24.\n")

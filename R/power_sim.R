library(truncnorm)
library(survival)

###########################
########## Aim 1 ##########
###########################

# ~ 20% deaths observed, ~ 80% censoring
generateData <- function(n = 4500, tau = 7.5, beta0 = 2.5, beta1 = 0.5, lambda = 0.7){
  Z <- rtruncnorm(n, a = -3, b = 3, mean = 0, sd = 1)
  e <- rnorm(n, mean = 0, sd = 0.5)
  D <- beta0 + beta1 * Z + e
  C <- rexp(n, rate = lambda)
  X <- pmin(D, C, tau)
  deltaX <- ifelse(D < pmin(C, tau), 1, 0)
  print(paste0(signif(sum(deltaX)/n*100, digits = 2), "% deaths observed"), sep = "")
  dat <- cbind.data.frame(Z, D, C, X, deltaX)
  return(dat)
}

analyze <- function(dat, L = 5){
  dat$Y <- pmin(dat$X, L)
  dat$deltaY <- ifelse(dat$deltaX == 1, 1, ifelse(L < dat$X, 1, 0))
  modC <- coxph(Surv(X, 1 - deltaX) ~ Z, data = dat, ties = "breslow")
  dat.new = dat
  dat.new$X <- dat.new$Y
  dat$W <- dat$deltaY/predict(modC, newdata = dat.new, type = "survival")
  mod <- lm(Y ~ Z, data = dat, weights = W)
  p <- summary(mod)$coefficients[2,4]
  return(p)
}

calcPower <- function(n.sim = 10000, n = 4500, tau = 7.5, L = 5, beta0 = 2.5, beta1 = 0.05, lambda = 0.7){
  res <- rep(NA, n.sim)
  for (i in 1:n.sim){
    dat <- generateData(n = n, tau = tau, beta0 = beta0, beta1 = beta1, lambda = lambda)
    p <- analyze(dat, L = L)
    res[i] <- ifelse(p < 0.05, 1, 0)
  }
  return(sum(res)/n.sim)
}

calcPower(n.sim = 10000, n = 4000, L = 5, beta1 = 0.05) # 0.75
calcPower(n.sim = 10000, n = 4000, L = 5, beta1 = 0.1) # 0.99
calcPower(n.sim = 10000, n = 4000, L = 5, beta1 = 0.15) # 1
calcPower(n.sim = 10000, n = 4000, L = 5, beta1 = 0.2) # 1

###########################
########## Aim 3 ##########
###########################

# ~ 40% deaths observed, ~ 60% censoring
dat <- generateData(n = 3000, tau = 7.5, beta0 = 3.5, beta1 = 0.05, lambda = 0.27)
calcPower(n.sim = 10000, n = 3000, L = 5, beta0 = 3.5, beta1 = 0.05, lambda = 0.27) # 0.9187
calcPower(n.sim = 10000, n = 3000, L = 5, beta0 = 3.5, beta1 = 0.1, lambda = 0.27) # 1
calcPower(n.sim = 10000, n = 3000, L = 5, beta0 = 3.5, beta1 = 0.15, lambda = 0.27) # 1
calcPower(n.sim = 10000, n = 3000, L = 5, beta0 = 3.5, beta1 = 0.2, lambda = 0.27) # 1


#' Estimate Restricted Mean Survival Time (RMST)
#'
#' Compute the RMST estimate from observed survival times and censoring indicators
#' up to a specified truncation time tau.
#'
#' @param time Numeric vector of observed survival or censoring times.
#' @param status Numeric vector of event indicators (1 = event, 0 = censored).
#' @param tau Numeric truncation time (max time to integrate survival curve).
#'
#' @return Numeric estimate of RMST up to time tau.
#' @importFrom survival Surv survfit
#' @export
estimate_rmst <- function(time, status, tau) {
   requireNamespace("survival", quietly = TRUE)

   fit <- survival::survfit(survival::Surv(time, status) ~ 1)
   surv_times <- fit$time
   surv_prob <- fit$surv

   # Add zero time point
   surv_times <- c(0, surv_times)
   surv_prob <- c(1, surv_prob)

   # Truncate survival times at tau
   surv_times_trunc <- c(surv_times[surv_times <= tau], tau)
   surv_prob_trunc <- surv_prob[c(surv_times <= tau, TRUE)]

   # Calculate area under survival curve by trapezoidal rule
   dt <- diff(surv_times_trunc)
   rmst <- sum(dt * surv_prob_trunc[-length(surv_prob_trunc)])
   return(rmst)
}

#' Simulate Stratified Survival Data with Exponential Hazards and Uniform Censoring
#'
#' Generate a stratified survival dataset with specified sample sizes, hazard rates,
#' and censoring rate.
#'
#' @param n_per_stratum Numeric vector with sample sizes for each stratum.
#' @param hazard_rates Numeric vector of hazard rates per stratum.
#' @param censoring_rate Numeric scalar for maximum censoring time (uniform censoring).
#' @param strata_names Optional character vector naming strata.
#'
#' @return A data.frame with columns: time, status, stratum.
#' @export


simulate_stratified_survival <- function(n_per_stratum, hazard_rates, censoring_rate, strata_names = NULL) {
   if(length(n_per_stratum) != length(hazard_rates)) {
      stop("Lengths of n_per_stratum and hazard_rates must be equal.")
   }

   n_strata <- length(n_per_stratum)
   if (is.null(strata_names)) {
      strata_names <- paste0("Stratum", seq_len(n_strata))
   }

   df_list <- vector("list", n_strata)
   for (i in seq_len(n_strata)) {
      n <- n_per_stratum[i]
      hazard <- hazard_rates[i]

      true_times <- rexp(n, rate = hazard)
      censoring_times <- runif(n, 0, censoring_rate)

      observed_times <- pmin(true_times, censoring_times)
      status <- as.numeric(true_times <= censoring_times)

      df_list[[i]] <- data.frame(
         time = observed_times,
         status = status,
         stratum = factor(rep(strata_names[i], n))
      )
   }

   df <- do.call(rbind, df_list)
   rownames(df) <- NULL
   return(df)
}

#' Calculate Required Sample Size per Stratum for Detecting RMST Difference
#'
#' Calculate sample sizes needed per stratum to detect a specified RMST effect size
#' with given alpha and power, using stratified survival data.
#'
#' @param data Data frame with columns: time, status, stratum (factor).
#' @param tau Numeric truncation time for RMST calculation.
#' @param effect_size Numeric expected difference in RMST between groups.
#' @param alpha Numeric Type I error rate (default 0.05).
#' @param power Numeric desired power (default 0.8).
#' @param n_per_stratum Optional numeric vector of current sample sizes per stratum (to estimate variance).
#'
#' @return Numeric vector of required sample sizes per stratum.
#' @importFrom stats qnorm var
#' @export
rmst_sample_size <- function(data, tau, effect_size, alpha = 0.05, power = 0.8, n_per_stratum = NULL) {
   strata_levels <- levels(data$stratum)
   if (is.null(strata_levels)) {
      stop("Data must have a factor 'stratum' column with levels.")
   }

   n_strata <- length(strata_levels)

   # Estimate RMST variance per stratum via bootstrap or analytical approx
   var_vec <- numeric(n_strata)
   n_vec <- numeric(n_strata)

   for (i in seq_along(strata_levels)) {
      sub_data <- subset(data, stratum == strata_levels[i])
      n_i <- nrow(sub_data)
      n_vec[i] <- ifelse(is.null(n_per_stratum), n_i, n_per_stratum[i])

      # Bootstrap variance estimate for RMST
      boot_rmst <- replicate(100, {
         idx <- sample(seq_len(n_i), n_i, replace = TRUE)
         estimate_rmst(sub_data$time[idx], sub_data$status[idx], tau)
      })
      var_vec[i] <- var(boot_rmst)
   }

   # Compute pooled variance
   pooled_var <- pooled_variance(var_vec, n_vec)

   # Calculate z values
   z_alpha <- qnorm(1 - alpha / 2)
   z_beta <- qnorm(power)

   # Sample size formula for detecting difference in RMST effect size
   n_required <- (z_alpha + z_beta)^2 * pooled_var / effect_size^2

   # Allocate sample size evenly across strata (can be modified)
   n_per_stratum_req <- rep(ceiling(n_required / n_strata), n_strata)
   names(n_per_stratum_req) <- strata_levels

   return(n_per_stratum_req)
}


#' Power Calculation for Difference in RMST
#'
#' Computes the power of a two-sample test comparing RMST between two groups,
#' assuming normal approximation.
#'
#' @param delta Expected difference in RMST between treatment and control groups.
#' @param sigma Estimated common standard deviation of RMST.
#' @param n_total Total sample size (evenly split between groups).
#' @param alpha Type I error rate (default is 0.05).
#'
#' @return Power of the test (probability of rejecting null if true difference is `delta`)
#' @examples
#' RMSt.power(delta = 1.2, sigma = 2, n_total = 100, alpha = 0.05)
#' @export
RMSt.power <- function(delta, sigma, n_total, alpha = 0.05) {
   if (n_total %% 2 != 0) {
      warning("Total sample size n_total should be even. Rounding down.")
      n_total <- floor(n_total / 2) * 2
   }

   n_per_group <- n_total / 2
   se <- sqrt(2) * sigma / sqrt(n_per_group)
   z_alpha <- qnorm(1 - alpha / 2)
   z_power <- (abs(delta) / se) - z_alpha

   power <- pnorm(z_power)
   return(power)
}
#' Sample Size Calculation Based on RMST
#'
#' Computes the required sample size per group to achieve a specified power for detecting a given difference in RMST.
#'
#' @param delta The effect size (difference in RMST between groups).
#' @param sigma2 The variance of the RMST difference estimator.
#' @param alpha Significance level (default: 0.05).
#' @param power Desired power (default: 0.80).
#'
#' @return Required sample size per group.
#' @examples
#' RMSt.sample.size(delta = 0.5, sigma2 = 0.25, alpha = 0.05, power = 0.8)
#' @export
RMSt.sample.size <- function(delta, sigma2, alpha = 0.05, power = 0.8) {
   if (delta == 0) stop("Effect size delta cannot be zero.")

   z_alpha <- qnorm(1 - alpha / 2)
   z_beta <- qnorm(power)

   n_per_group <- (2 * sigma2 * (z_alpha + z_beta)^2) / delta^2
   return(ceiling(n_per_group))
}



#' Pooled Variance Calculation
#'
#' Calculate the pooled variance from a vector of variances and corresponding sample sizes.
#' Used internally for stratified sample size estimation.
#'
#' @param var_vec Numeric vector of variances for each stratum.
#' @param n_vec Numeric vector of sample sizes for each stratum.
#'
#' @return Numeric pooled variance.
#' @keywords internal
#' @export
pooled_variance <- function(var_vec, n_vec) {
   if (length(var_vec) != length(n_vec)) {
      stop("Length of var_vec and n_vec must be equal.")
   }
   numerator <- sum((n_vec - 1) * var_vec)
   denominator <- sum(n_vec) - length(n_vec)
   pooled_var <- numerator / denominator
   return(pooled_var)
}

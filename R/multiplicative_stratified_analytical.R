# Power Calculation -------------------------------------------------------

#' @title Analyze Power for a Multiplicative Stratified RMST Model (Analytic)
#' @description Performs power analysis for a multiplicative, stratified RMST model using an
#'   analytic method based on the work of Wang et al. (2019).
#'
#' @details
#' This function is based on the method for evaluating center (stratum) effects
#' using a multiplicative model for RMST: \eqn{\mu_{ij} = \mu_{0j} \exp\{\beta'Z_i\}}.
#' The method uses IPCW with a stratified Cox model for the censoring distribution.
#'
#' The formal estimation of \eqn{\beta} requires an iterative solver for a complex
#' estimating equation (Equation 8 in the paper). As this is computationally intensive,
#' this implementation uses a practical approximation by fitting a weighted log-linear
#' model (`lm(log(Y_rmst) ~ ...)`), which is conceptually similar and provides robust
#' estimates for the effect size and its variance.
#'
#' The power calculation relies on the asymptotic variance of the log-RMST ratio
#' estimator, \eqn{\hat{\beta}}. The variance is derived from the robust variance-covariance
#' matrix of the `lm` fit, which serves as a proxy for the formal sandwich estimator
#' \eqn{A^{-1}B(A^{-1})'} described in Theorem 1 of the paper.
#'
#' @param pilot_data A `data.frame` with pilot study data.
#' @param time_var A character string for the time-to-event variable.
#' @param status_var A character string for the event status variable (1=event, 0=censored).
#' @param arm_var A character string for the treatment arm variable (1=treatment, 0=control).
#' @param strata_var A character string for the stratification variable.
#' @param sample_sizes A numeric vector of sample sizes *per stratum* to calculate power for.
#' @param linear_terms An optional character vector of other covariate names.
#' @param tau The numeric value for the RMST truncation time.
#' @param alpha The significance level (Type I error rate).
#'
#' @return A `list` containing:
#' \item{results_data}{A `data.frame` with sample sizes and corresponding powers.}
#' \item{results_plot}{A `ggplot` object visualizing the power curve.}
#' @export
#' @importFrom survival Surv coxph basehaz
#' @importFrom stats lm as.formula complete.cases na.omit sd quantile pnorm qnorm coef vcov predict
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline labs theme_minimal ylim
#' @examples
#' set.seed(123)
#' pilot_df_strat <- data.frame(
#'  time = rexp(120, 0.15),
#'  status = rbinom(120, 1, 0.6),
#'  arm = rep(0:1, each = 60),
#'  region = factor(rep(c("A", "B", "C"), each = 40))
#' )
#' pilot_df_strat$time[pilot_df_strat$arm == 1] <- pilot_df_strat$time[pilot_df_strat$arm == 1] * 1.5
#'
#' power_results <- MS.power.analytical(
#'  pilot_data = pilot_df_strat,
#'  time_var = "time", status_var = "status", arm_var = "arm", strata_var = "region",
#'  sample_sizes = c(50, 75, 100),
#'  tau = 10, alpha = 0.05
#' )
#' print(power_results$results_data)
MS.power.analytical <- function(pilot_data, time_var, status_var, arm_var, strata_var,
                                sample_sizes, linear_terms = NULL, tau, alpha = 0.05) {
   # --- 1. Estimate Parameters from Pilot Data ---
   cat("--- Estimating parameters from pilot data (log-linear approximation)... ---\n")

   core_vars <- c(time_var, status_var, arm_var, strata_var)
   covariates <- c(arm_var, linear_terms)
   all_vars <- c(core_vars, linear_terms)
   df <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]
   n_pilot <- nrow(df)

   # Prepare data for IPCW
   df$Y_rmst <- pmin(df[[time_var]], tau)
   df$is_event <- df[[status_var]] == 1

   # Model for censoring (stratified by strata_var)
   cens_formula <- stats::as.formula(paste0("survival::Surv(Y_rmst, is_event == 0) ~ ",
                                            paste(covariates, collapse = " + "),
                                            " + strata(", strata_var, ")"))
   fit_cens <- survival::coxph(cens_formula, data = df, ties = "breslow")

   # Calculate IPCW weights
   bh_cens <- survival::basehaz(fit_cens, centered = FALSE)
   df$H_cens <- 0
   unique_strata_from_bh <- unique(bh_cens$strata)
   for(st in unique(df[[strata_var]])){
      st_label <- paste0(strata_var, "=", st)
      is_stratum <- df[[strata_var]] == st
      if (st_label %in% unique_strata_from_bh) {
         is_bh_stratum <- bh_cens$strata == st_label
         if(sum(is_bh_stratum) > 0){
            H_st <- stats::stepfun(bh_cens$time[is_bh_stratum], c(0, bh_cens$hazard[is_bh_stratum]))(df$Y_rmst[is_stratum])
            df$H_cens[is_stratum] <- H_st
         }
      }
   }
   df$weights <- exp(df$H_cens * exp(predict(fit_cens, newdata = df, type="lp", reference="zero")))

   # Stabilize weights
   finite_weights <- df$weights[is.finite(df$weights)]
   if (length(finite_weights) > 0) {
      weight_cap <- stats::quantile(finite_weights, probs = 0.99, na.rm = TRUE)
      df$weights[df$weights > weight_cap] <- weight_cap
   }
   df$weights[!is.finite(df$weights)] <- 0
   df$w_delta <- df$weights * df$is_event
   df$w_delta[is.na(df$w_delta)] <- 0

   # Use a log-linear model as a practical approximation to the multiplicative model
   fit_data <- df[df$w_delta > 0 & df$Y_rmst > 0, ]
   fit_weights <- fit_data$w_delta

   log_model_formula <- stats::as.formula(paste("log(Y_rmst) ~", paste(covariates, collapse=" + "), "+", strata_var))
   message("Approximation Model: ", deparse(log_model_formula))

   if(nrow(fit_data) < (length(covariates) + length(unique(fit_data[[strata_var]])))) {
      stop("Not enough data points to fit the approximation model after filtering.", call. = FALSE)
   }

   fit_log_lm <- lm(log_model_formula, data = fit_data, weights = fit_weights)

   beta_summary <- coef(summary(fit_log_lm))
   beta_effect <- beta_summary[arm_var, "Estimate"]

   # --- 2. Calculate Asymptotic Variance from the approximate model ---
   V_hat <- vcov(fit_log_lm) * nrow(fit_data) # Scale to Var(sqrt(n)*beta)
   var_beta_n1 <- V_hat[arm_var, arm_var]
   se_beta_n1 <- sqrt(var_beta_n1)

   # --- 3. Calculate Power ---
   cat("--- Calculating power for specified sample sizes... ---\n")
   z_alpha <- stats::qnorm(1 - alpha / 2)
   power_values <- sapply(sample_sizes, function(n_per_stratum) {
      n_strata <- length(unique(df[[strata_var]]))
      total_n <- n_per_stratum * n_strata
      se_final <- se_beta_n1 / sqrt(total_n)
      power <- stats::pnorm( (abs(beta_effect) / se_final) - z_alpha )
      return(power)
   })

   results_df <- data.frame(N_per_Stratum = sample_sizes, Power = power_values)

   # --- 4. Plot and Return ---
   p <- ggplot2::ggplot(results_df, ggplot2::aes(x = N_per_Stratum, y = Power)) +
      ggplot2::geom_line(color = "#E69F00", linewidth = 1) +
      ggplot2::geom_point(color = "#E69F00", size = 3) +
      ggplot2::labs(
         title = "Analytic Power Curve: Multiplicative Stratified RMST Model",
         subtitle = "Approximate method based on Wang et al. (2019).",
         x = "Sample Size Per Stratum", y = "Estimated Power"
      ) +
      ggplot2::ylim(0, 1) + ggplot2::theme_minimal()

   return(list(results_data = results_df, results_plot = p))
}

# Sample Size Search ------------------------------------------------------

#' @title Find Sample Size for a Multiplicative Stratified RMST Model (Analytic)
#' @description Calculates the required sample size for a target power using the analytic
#'   (approximate) method from Wang et al. (2019).
#'
#' @details
#' This function performs an iterative search for the sample size required to
#' achieve a specified `target_power`. It uses the same underlying theory and
#' log-linear approximation as `MS.power.analytical`. It performs a one-time
#' estimation of the log-RMST ratio and its asymptotic variance from the pilot data,
#' then uses these parameters in an analytic formula to efficiently find the
#' required sample size.
#'
#' @param pilot_data A `data.frame` containing pilot study data.
#' @param time_var A character string for the time-to-event variable.
#' @param status_var A character string for the event status variable (1=event, 0=censored).
#' @param arm_var A character string for the treatment arm variable (1=treatment, 0=control).
#' @param strata_var A character string for the stratification variable.
#' @param target_power A single numeric value for the desired power.
#' @param linear_terms An optional character vector of other covariate names.
#' @param tau The numeric value for the RMST truncation time.
#' @param alpha The significance level (Type I error rate).
#' @param n_start The starting sample size *per stratum* for the search.
#' @param n_step The increment in sample size at each step of the search.
#' @param max_n_per_arm The maximum sample size *per stratum* to search up to.
#'
#' @return A `list` containing:
#' \item{results_data}{A `data.frame` with the target power and required sample size.}
#' \item{results_plot}{A `ggplot` object visualizing the search path.}
#' \item{results_summary}{A `data.frame` summarizing the estimated log(RMST Ratio).}
#' @export
#' @importFrom survival Surv coxph basehaz
#' @importFrom stats lm as.formula complete.cases na.omit sd quantile pnorm qnorm coef vcov predict
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline geom_vline labs theme_minimal
#' @importFrom knitr kable
#' @examples
#' set.seed(456)
#' pilot_df_strat_effect <- data.frame(
#'  time = c(rexp(60, 0.15), rexp(60, 0.08)), # Effect
#'  status = rbinom(120, 1, 0.7),
#'  arm = rep(0:1, each = 60),
#'  region = factor(rep(c("A", "B"), each = 60))
#' )
#'
#' ss_results <- MS.ss.analytical(
#'  pilot_data = pilot_df_strat_effect,
#'  time_var = "time", status_var = "status", arm_var = "arm", strata_var = "region",
#'  target_power = 0.80, tau = 10,
#'  n_start = 100, n_step = 50
#' )
#' print(ss_results$results_data)
MS.ss.analytical <- function(pilot_data, time_var, status_var, arm_var, strata_var,
                             target_power, linear_terms = NULL, tau, alpha = 0.05,
                             n_start = 50, n_step = 25, max_n_per_arm = 2000) {

   # --- 1. Estimate Parameters and Variance from Pilot Data (One Time) ---
   cat("--- Estimating parameters from pilot data (log-linear approximation)... ---\n")

   core_vars <- c(time_var, status_var, arm_var, strata_var)
   covariates <- c(arm_var, linear_terms)
   all_vars <- c(core_vars, linear_terms)
   df <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]
   n_pilot <- nrow(df)
   df$Y_rmst <- pmin(df[[time_var]], tau)
   df$is_event <- df[[status_var]] == 1
   cens_formula <- stats::as.formula(paste0("survival::Surv(Y_rmst, is_event == 0) ~ ",
                                            paste(covariates, collapse = " + "),
                                            " + strata(", strata_var, ")"))
   fit_cens <- survival::coxph(cens_formula, data = df, ties = "breslow")
   bh_cens <- survival::basehaz(fit_cens, centered = FALSE)
   df$H_cens <- 0
   unique_strata_from_bh <- unique(bh_cens$strata)
   for(st in unique(df[[strata_var]])){
      st_label <- paste0(strata_var, "=", st)
      is_stratum <- df[[strata_var]] == st
      if (st_label %in% unique_strata_from_bh) {
         is_bh_stratum <- bh_cens$strata == st_label
         if(sum(is_bh_stratum) > 0){
            H_st <- stats::stepfun(bh_cens$time[is_bh_stratum], c(0, bh_cens$hazard[is_bh_stratum]))(df$Y_rmst[is_stratum])
            df$H_cens[is_stratum] <- H_st
         }
      }
   }
   df$weights <- exp(df$H_cens * exp(predict(fit_cens, newdata = df, type="lp", reference="zero")))
   finite_weights <- df$weights[is.finite(df$weights)]
   if (length(finite_weights) > 0) {
      weight_cap <- stats::quantile(finite_weights, probs = 0.99, na.rm = TRUE)
      df$weights[df$weights > weight_cap] <- weight_cap
   }
   df$weights[!is.finite(df$weights)] <- 0
   df$w_delta <- df$weights * df$is_event
   df$w_delta[is.na(df$w_delta)] <- 0

   fit_data <- df[df$w_delta > 0 & df$Y_rmst > 0, ]
   fit_weights <- fit_data$w_delta
   log_model_formula <- stats::as.formula(paste("log(Y_rmst) ~", paste(covariates, collapse=" + "), "+", strata_var))
   if(nrow(fit_data) < (length(covariates) + length(unique(fit_data[[strata_var]])))) {
      stop("Not enough data points to fit the approximation model after filtering.", call. = FALSE)
   }
   fit_log_lm <- lm(log_model_formula, data = fit_data, weights = fit_weights)

   beta_summary <- coef(summary(fit_log_lm))
   beta_effect <- beta_summary[arm_var, "Estimate"]
   V_hat <- vcov(fit_log_lm) * nrow(fit_data)
   var_beta_n1 <- V_hat[arm_var, arm_var]
   se_beta_n1 <- sqrt(var_beta_n1)

   # --- 2. Iterative Search for Sample Size ---
   cat("--- Searching for Sample Size (Method: Analytic/Approximation) ---\n")
   current_n <- n_start
   search_path <- list()
   final_n <- NA_integer_
   z_alpha <- stats::qnorm(1 - alpha / 2)
   n_strata <- length(unique(df[[strata_var]]))

   while (current_n <= max_n_per_arm) {
      total_n <- current_n * n_strata
      se_final <- se_beta_n1 / sqrt(total_n)
      calculated_power <- stats::pnorm((abs(beta_effect) / se_final) - z_alpha)
      if (!is.finite(calculated_power)) calculated_power <- 0

      search_path[[as.character(current_n)]] <- calculated_power
      cat(paste0("  N = ", current_n, "/stratum, Calculated Power = ", round(calculated_power, 3), "\n"))

      if (calculated_power >= target_power) {
         final_n <- current_n
         break
      }
      current_n <- current_n + n_step
   }

   if (is.na(final_n)) {
      warning(paste("Target power", target_power, "not achieved by max N of", max_n_per_arm), call. = FALSE)
      final_n <- max_n_per_arm
   }

   # --- 3. Finalize and Return Results ---
   results_summary <- data.frame(Statistic = "Assumed log(RMST Ratio) (from pilot)", Value = beta_effect)
   results_df <- data.frame(Target_Power = target_power, Required_N_per_Stratum = final_n)
   search_path_df <- data.frame(N_per_Stratum = as.integer(names(search_path)), Power = unlist(search_path))

   p <- ggplot2::ggplot(na.omit(search_path_df), ggplot2::aes(x = N_per_Stratum, y = Power)) +
      ggplot2::geom_line(color = "#009E73", linewidth = 1) +
      ggplot2::geom_point(color = "#009E73", size = 3) +
      ggplot2::geom_hline(yintercept = target_power, linetype = "dashed", color = "red") +
      ggplot2::geom_vline(xintercept = final_n, linetype = "dotted", color = "blue") +
      ggplot2::labs(
         title = "Analytic Sample Size Search: Multiplicative Stratified RMST Model",
         subtitle = "Power calculated from formula at each step (approximate method).",
         x = "Sample Size Per Stratum", y = "Calculated Power"
      ) + ggplot2::theme_minimal()

   cat("\n--- Calculation Summary ---\n")
   print(knitr::kable(results_df, caption = "Required Sample Size"))

   return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}

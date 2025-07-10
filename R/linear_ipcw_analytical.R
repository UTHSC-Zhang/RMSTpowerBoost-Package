# Power Calculation -------------------------------------------------------------------


#' @title Analyze Power for a Linear RMST Model (Analytic)
#' @description Performs power analysis using a direct formula based on the
#'   asymptotic variance estimator for the linear RMST model.
#'
#' @details
#' This function implements the analytic power calculation for the direct linear
#' regression model of the Restricted Mean Survival Time (RMST) proposed by Tian et al. (2014).
#' The core of the method is a weighted linear model of the form
#' \deqn{E[Y_i|Z_i] = \beta_0 + \beta_1 Z_i}
#' where \eqn{Y_i = \min(T_i, \L)} is the event time truncated at \eqn{\L}.
#'
#' To handle right-censoring, the method uses Inverse Probability of Censoring
#' Weighting (IPCW). The weight for an uncensored individual `i` is the inverse
#' of the probability of remaining uncensored until their event time, \eqn{w_i = \delta_i / \hat{G}(Y_i)},
#' where \eqn{\hat{G}(t) = P(C > t)} is the Kaplan-Meier estimate of the censoring distribution.
#'
#' Power is calculated analytically based on the asymptotic properties of the
#' coefficient estimators. The variance of the treatment effect estimator, \eqn{\hat{\beta}}, is derived from a
#' robust sandwich variance estimator of the form \eqn{A^{-1}B(A^{-1})'}. In this implementation,
#' `A` is the scaled information matrix \eqn{(X'WX)/n}, and `B` is the empirical second moment of the
#' influence functions, \eqn{(\sum \epsilon_i \epsilon_i')/n}, where \eqn{\epsilon_i} is the influence curve
#' for observation `i`. The resulting variance is used to calculate the standard error for a
#' given sample size, which in turn is used in the power formula.
#'
#' @param pilot_data A `data.frame` containing pilot study data.
#' @param time_var A character string specifying the name of the time-to-event variable.
#' @param status_var A character string specifying the name of the event status variable (1=event, 0=censored).
#' @param arm_var A character string specifying the name of the treatment arm variable (1=treatment, 0=control).
#' @param sample_sizes A numeric vector of sample sizes *per arm* to calculate power for.
#' @param linear_terms An optional character vector of other covariate names to include in the model.
#' @param L The numeric value for the RMST truncation time.
#' @param alpha The significance level for the power calculation (Type I error rate).
#'
#' @return A `list` containing:
#' \item{results_data}{A `data.frame` with the specified sample sizes and their corresponding calculated power.}
#' \item{results_plot}{A `ggplot` object visualizing the power curve.}
#' \item{results_summary}{A `data.frame` summarizing the treatment effect from the pilot data used for the calculation.}
#'
#' @importFrom survival Surv survfit
#' @importFrom stats lm as.formula complete.cases na.omit sd quantile pnorm qnorm model.matrix coef vcov predict
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline labs theme_minimal ylim
#' @importFrom knitr kable
#' @export
#' @examples
#' pilot_df <- data.frame(
#'   time = rexp(100, 0.1),
#'   status = rbinom(100, 1, 0.7),
#'   arm = rep(0:1, each = 50),
#'   age = rnorm(100, 55, 10)
#' )
#' power_results <- linear.power.analytical(
#'   pilot_data = pilot_df,
#'   time_var = "time",
#'   status_var = "status",
#'   arm_var = "arm",
#'   linear_terms = "age",
#'   sample_sizes = c(100, 200, 300),
#'   L = 10
#' )
#' print(power_results$results_data)
#' print(power_results$results_plot)
linear.power.analytical <- function(pilot_data, time_var, status_var, arm_var,
                                    sample_sizes, linear_terms = NULL, L, alpha = 0.05) {

   # --- 1. Estimate Nuisance Parameters from Pilot Data ---
   cat("--- Estimating parameters from pilot data for analytic calculation... ---\n")

   core_vars <- c(time_var, status_var, arm_var)
   all_vars <- c(core_vars, linear_terms)
   df <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]
   n_pilot <- nrow(df)

   # Define model formulas
   factor_arm_str <- paste0("factor(", arm_var, ")")
   model_rhs <- paste(c(factor_arm_str, linear_terms), collapse = " + ")
   model_formula <- stats::as.formula(paste("Y_rmst ~", model_rhs))
   message("Model: Y_rmst ~ ", model_rhs)

   # Prepare data for IPCW
   df$Y_rmst <- pmin(df[[time_var]], L)
   df$is_censored <- df[[status_var]] == 0
   df$is_event <- df[[status_var]] == 1

   # Fit censoring model (Kaplan-Meier for G(t))
   cens_fit <- survival::survfit(survival::Surv(Y_rmst, is_censored) ~ 1, data = df)
   cens_surv_prob <- stats::stepfun(cens_fit$time, c(1, cens_fit$surv))(df$Y_rmst)
   df$weights <- df$is_event / cens_surv_prob

   # Stabilize weights
   finite_weights <- df$weights[is.finite(df$weights) & df$weights > 0]
   if (length(finite_weights) > 0) {
      weight_cap <- stats::quantile(finite_weights, probs = 0.99, na.rm = TRUE)
      df$weights[df$weights > weight_cap] <- weight_cap
   }
   df$weights[!is.finite(df$weights)] <- 0

   # Filter for model fitting
   fit_data <- df[df$weights > 0, ]
   fit_weights <- fit_data$weights

   if (length(unique(fit_data[[arm_var]])) < 2) {
      stop("Pilot data contains events in only one arm after filtering.", call. = FALSE)
   }

   # Fit the primary linear model
   fit_lm <- stats::lm(model_formula, data = fit_data, weights = fit_weights)
   beta_hat <- stats::coef(fit_lm)

   arm_pattern <- paste0("^factor\\(", arm_var, "\\)1$")
   arm_coeff_name <- names(beta_hat)[grep(arm_pattern, names(beta_hat))]
   if (length(arm_coeff_name) == 0) stop("Could not find treatment effect coefficient.")
   beta_effect <- beta_hat[arm_coeff_name]

   # --- 2. Calculate Asymptotic Sandwich Variance Components ---
   cat("--- Calculating asymptotic variance... ---\n")
   X <- stats::model.matrix(model_formula, data = df)

   A_hat <- crossprod(X * sqrt(df$weights), X * sqrt(df$weights)) / n_pilot

   A_hat_inv <- tryCatch({
      solve(A_hat)
   }, error = function(e) {
      stop("The design matrix (A_hat) is singular. This can happen with small pilot datasets or perfect separation.", call. = FALSE)
   })

   df$predicted_rmst <- predict(fit_lm, newdata = df)
   residuals <- df$Y_rmst - df$predicted_rmst

   epsilon <- X * residuals * df$weights
   epsilon[is.na(epsilon)] <- 0
   B_hat <- crossprod(epsilon) / n_pilot

   V_hat_n <- A_hat_inv %*% B_hat %*% t(A_hat_inv)

   # Scale to get the variance for n=1
   var_beta_n1 <- V_hat_n[arm_coeff_name, arm_coeff_name]
   se_beta_n1 <- sqrt(var_beta_n1)

   # --- 4. Calculate Power for Each Sample Size ---
   cat("--- Calculating power for specified sample sizes... ---\n")
   z_alpha <- stats::qnorm(1 - alpha / 2)
   power_values <- sapply(sample_sizes, function(n_per_arm) {
      total_n <- n_per_arm * 2
      se_final <- se_beta_n1 / sqrt(total_n)
      power <- stats::pnorm( (abs(beta_effect) / se_final) - z_alpha )
      return(power)
   })

   results_df <- data.frame(N_per_Arm = sample_sizes, Power = power_values)

   # --- CORRECTED: Add the summary object to the return list ---
   results_summary <- data.frame(
      Statistic = "Assumed RMST Difference (from pilot)",
      Value = beta_effect
   )

   # --- 5. Create Plot and Return ---
   p <- ggplot2::ggplot(results_df, ggplot2::aes(x = N_per_Arm, y = Power)) +
      ggplot2::geom_line(color = "#D55E00", linewidth = 1) +
      ggplot2::geom_point(color = "#D55E00", size = 3) +
      ggplot2::labs(
         title = "Analytic Power Curve: Linear IPCW RMST Model",
         subtitle = "Based on the asymptotic variance from Tian et al. (2014).",
         x = "Sample Size Per Arm", y = "Estimated Power"
      ) +
      ggplot2::ylim(0, 1) + ggplot2::theme_minimal()

   return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}

# Sample Size Search ------------------------------------------------------


#' @title Find Sample Size for a Linear RMST Model (Analytic)
#' @description Calculates the required sample size for a target power using an
#'   analytic formula based on the methods of Tian et al. (2014).
#'
#' @details
#' This function performs an iterative search to find the sample size needed to
#' achieve a specified `target_power`. It uses the same underlying theory as
#' `linear.power.analytical`. First, it estimates the treatment effect size and its
#' asymptotic variance from the pilot data. Then, it iteratively calculates the
#' power for increasing sample sizes using the analytic formula until the
#' target power is achieved.
#'
#' @param pilot_data A `data.frame` containing pilot study data.
#' @param time_var A character string specifying the name of the time-to-event variable.
#' @param status_var A character string specifying the name of the event status variable (1=event, 0=censored).
#' @param arm_var A character string specifying the name of the treatment arm variable (1=treatment, 0=control).
#' @param target_power A single numeric value for the desired power (e.g., 0.80 or 0.90).
#' @param linear_terms An optional character vector of other covariate names to include in the model.
#' @param L The numeric value for the RMST truncation time.
#' @param alpha The significance level (Type I error rate).
#' @param n_start The starting sample size *per arm* for the search.
#' @param n_step The increment in sample size at each step of the search.
#' @param max_n_per_arm The maximum sample size *per arm* to search up to.
#'
#' @return A `list` containing:
#' \item{results_data}{A `data.frame` with the target power and the required sample size per arm.}
#' \item{results_plot}{A `ggplot` object visualizing the sample size search path.}
#' \item{results_summary}{A `data.frame` summarizing the treatment effect from the pilot data used for the calculation.}
#'
#' @importFrom survival Surv survfit
#' @importFrom stats lm as.formula complete.cases na.omit sd quantile pnorm qnorm model.matrix coef vcov predict
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline geom_vline labs theme_minimal
#' @importFrom knitr kable
#' @export
#' @examples
#' pilot_df <- data.frame(
#'   time = c(rexp(50, 0.1), rexp(50, 0.07)), # Introduce an effect
#'   status = rbinom(100, 1, 0.8),
#'   arm = rep(0:1, each = 50),
#'   age = rnorm(100, 55, 10)
#' )
#' ss_results <- linear.ss.analytical(
#'   pilot_data = pilot_df,
#'   time_var = "time",
#'   status_var = "status",
#'   arm_var = "arm",
#'   target_power = 0.80,
#'   L = 10
#' )
#' print(ss_results$results_data)
#' print(ss_results$results_plot)
linear.ss.analytical <- function(pilot_data, time_var, status_var, arm_var,
                                 target_power, linear_terms = NULL, L, alpha = 0.05,
                                 n_start = 50, n_step = 25, max_n_per_arm = 2000) {

   # --- 1. Estimate Parameters and Variance from Pilot Data (One Time) ---
   cat("--- Estimating parameters from pilot data for analytic search... ---\n")
   # This part is identical to the power function above
   core_vars <- c(time_var, status_var, arm_var)
   all_vars <- c(core_vars, linear_terms)
   df <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]
   n_pilot <- nrow(df)

   factor_arm_str <- paste0("factor(", arm_var, ")")
   model_rhs <- paste(c(factor_arm_str, linear_terms), collapse = " + ")
   model_formula <- stats::as.formula(paste("Y_rmst ~", model_rhs))
   message("Model: Y_rmst ~ ", model_rhs)

   df$Y_rmst <- pmin(df[[time_var]], L)
   df$is_censored <- df[[status_var]] == 0
   df$is_event <- df[[status_var]] == 1

   cens_fit <- survival::survfit(survival::Surv(Y_rmst, is_censored) ~ 1, data = df)
   cens_surv_prob <- stats::stepfun(cens_fit$time, c(1, cens_fit$surv))(df$Y_rmst)
   df$weights <- df$is_event / cens_surv_prob
   finite_weights <- df$weights[is.finite(df$weights) & df$weights > 0]
   if (length(finite_weights) > 0) {
      weight_cap <- stats::quantile(finite_weights, probs = 0.99, na.rm = TRUE)
      df$weights[df$weights > weight_cap] <- weight_cap
   }
   df$weights[!is.finite(df$weights)] <- 0

   fit_data <- df[df$weights > 0, ]
   fit_weights <- fit_data$weights
   if (length(unique(fit_data[[arm_var]])) < 2) stop("Pilot data events in only one arm.")

   fit_lm <- stats::lm(model_formula, data = fit_data, weights = fit_weights)
   beta_hat <- stats::coef(fit_lm)
   arm_pattern <- paste0("^factor\\(", arm_var, "\\)1$")
   arm_coeff_name <- names(beta_hat)[grep(arm_pattern, names(beta_hat))]
   if (length(arm_coeff_name) == 0) stop("Could not find treatment effect coefficient.")
   beta_effect <- beta_hat[arm_coeff_name]

   X <- stats::model.matrix(model_formula, data = df)
   A_hat <- crossprod(X * sqrt(df$weights), X * sqrt(df$weights)) / n_pilot

   A_hat_inv <- tryCatch({
      solve(A_hat)
   }, error = function(e) {
      stop("The design matrix (A_hat) is singular. This can happen with small pilot datasets or perfect separation.", call. = FALSE)
   })

   df$predicted_rmst <- predict(fit_lm, newdata = df)
   residuals <- df$Y_rmst - df$predicted_rmst
   epsilon <- X * residuals * df$weights
   epsilon[is.na(epsilon)] <- 0
   B_hat <- crossprod(epsilon) / n_pilot

   V_hat_n <- A_hat_inv %*% B_hat %*% t(A_hat_inv)
   var_beta_n1 <- V_hat_n[arm_coeff_name, arm_coeff_name]
   se_beta_n1 <- sqrt(var_beta_n1)

   # --- 2. Iterative Search for Sample Size using Analytic Formula ---
   cat("--- Searching for Sample Size (Method: Analytic) ---\n")
   current_n <- n_start
   search_path <- list()
   final_n <- NA_integer_
   z_alpha <- stats::qnorm(1 - alpha / 2)

   while (current_n <= max_n_per_arm) {
      total_n <- current_n * 2
      se_final <- se_beta_n1 / sqrt(total_n)
      calculated_power <- stats::pnorm((abs(beta_effect) / se_final) - z_alpha)
      if (!is.finite(calculated_power)) calculated_power <- 0

      search_path[[as.character(current_n)]] <- calculated_power
      cat(paste0("  N = ", current_n, "/arm, Calculated Power = ", round(calculated_power, 3), "\n"))

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
   results_summary <- data.frame(
      Statistic = "Assumed RMST Difference (from pilot)",
      Value = beta_effect
   )
   results_df <- data.frame(Target_Power = target_power, Required_N_per_Arm = final_n)
   search_path_df <- data.frame(N_per_Arm = as.integer(names(search_path)), Power = unlist(search_path))

   p <- ggplot2::ggplot(na.omit(search_path_df), ggplot2::aes(x = N_per_Arm, y = Power)) +
      ggplot2::geom_line(color = "#009E73", linewidth = 1) +
      ggplot2::geom_point(color = "#009E73", size = 3) +
      ggplot2::geom_hline(yintercept = target_power, linetype = "dashed", color = "red") +
      ggplot2::geom_vline(xintercept = final_n, linetype = "dotted", color = "blue") +
      ggplot2::labs(
         title = "Analytic Sample Size Search: Linear IPCW RMST Model",
         subtitle = "Power calculated from formula at each step.",
         x = "Sample Size Per Arm", y = "Calculated Power"
      ) + ggplot2::theme_minimal()

   cat("\n--- Calculation Summary ---\n")
   print(knitr::kable(results_df, caption = "Required Sample Size"))

   return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}


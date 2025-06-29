#' @title Analyze Power for a Dependent Censoring RMST Model via Analytic Formula
#' @description Performs power analysis using a direct formula-based method based
#'   on the asymptotic variance of the estimator from Wang & Schaubel (2018).
#'   This version uses a simplified ("ASE1") sandwich variance estimator.
#'
#' @note This function does not use bootstrapping and is therefore much faster
#'   but requires more complex calculations based on statistical theory.
#'
#' @param pilot_data A data.frame with pilot study data.
#' @param time_var,status_var,arm_var,dep_cens_status_var Strings for column names.
#' @param sample_sizes A numeric vector of sample sizes per arm to calculate power for.
#' @param linear_terms Optional vector of covariates for the model.
#' @param tau The truncation time for RMST.
#' @param alpha The significance level. Default is 0.05.
#'
#' @return A list containing `results_data` (a data.frame of sample sizes and
#'   corresponding powers) and `results_plot` (a ggplot object).
#'
#' @references Wang, X., & Schaubel, D. E. (2018). Modeling Restricted Mean Survival
#'   Time Under General Censoring Mechanisms. *Lifetime Data Analysis*, 24, 176-199.
#'
design_rmst_dc_power<- function(pilot_data,
                                          time_var,
                                          status_var,
                                          arm_var,
                                          dep_cens_status_var,
                                          sample_sizes,
                                          linear_terms = NULL,
                                          tau,
                                          alpha = 0.05) {

   # --- 1. Estimate Nuisance Parameters from Pilot Data ---
   cat("--- Estimating parameters from pilot data... ---\n")

   core_vars <- c(time_var, status_var, arm_var, dep_cens_status_var)
   if (is.null(linear_terms)) {
      linear_terms <- setdiff(names(pilot_data), core_vars)
      if (length(linear_terms) > 0) message(paste("Using unspecified columns as linear terms:", paste(linear_terms, collapse = ", ")))
   }

   all_vars <- c(core_vars, linear_terms)
   df <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]

   # Define model formulas
   factor_arm_str <- paste0("factor(", arm_var, ")")
   model_rhs <- paste(c(factor_arm_str, linear_terms), collapse = " + ")

   cens_formula <- stats::as.formula(paste("survival::Surv(", time_var, ", cens_ind) ~", model_rhs))
   dep_cens_formula <- stats::as.formula(paste("survival::Surv(", time_var, ", dep_cens_ind) ~", model_rhs))
   model_formula <- stats::as.formula(paste("Y_rmst ~", model_rhs))

   message("Model: Y_rmst ~ ", model_rhs)

   # Define censoring indicators
   df$cens_ind <- df[[status_var]] == 0 & df[[dep_cens_status_var]] == 0
   df$dep_cens_ind <- df[[dep_cens_status_var]] == 1
   df$Y_rmst <- pmin(df[[time_var]], tau)

   # Fit censoring models
   fit_cens <- tryCatch(survival::coxph(cens_formula, data = df, ties = "breslow"), error = function(e) NULL)
   fit_dep_cens <- tryCatch(survival::coxph(dep_cens_formula, data = df, ties = "breslow"), error = function(e) NULL)

   if(is.null(fit_cens) || is.null(fit_dep_cens)) {
      stop("One or both of the censoring models failed to fit on the pilot data.", call. = FALSE)
   }

   # Calculate IPCW weights
   bh_cens <- survival::basehaz(fit_cens, centered = FALSE)
   H_cens <- stats::stepfun(bh_cens$time, c(0, bh_cens$hazard))(df$Y_rmst) * exp(predict(fit_cens, type = "lp"))

   bh_dep_cens <- survival::basehaz(fit_dep_cens, centered = FALSE)
   H_dep_cens <- stats::stepfun(bh_dep_cens$time, c(0, bh_dep_cens$hazard))(df$Y_rmst) * exp(predict(fit_dep_cens, type = "lp"))

   df$weights <- exp(H_cens + H_dep_cens)

   # Stabilize weights
   finite_weights <- df$weights[is.finite(df$weights)]
   if (length(finite_weights) > 0) {
      weight_cap <- stats::quantile(finite_weights, probs = 0.99, na.rm = TRUE)
      df$weights[df$weights > weight_cap] <- weight_cap
   }
   df$weights[!is.finite(df$weights)] <- NA

   # Filter data for final model
   event_indices <- which(df[[status_var]] == 1 & !is.na(df$weights))
   fit_data <- df[event_indices, ]
   fit_weights <- df$weights[event_indices]

   if (length(unique(fit_data[[arm_var]])) < 2) {
      stop("Pilot data contains primary events in only one arm after filtering. Cannot estimate treatment effect.", call. = FALSE)
   }

   fit_lm <- stats::lm(model_formula, data = fit_data, weights = fit_weights)
   beta_hat <- stats::coef(fit_lm)

   arm_pattern <- paste0("^factor\\(", arm_var, "\\)")
   arm_coeff_name <- names(beta_hat)[grep(arm_pattern, names(beta_hat))]

   if (length(arm_coeff_name) == 0) {
      stop("Could not find the treatment effect coefficient after fitting the model.", call. = FALSE)
   }
   beta_effect <- beta_hat[arm_coeff_name]

   # --- 2. Calculate Asymptotic Variance Components ---
   cat("--- Calculating asymptotic variance... ---\n")

   X <- stats::model.matrix(model_formula, data = df)
   A_hat <- crossprod(X) / nrow(df)

   df$predicted_rmst <- predict(fit_lm, newdata = df)
   epsilon <- X * (df$Y_rmst - df$predicted_rmst) * df$weights * (df[[status_var]] == 1)
   epsilon[is.na(epsilon)] <- 0
   B_star_hat <- crossprod(epsilon) / nrow(df)

   # --- 3. Assemble the Sandwich Variance Estimator ---
   A_hat_inv <- solve(A_hat)
   V_hat <- A_hat_inv %*% B_star_hat %*% t(A_hat_inv)

   var_beta_effect <- diag(V_hat)[arm_coeff_name]
   sigma_beta_effect <- sqrt(var_beta_effect)

   # --- 4. Calculate Power for Each Sample Size ---
   cat("--- Calculating power for specified sample sizes... ---\n")

   total_n <- sample_sizes * 2
   z_alpha <- stats::qnorm(1 - alpha / 2)
   power_values <- stats::pnorm( (abs(beta_effect) * sqrt(total_n) / sigma_beta_effect) - z_alpha )

   results_df <- data.frame(N_per_Arm = sample_sizes, Power = power_values)

   # --- 5. Create Plot and Return ---
   p <- ggplot2::ggplot(results_df, ggplot2::aes(x = N_per_Arm, y = Power)) +
      ggplot2::geom_line(color = "#0072B2", linewidth = 1) +
      ggplot2::geom_point(color = "#0072B2", size = 3) +
      ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
      ggplot2::labs(
         title = "Analytic Power Curve: Dependent Censoring RMST Model",
         subtitle = "Based on the asymptotic variance of the estimator.",
         x = "Sample Size Per Arm", y = "Estimated Power"
      ) +
      ggplot2::ylim(0, 1) + ggplot2::theme_minimal()

   return(list(results_data = results_df, results_plot = p))
}

#' @title Find Sample Size for a Dependent Censoring RMST Model via Analytic Search
#' @description Calculates the required sample size by performing an iterative search
#'   that uses a direct formula for power at each step.
#'
#' @note This function does not use bootstrapping. It calculates power analytically
#'   at each step of the search to produce a power curve.
#'
#' @param pilot_data A data.frame with pilot study data.
#' @param time_var,status_var,arm_var,dep_cens_status_var Strings for column names.
#' @param target_power A single numeric value for the target power (e.g., 0.80).
#' @param linear_terms Optional vector of covariates for the model.
#' @param tau The truncation time for RMST.
#' @param alpha The significance level. Default is 0.05.
#' @param n_start The starting sample size per arm for the search.
#' @param n_step The increment in sample size at each step of the search.
#' @param max_n_per_arm The maximum sample size per arm to search up to.
#'
#' @return A list containing `results_data`, `results_plot`, and `results_summary`.
#'
design_rmst_dc_ss<- function(pilot_data,
                                       time_var,
                                       status_var,
                                       arm_var,
                                       dep_cens_status_var,
                                       target_power,
                                       linear_terms = NULL,
                                       tau,
                                       alpha = 0.05,
                                       n_start = 50,
                                       n_step = 25,
                                       max_n_per_arm = 2000) {

   # --- 1. Estimate Parameters and Variance from Pilot Data (One Time) ---
   cat("--- Estimating parameters from pilot data for analytic calculation... ---\n")

   core_vars <- c(time_var, status_var, arm_var, dep_cens_status_var)
   if (is.null(linear_terms)) {
      linear_terms <- setdiff(names(pilot_data), core_vars)
      if (length(linear_terms) > 0) message(paste("Using unspecified columns as linear terms:", paste(linear_terms, collapse = ", ")))
   }

   all_vars <- c(core_vars, linear_terms)
   df <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]

   factor_arm_str <- paste0("factor(", arm_var, ")")
   model_rhs <- paste(c(factor_arm_str, linear_terms), collapse = " + ")

   cens_formula <- stats::as.formula(paste("survival::Surv(", time_var, ", cens_ind) ~", model_rhs))
   dep_cens_formula <- stats::as.formula(paste("survival::Surv(", time_var, ", dep_cens_ind) ~", model_rhs))
   model_formula <- stats::as.formula(paste("Y_rmst ~", model_rhs))

   message("Model: Y_rmst ~ ", model_rhs)

   df$cens_ind <- df[[status_var]] == 0 & df[[dep_cens_status_var]] == 0
   df$dep_cens_ind <- df[[dep_cens_status_var]] == 1
   df$Y_rmst <- pmin(df[[time_var]], tau)

   fit_cens <- tryCatch(survival::coxph(cens_formula, data = df, ties = "breslow"), error = function(e) NULL)
   fit_dep_cens <- tryCatch(survival::coxph(dep_cens_formula, data = df, ties = "breslow"), error = function(e) NULL)

   if(is.null(fit_cens) || is.null(fit_dep_cens)) {
      stop("One or both of the censoring models failed to fit on the pilot data.", call. = FALSE)
   }

   bh_cens <- survival::basehaz(fit_cens, centered = FALSE)
   H_cens <- stats::stepfun(bh_cens$time, c(0, bh_cens$hazard))(df$Y_rmst) * exp(predict(fit_cens, type = "lp"))
   bh_dep_cens <- survival::basehaz(fit_dep_cens, centered = FALSE)
   H_dep_cens <- stats::stepfun(bh_dep_cens$time, c(0, bh_dep_cens$hazard))(df$Y_rmst) * exp(predict(fit_dep_cens, type = "lp"))
   df$weights <- exp(H_cens + H_dep_cens)

   finite_weights <- df$weights[is.finite(df$weights)]
   if (length(finite_weights) > 0) {
      weight_cap <- stats::quantile(finite_weights, probs = 0.99, na.rm = TRUE)
      df$weights[df$weights > weight_cap] <- weight_cap
   }
   df$weights[!is.finite(df$weights)] <- NA

   event_indices <- which(df[[status_var]] == 1 & !is.na(df$weights))
   fit_data <- df[event_indices, ]
   fit_weights <- df$weights[event_indices]

   if (length(unique(fit_data[[arm_var]])) < 2) {
      stop("Pilot data contains primary events in only one arm after filtering. Cannot estimate treatment effect.", call. = FALSE)
   }

   fit_lm <- stats::lm(model_formula, data = fit_data, weights = fit_weights)
   beta_hat <- stats::coef(fit_lm)

   arm_pattern <- paste0("^factor\\(", arm_var, "\\)")
   arm_coeff_name <- names(beta_hat)[grep(arm_pattern, names(beta_hat))]
   if (length(arm_coeff_name) == 0) {
      stop("Could not find the treatment effect coefficient after fitting the model.", call. = FALSE)
   }
   beta_effect <- beta_hat[arm_coeff_name]

   X <- stats::model.matrix(model_formula, data = df)
   A_hat <- crossprod(X) / nrow(df)
   df$predicted_rmst <- predict(fit_lm, newdata = df)
   epsilon <- X * (df$Y_rmst - df$predicted_rmst) * df$weights * (df[[status_var]] == 1)
   epsilon[is.na(epsilon)] <- 0
   B_star_hat <- crossprod(epsilon) / nrow(df)
   A_hat_inv <- solve(A_hat)
   V_hat <- A_hat_inv %*% B_star_hat %*% t(A_hat_inv)
   var_beta_effect <- diag(V_hat)[arm_coeff_name]
   sigma_beta_effect <- sqrt(var_beta_effect)

   # --- 2. Iterative Search for Sample Size using Analytic Formula ---
   cat("--- Searching for Sample Size (Method: Analytic) ---\n")
   current_n <- n_start
   search_path <- list()
   final_n <- NA_integer_
   z_alpha <- stats::qnorm(1 - alpha / 2)

   while (current_n <= max_n_per_arm) {
      total_n <- current_n * 2
      calculated_power <- stats::pnorm((abs(beta_effect) * sqrt(total_n) / sigma_beta_effect) - z_alpha)
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

   # Create data for the plot from the search path
   search_path_df <- data.frame(
      N_per_Arm = as.integer(names(search_path)),
      Power = as.numeric(unlist(search_path))
   )

   # Generate the plot
   p <- ggplot2::ggplot(na.omit(search_path_df), ggplot2::aes(x = N_per_Arm, y = Power)) +
      ggplot2::geom_line(color = "#009E73", linewidth = 1) +
      ggplot2::geom_point(color = "#009E73", size = 3) +
      ggplot2::geom_hline(yintercept = target_power, linetype = "dashed", color = "red") +
      ggplot2::geom_vline(xintercept = final_n, linetype = "dotted", color = "blue") +
      ggplot2::labs(
         title = "Analytic Sample Size Search Path",
         subtitle = "Power is calculated directly from formula at each step.",
         x = "Sample Size Per Arm", y = "Calculated Power"
      ) +
      ggplot2::theme_minimal()

   cat("\n--- Calculation Summary ---\n")
   # Print the results table
   print(knitr::kable(results_df, caption = "Required Sample Size"))

   return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}

# Power Calculation -------------------------------------------------------

#' @title Analyze Power for RMST Model with Dependent Censoring (Analytic)
#' @description Performs power analysis for an RMST model with multiple censoring causes
#'   using a direct formula-based method.
#'
#' @details
#' This function calculates power based on a linear model for RMST in the presence
#' of both independent and dependent censoring (or competing risks). The method uses
#' cause-specific IPCW to account for the censoring events.
#'
#' Specifically, it fits a separate Cox proportional hazards model for each censoring
#' cause to estimate the cause-specific hazards. The final IPCW weight for each
#' subject is calculated by combining the cumulative hazards from all censoring causes:
#' \deqn{W_i = \exp(\sum_{k=1}^{K} \hat{\Lambda}_{k}(Y_i))}
#' where \eqn{\hat{\Lambda}_{k}} is the estimated cumulative hazard for censoring cause `k`.
#'
#' The power is then derived from the asymptotic sandwich variance of the treatment
#' effect in the resulting weighted linear model. This implementation uses a robust
#' variance estimator that provides a good approximation but does not account for the
#' variability from estimating the IPCW weights themselves.
#'
#' @param pilot_data A `data.frame` containing pilot study data.
#' @param time_var A character string for the time-to-event variable.
#' @param status_var A character string for the primary event status (1=event, 0=otherwise).
#' @param arm_var A character string for the treatment arm variable (1=treatment, 0=control).
#' @param dep_cens_status_var A character string for the dependent censoring status (1=dependent event, 0=otherwise).
#' @param sample_sizes A numeric vector of sample sizes *per arm* to calculate power for.
#' @param linear_terms An optional character vector of other covariate names.
#' @param L The numeric value for the RMST truncation time.
#' @param alpha The significance level (Type I error rate).
#'
#' @return A `list` containing:
#' \item{results_data}{A `data.frame` with sample sizes and corresponding powers.}
#' \item{results_plot}{A `ggplot` object visualizing the power curve.}
#' @export
#' @importFrom stats pnorm qnorm
#' @importFrom survival Surv coxph basehaz
#' @importFrom stats lm as.formula complete.cases na.omit sd quantile pnorm qnorm model.matrix coef vcov predict
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline labs theme_minimal ylim
#' @examples
#' # Generate sample pilot data
#' set.seed(123)
#' n_pilot <- 150
#' pilot_df <- data.frame(
#'   time = rexp(n_pilot, rate = 0.1),
#'   arm = rep(0:1, each = n_pilot / 2),
#'   age = rnorm(n_pilot, mean = 60, sd = 10)
#' )
#' # Introduce a treatment effect
#' pilot_df$time[pilot_df$arm == 1] <- pilot_df$time[pilot_df$arm == 1] * 1.2
#'
#' # Create competing event indicators
#' # Assume 70% primary event, 15% dependent censoring, 15% independent censoring
#' event_type <- sample(0:2, n_pilot, replace = TRUE, prob = c(0.7, 0.15, 0.15))
#' pilot_df$status <- ifelse(event_type == 0, 1, 0)
#' pilot_df$dep_cens_status <- ifelse(event_type == 1, 1, 0)
#' pilot_df$time[event_type != 0] <- pilot_df$time[event_type != 0] * 0.8
#'
#' # Run the power analysis
#' dc_power_results <- DC.power.analytical(
#'   pilot_data = pilot_df,
#'   time_var = "time",
#'   status_var = "status",
#'   arm_var = "arm",
#'   dep_cens_status_var = "dep_cens_status",
#'   sample_sizes = c(200, 300, 400),
#'   linear_terms = "age",
#'   L = 20,
#'   alpha = 0.05
#' )
#' print(dc_power_results$results_data)
#' print(dc_power_results$results_plot)
DC.power.analytical <- function(pilot_data,
                                time_var,
                                status_var,
                                arm_var,
                                dep_cens_status_var,
                                sample_sizes,
                                linear_terms = NULL,
                                L,
                                alpha = 0.05) {

   # --- 1. Estimate Nuisance Parameters from Pilot Data ---
   cat("--- Estimating parameters from pilot data... ---\n")

   core_vars <- c(time_var, status_var, arm_var, dep_cens_status_var)
   covariates <- c(arm_var, linear_terms)
   all_vars <- c(core_vars, covariates)
   df <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]
   n_pilot <- nrow(df)

   # Define model formulas
   model_rhs <- paste(covariates, collapse = " + ")
   cens_formula <- stats::as.formula(paste("survival::Surv(", time_var, ", cens_ind) ~", model_rhs))
   dep_cens_formula <- stats::as.formula(paste("survival::Surv(", time_var, ", dep_cens_ind) ~", model_rhs))
   model_formula <- stats::as.formula(paste("Y_rmst ~", model_rhs))

   message("Model: Y_rmst ~ ", model_rhs)

   # Define censoring indicators and truncated time
   df$cens_ind <- df[[status_var]] == 0 & df[[dep_cens_status_var]] == 0
   df$dep_cens_ind <- df[[dep_cens_status_var]] == 1
   df$Y_rmst <- pmin(df[[time_var]], L)

   # Fit censoring models
   fit_cens <- tryCatch(survival::coxph(cens_formula, data = df, ties = "breslow"), error = function(e) NULL)
   fit_dep_cens <- tryCatch(survival::coxph(dep_cens_formula, data = df, ties = "breslow"), error = function(e) NULL)

   if(is.null(fit_cens) || is.null(fit_dep_cens)) {
      stop("One or both of the censoring models failed to fit on the pilot data.", call. = FALSE)
   }

   # Calculate IPCW weights
   bh_cens <- survival::basehaz(fit_cens, centered = FALSE)
   H_cens <- stats::stepfun(bh_cens$time, c(0, bh_cens$hazard))(df$Y_rmst) * exp(predict(fit_cens, newdata = df, type = "lp"))

   bh_dep_cens <- survival::basehaz(fit_dep_cens, centered = FALSE)
   H_dep_cens <- stats::stepfun(bh_dep_cens$time, c(0, bh_dep_cens$hazard))(df$Y_rmst) * exp(predict(fit_dep_cens, newdata = df, type = "lp"))

   df$weights <- exp(H_cens + H_dep_cens)
   df$weights[df[[status_var]] != 1] <- 0 # Weights are only applied to uncensored subjects

   # Stabilize weights
   finite_weights <- df$weights[is.finite(df$weights) & df$weights > 0]
   if (length(finite_weights) > 0) {
      weight_cap <- stats::quantile(finite_weights, probs = 0.99, na.rm = TRUE)
      df$weights[df$weights > weight_cap] <- weight_cap
   }
   df$weights[!is.finite(df$weights)] <- 0

   # Filter data for final model
   fit_data <- df[df$weights > 0, ]
   fit_weights <- fit_data$weights

   if (length(unique(fit_data[[arm_var]])) < 2) {
      stop("Pilot data contains primary events in only one arm after filtering. Cannot estimate treatment effect.", call. = FALSE)
   }

   # Fit the primary linear model
   fit_lm <- stats::lm(model_formula, data = fit_data, weights = fit_weights)
   beta_hat <- stats::coef(fit_lm)
   arm_coeff_name <- arm_var # Since arm_var is a covariate directly
   beta_effect <- beta_hat[arm_coeff_name]

   # --- 2. Calculate Asymptotic Sandwich Variance Components ---
   cat("--- Calculating asymptotic variance... ---\n")
   X <- stats::model.matrix(model_formula, data = df)

   # A matrix (derivative of estimating function)
   # For weighted least squares, it's X'WX / n, but here we use X'X/n as per some formulations
   A_hat <- crossprod(X) / n_pilot

   # B matrix (outer product of influence functions)
   # Note: This is an approximation that ignores the influence from estimating the censoring models.
   df$predicted_rmst <- predict(fit_lm, newdata = df)
   residuals <- df$Y_rmst - df$predicted_rmst
   epsilon <- X * residuals * df$weights
   epsilon[is.na(epsilon)] <- 0
   B_hat <- crossprod(epsilon) / n_pilot

   # --- 3. Assemble the Sandwich Variance Estimator ---
   # This gives Var(sqrt(n) * (beta_hat - beta))
   A_hat_inv <- solve(A_hat)
   V_hat_n <- A_hat_inv %*% B_hat %*% t(A_hat_inv)

   # This is the variance of beta_hat for the pilot sample size n_pilot
   var_beta_pilot <- V_hat_n[arm_coeff_name, arm_coeff_name] / n_pilot
   # Scale to get the variance for a sample size of 1
   var_beta_n1 <- var_beta_pilot * n_pilot
   se_beta_n1 <- sqrt(var_beta_n1)

   # --- 4. Calculate Power for Each Sample Size ---
   cat("--- Calculating power for specified sample sizes... ---\n")
   z_alpha <- stats::qnorm(1 - alpha / 2)
   power_values <- sapply(sample_sizes, function(n_per_arm) {
      total_n <- n_per_arm * 2 # Assuming two arms
      se_final <- se_beta_n1 / sqrt(total_n)
      power <- stats::pnorm( (abs(beta_effect) / se_final) - z_alpha )
      return(power)
   })

   results_df <- data.frame(N_per_Arm = sample_sizes, Power = power_values)

   # --- 5. Create Plot and Return ---
   p <- ggplot2::ggplot(results_df, ggplot2::aes(x = N_per_Arm, y = Power)) +
      ggplot2::geom_line(color = "#0072B2", linewidth = 1) +
      ggplot2::geom_point(color = "#0072B2", size = 3) +
      ggplot2::labs(
         title = "Analytic Power Curve: Dependent Censoring RMST Model",
         subtitle = "Based on the asymptotic variance of the estimator.",
         x = "Sample Size Per Arm", y = "Estimated Power"
      ) +
      ggplot2::ylim(0, 1) + ggplot2::theme_minimal()

   return(list(results_data = results_df, results_plot = p))
}



# Sample Size Search ------------------------------------------------------


#' @title Find Sample Size for RMST Model with Dependent Censoring (Analytic)
#' @description Calculates the required sample size for a target power for an RMST model
#'  with dependent censoring.
#'
#' @details
#' This function performs an iterative search for the sample size needed to
#' achieve a specified `target_power`. It uses the same underlying theory as
#' `DC.power.analytical`. It performs a one-time estimation of the treatment
#' effect and its asymptotic variance from the pilot data, then uses these
#' parameters in an analytic formula to efficiently search for the required sample size.
#'
#' @param pilot_data A `data.frame` containing pilot study data.
#' @param time_var A character string for the time-to-event variable.
#' @param status_var A character string for the primary event status (1=event, 0=otherwise).
#' @param arm_var A character string for the treatment arm variable (1=treatment, 0=control).
#' @param dep_cens_status_var A character string for the dependent censoring status (1=dependent event, 0=otherwise).
#' @param target_power A single numeric value for the desired power.
#' @param linear_terms An optional character vector of other covariate names.
#' @param L The numeric value for the RMST truncation time.
#' @param alpha The significance level (Type I error rate).
#' @param n_start The starting sample size *per arm* for the search.
#' @param n_step The increment in sample size at each step of the search.
#' @param max_n_per_arm The maximum sample size *per arm* to search up to.
#'
#' @return A `list` containing:
#' \item{results_data}{A `data.frame` with the target power and required sample size.}
#' \item{results_plot}{A `ggplot` object visualizing the search path.}
#' \item{results_summary}{A `data.frame` summarizing the estimated treatment effect.}
#' @export
#' @importFrom survival Surv coxph basehaz
#' @importFrom stats lm as.formula complete.cases na.omit sd quantile pnorm qnorm model.matrix coef vcov predict
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline geom_vline labs theme_minimal
#' @importFrom knitr kable
#' @examples
#' # Generate sample pilot data with a clear treatment effect
#' set.seed(456)
#' n_pilot <- 200
#' pilot_df_ss <- data.frame(
#'   time = rexp(n_pilot, rate = 0.2),
#'   arm = rep(0:1, each = n_pilot / 2),
#'   age = rnorm(n_pilot, mean = 55, sd = 8)
#' )
#' # Introduce a treatment effect
#' pilot_df_ss$time[pilot_df_ss$arm == 1] <- pilot_df_ss$time[pilot_df_ss$arm == 1] * 1.5
#'
#' # Create competing event indicators
#' event_type <- sample(0:2, n_pilot, replace = TRUE, prob = c(0.6, 0.2, 0.2))
#' pilot_df_ss$status <- ifelse(event_type == 0, 1, 0)
#' pilot_df_ss$dep_cens_status <- ifelse(event_type == 1, 1, 0)
#' pilot_df_ss$time[event_type != 0] <- pilot_df_ss$time[event_type != 0] * 0.7
#'
#' # Run the sample size search
#' dc_ss_results <- DC.ss.analytical(
#'   pilot_data = pilot_df_ss,
#'   time_var = "time",
#'   status_var = "status",
#'   arm_var = "arm",
#'   dep_cens_status_var = "dep_cens_status",
#'   target_power = 0.80,
#'   linear_terms = "age",
#'   L = 15,
#'   alpha = 0.05,
#'   n_start = 100,
#'   n_step = 50
#' )
#' print(dc_ss_results$results_data)
#' print(dc_ss_results$results_plot)
DC.ss.analytical <- function(pilot_data,
                             time_var,
                             status_var,
                             arm_var,
                             dep_cens_status_var,
                             target_power,
                             linear_terms = NULL,
                             L,
                             alpha = 0.05,
                             n_start = 50,
                             n_step = 25,
                             max_n_per_arm = 2000) {

   # --- 1. Estimate Parameters and Variance from Pilot Data (One Time) ---
   cat("--- Estimating parameters from pilot data for analytic calculation... ---\n")

   core_vars <- c(time_var, status_var, arm_var, dep_cens_status_var)
   covariates <- c(arm_var, linear_terms)
   all_vars <- c(core_vars, covariates)
   df <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]
   n_pilot <- nrow(df)

   model_rhs <- paste(covariates, collapse = " + ")
   cens_formula <- stats::as.formula(paste("survival::Surv(", time_var, ", cens_ind) ~", model_rhs))
   dep_cens_formula <- stats::as.formula(paste("survival::Surv(", time_var, ", dep_cens_ind) ~", model_rhs))
   model_formula <- stats::as.formula(paste("Y_rmst ~", model_rhs))

   message("Model: Y_rmst ~ ", model_rhs)

   df$cens_ind <- df[[status_var]] == 0 & df[[dep_cens_status_var]] == 0
   df$dep_cens_ind <- df[[dep_cens_status_var]] == 1
   df$Y_rmst <- pmin(df[[time_var]], L)

   fit_cens <- tryCatch(survival::coxph(cens_formula, data = df, ties = "breslow"), error = function(e) NULL)
   fit_dep_cens <- tryCatch(survival::coxph(dep_cens_formula, data = df, ties = "breslow"), error = function(e) NULL)

   if(is.null(fit_cens) || is.null(fit_dep_cens)) {
      stop("One or both of the censoring models failed to fit on the pilot data.", call. = FALSE)
   }

   bh_cens <- survival::basehaz(fit_cens, centered = FALSE)
   H_cens <- stats::stepfun(bh_cens$time, c(0, bh_cens$hazard))(df$Y_rmst) * exp(predict(fit_cens, newdata = df, type = "lp"))
   bh_dep_cens <- survival::basehaz(fit_dep_cens, centered = FALSE)
   H_dep_cens <- stats::stepfun(bh_dep_cens$time, c(0, bh_dep_cens$hazard))(df$Y_rmst) * exp(predict(fit_dep_cens, newdata = df, type = "lp"))
   df$weights <- exp(H_cens + H_dep_cens)
   df$weights[df[[status_var]] != 1] <- 0

   finite_weights <- df$weights[is.finite(df$weights) & df$weights > 0]
   if (length(finite_weights) > 0) {
      weight_cap <- stats::quantile(finite_weights, probs = 0.99, na.rm = TRUE)
      df$weights[df$weights > weight_cap] <- weight_cap
   }
   df$weights[!is.finite(df$weights)] <- 0

   fit_data <- df[df$weights > 0, ]
   fit_weights <- fit_data$weights

   if (length(unique(fit_data[[arm_var]])) < 2) {
      stop("Pilot data contains primary events in only one arm after filtering. Cannot estimate treatment effect.", call. = FALSE)
   }

   fit_lm <- stats::lm(model_formula, data = fit_data, weights = fit_weights)
   beta_hat <- stats::coef(fit_lm)
   arm_coeff_name <- arm_var
   beta_effect <- beta_hat[arm_coeff_name]

   X <- stats::model.matrix(model_formula, data = df)
   A_hat <- crossprod(X) / n_pilot
   df$predicted_rmst <- predict(fit_lm, newdata = df)
   residuals <- df$Y_rmst - df$predicted_rmst
   epsilon <- X * residuals * df$weights
   epsilon[is.na(epsilon)] <- 0
   B_hat <- crossprod(epsilon) / n_pilot
   A_hat_inv <- solve(A_hat)
   V_hat_n <- A_hat_inv %*% B_hat %*% t(A_hat_inv)
   var_beta_pilot <- V_hat_n[arm_coeff_name, arm_coeff_name] / n_pilot
   var_beta_n1 <- var_beta_pilot * n_pilot
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

   search_path_df <- data.frame(
      N_per_Arm = as.integer(names(search_path)),
      Power = as.numeric(unlist(search_path))
   )

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
   print(knitr::kable(results_df, caption = "Required Sample Size"))

   return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}


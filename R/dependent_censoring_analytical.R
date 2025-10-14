# Power Calculation -------------------------------------------------------

#' @title Analyze Power for RMST Model with Covariate-Dependent Censoring (Analytic)
#' @description Performs power analysis for an RMST model when the censoring mechanism
#'   depends on observed covariates. Competing risks are not modeled here.
#'
#' @details
#' This function assumes a single censoring process whose hazard depends on covariates
#' (but not on treatment by default). It fits a Cox model for censoring
#' \deqn{\Pr(\text{censoring by } t \mid X) = 1 - G(t \mid X)}
#' using \code{Surv(time, status==0) ~ linear\_terms}, then forms inverse-probability-of-censoring
#' weights (IPCW) \eqn{w_i = 1/\hat G(Y_i\mid X_i)} evaluated at \eqn{Y_i=\min(T_i,L)}.
#' The RMST regression \eqn{E[Y_i \mid A_i,X_i]} is then fit by weighted least squares,
#' and power is derived from a sandwich variance that ignores uncertainty from
#' estimating \eqn{\hat G}.
#'
#' Note: \code{dep_cens_status_var} is accepted for API compatibility but is
#' ignored under this setting (no competing risks are modeled).
#'
#' @param pilot_data A `data.frame` with pilot data.
#' @param time_var Name of the time-to-event variable.
#' @param status_var Name of the primary event indicator (1=event, 0=censored).
#' @param arm_var Name of the treatment indicator (1=treatment, 0=control).
#' @param sample_sizes Numeric vector of per-arm sample sizes for power.
#' @param linear_terms Optional character vector of additional covariate names (used in both models).
#' @param L RMST truncation time.
#' @param alpha Two-sided Type I error.
#'
#' @return A `list` with:
#' \item{results_data}{A data.frame with \code{N_per_Arm} and \code{Power}.}
#' \item{results_plot}{A ggplot object of the power curve.}
#' @export
#' @importFrom stats pnorm qnorm lm as.formula complete.cases quantile model.matrix coef predict
#' @importFrom survival Surv coxph basehaz
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_minimal ylim
DC.power.analytical <- function(pilot_data,
                                time_var,
                                status_var,
                                arm_var,
                                sample_sizes,
                                linear_terms = NULL,
                                L,
                                alpha = 0.05) {

   # --- 1) Prep & outcome ---
   covariates <- c(arm_var, linear_terms)
   all_vars <- unique(c(time_var, status_var, covariates))
   df <- pilot_data[stats::complete.cases(pilot_data[, all_vars, drop = FALSE]), ]
   n_pilot <- nrow(df)
   if (n_pilot < 10) stop("Too few complete cases after filtering.", call. = FALSE)

   df$Y_rmst <- pmin(df[[time_var]], L)

   # --- 2) Censoring model: covariate-dependent only (no arm, no competing risks) ---
   cens_rhs <- if (is.null(linear_terms) || length(linear_terms) == 0) "1" else paste(linear_terms, collapse = " + ")
   cens_formula <- stats::as.formula(
      paste0("survival::Surv(", time_var, ", ", status_var, "==0) ~ ", cens_rhs)
   )
   fit_cens <- survival::coxph(cens_formula, data = df, ties = "breslow")

   bh <- survival::basehaz(fit_cens, centered = FALSE)
   H0_step <- stats::stepfun(bh$time, c(0, bh$hazard))
   lp <- if (cens_rhs == "1") 0 else stats::predict(fit_cens, newdata = df, type = "lp")
   Hc <- H0_step(df$Y_rmst) * exp(lp)
   Ghat <- exp(-Hc)

   # IPCW for all subjects (stabilize)
   eps <- 1e-6
   w <- 1 / pmax(Ghat, eps)
   w[!is.finite(w)] <- 0
   if (any(is.finite(w) & w > 0)) {
      cap <- stats::quantile(w[is.finite(w) & w > 0], 0.99, na.rm = TRUE)
      w[w > cap] <- cap
   }
   df$w <- w

   # --- 3) Weighted RMST regression ---
   model_rhs <- paste(c(arm_var, linear_terms), collapse = " + ")
   model_formula <- stats::as.formula(paste("Y_rmst ~", model_rhs))
   fit_wls <- stats::lm(model_formula, data = df, weights = df$w)
   beta <- stats::coef(fit_wls)
   if (!(arm_var %in% names(beta))) stop("Arm coefficient not found in model.", call. = FALSE)
   beta_arm <- beta[arm_var]

   # --- 4) Sandwich variance (A=X'WX/n, B = X'(w*r)^2 X / n) ---
   X <- stats::model.matrix(model_formula, data = df)
   r <- df$Y_rmst - as.numeric(X %*% beta)
   W <- df$w

   A_hat <- crossprod(X, X * W) / n_pilot
   meat <- X * (W * r)
   B_hat <- crossprod(meat) / n_pilot

   Ainv <- solve(A_hat)
   Vn <- Ainv %*% B_hat %*% Ainv
   var_beta_pilot <- Vn[colnames(X) == arm_var, colnames(X) == arm_var] / n_pilot
   var_beta_n1 <- as.numeric(var_beta_pilot * n_pilot)
   se_beta_n1 <- sqrt(var_beta_n1)

   # --- 5) Power across N ---
   z_alpha <- stats::qnorm(1 - alpha/2)
   power_values <- sapply(sample_sizes, function(n_per_arm) {
      N <- 2 * n_per_arm
      seN <- se_beta_n1 / sqrt(N)
      stats::pnorm((abs(beta_arm) / seN) - z_alpha)
   })

   results_df <- data.frame(N_per_Arm = sample_sizes, Power = power_values)

   p <- ggplot2::ggplot(results_df, ggplot2::aes(x = N_per_Arm, y = Power)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 3) +
      ggplot2::labs(
         title = "Analytic Power: Covariate-Dependent Censoring (IPCW)",
         subtitle = "Single censoring mechanism; no competing risks.",
         x = "Sample Size Per Arm", y = "Estimated Power"
      ) +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_minimal()

   list(results_data = results_df, results_plot = p)
}



# Sample Size Search ------------------------------------------------------

#' @title Find Sample Size for RMST with Covariate-Dependent Censoring (Analytic)
#' @description Iteratively finds required per-arm sample size for a target power,
#'   using the same IPCW-based analytic variance as \code{DC.power.analytical}.
#'
#' @details
#' Uses a single censoring Cox model \code{Surv(time, status==0) ~ linear_terms} to form IPCW
#' and fits a weighted RMST regression. Treatment is excluded from the censoring model
#' by default. Competing risks are not modeled. Variance ignores uncertainty in \eqn{\hat G}.
#'
#' Note: \code{dep_cens_status_var} is accepted for API compatibility but ignored here.
#'
#' @param pilot_data A `data.frame` containing pilot study data.
#' @param time_var Time variable.
#' @param status_var Event indicator (1=event, 0=censored).
#' @param arm_var Treatment indicator (1/0).
#' @param target_power Desired power.
#' @param linear_terms Optional covariates used in both models.
#' @param L RMST truncation time.
#' @param alpha Two-sided Type I error.
#' @param n_start Starting per-arm N.
#' @param n_step Step size for search.
#' @param max_n_per_arm Maximum per-arm N to search.
#'
#' @return A `list` with:
#' \item{results_data}{data.frame with \code{Target_Power} and \code{Required_N_per_Arm}.}
#' \item{results_plot}{ggplot showing the search path.}
#' \item{results_summary}{data.frame summarizing the pilot arm effect.}
#' @export
#' @importFrom survival Surv coxph basehaz
#' @importFrom stats lm as.formula complete.cases quantile model.matrix coef predict pnorm qnorm
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline geom_vline labs theme_minimal
#' @importFrom knitr kable
DC.ss.analytical <- function(pilot_data,
                             time_var,
                             status_var,
                             arm_var,
                             target_power,
                             linear_terms = NULL,
                             L,
                             alpha = 0.05,
                             n_start = 50,
                             n_step = 25,
                             max_n_per_arm = 2000) {

   # --- One-time estimation identical to DC.power.analytical ---
   covariates <- c(arm_var, linear_terms)
   all_vars <- unique(c(time_var, status_var, covariates))
   df <- pilot_data[stats::complete.cases(pilot_data[, all_vars, drop = FALSE]), ]
   n_pilot <- nrow(df)
   if (n_pilot < 10) stop("Too few complete cases after filtering.", call. = FALSE)

   df$Y_rmst <- pmin(df[[time_var]], L)

   cens_rhs <- if (is.null(linear_terms) || length(linear_terms) == 0) "1" else paste(linear_terms, collapse = " + ")
   cens_formula <- stats::as.formula(
      paste0("survival::Surv(", time_var, ", ", status_var, "==0) ~ ", cens_rhs)
   )
   fit_cens <- survival::coxph(cens_formula, data = df, ties = "breslow")

   bh <- survival::basehaz(fit_cens, centered = FALSE)
   H0_step <- stats::stepfun(bh$time, c(0, bh$hazard))
   lp <- if (cens_rhs == "1") 0 else stats::predict(fit_cens, newdata = df, type = "lp")
   Hc <- H0_step(df$Y_rmst) * exp(lp)
   Ghat <- exp(-Hc)

   eps <- 1e-6
   w <- 1 / pmax(Ghat, eps)
   w[!is.finite(w)] <- 0
   if (any(is.finite(w) & w > 0)) {
      cap <- stats::quantile(w[is.finite(w) & w > 0], 0.99, na.rm = TRUE)
      w[w > cap] <- cap
   }
   df$w <- w

   model_rhs <- paste(c(arm_var, linear_terms), collapse = " + ")
   model_formula <- stats::as.formula(paste("Y_rmst ~", model_rhs))
   fit_wls <- stats::lm(model_formula, data = df, weights = df$w)
   beta <- stats::coef(fit_wls)
   if (!(arm_var %in% names(beta))) stop("Arm coefficient not found in model.", call. = FALSE)
   beta_arm <- beta[arm_var]

   X <- stats::model.matrix(model_formula, data = df)
   r <- df$Y_rmst - as.numeric(X %*% beta)
   W <- df$w

   A_hat <- crossprod(X, X * W) / n_pilot
   meat <- X * (W * r)
   B_hat <- crossprod(meat) / n_pilot

   Ainv <- solve(A_hat)
   Vn <- Ainv %*% B_hat %*% Ainv
   var_beta_pilot <- Vn[colnames(X) == arm_var, colnames(X) == arm_var] / n_pilot
   var_beta_n1 <- as.numeric(var_beta_pilot * n_pilot)
   se_beta_n1 <- sqrt(var_beta_n1)

   # --- Iterative search ---
   z_alpha <- stats::qnorm(1 - alpha / 2)
   current_n <- n_start
   path <- list()
   final_n <- NA_integer_

   while (current_n <= max_n_per_arm) {
      N <- 2 * current_n
      seN <- se_beta_n1 / sqrt(N)
      pw <- stats::pnorm((abs(beta_arm) / seN) - z_alpha)
      path[[as.character(current_n)]] <- pw
      if (pw >= target_power) { final_n <- current_n; break }
      current_n <- current_n + n_step
   }
   if (is.na(final_n)) {
      warning(paste("Target power", target_power, "not achieved by max N of", max_n_per_arm), call. = FALSE)
      final_n <- max_n_per_arm
   }

   results_df <- data.frame(Target_Power = target_power, Required_N_per_Arm = final_n)
   search_path_df <- data.frame(
      N_per_Arm = as.integer(names(path)),
      Power = as.numeric(unlist(path))
   )

   p <- ggplot2::ggplot(na.omit(search_path_df), ggplot2::aes(x = N_per_Arm, y = Power)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_hline(yintercept = target_power, linetype = "dashed", color = "red") +
      ggplot2::geom_vline(xintercept = final_n, linetype = "dotted", color = "blue") +
      ggplot2::labs(
         title = "Analytic Sample Size Search (Covariate-Dependent Censoring; IPCW)",
         subtitle = "Single censoring mechanism; no competing risks.",
         x = "Sample Size Per Arm", y = "Calculated Power"
      ) +
      ggplot2::theme_minimal()

   cat("\n--- Calculation Summary ---\n")
   print(knitr::kable(results_df, caption = "Required Sample Size"))

   list(results_data = results_df, results_plot = p,
        results_summary = data.frame(Statistic = "Assumed RMST Difference (from pilot)",
                                     Value = beta_arm))
}

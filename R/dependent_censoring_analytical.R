# Power Calculation -------------------------------------------------------

#' @title Analyze Power for RMST Model with Covariate-Dependent Censoring (Analytic)
#' @description Performs power analysis for an RMST model when the censoring mechanism
#'   depends on observed covariates.
#'
#' @details
#' This function assumes a single censoring process whose hazard depends on covariates
#' (but not on treatment by default). It fits a Cox model for censoring
#' \deqn{\Pr(\text{censoring by } t \mid X) = 1 - G(t \mid X)}
#' using \code{Surv(time, status==0) ~ linear\_terms}, then forms inverse-probability-of-censoring
#' weights (IPCW) \eqn{w_i = \Delta_i^Y/\hat G(Y_i\mid X_i)} evaluated at
#' \eqn{Y_i=\min(T_i,L)}, where \eqn{\Delta_i^Y=1} if the event is observed before
#' \eqn{L} or follow-up reaches \eqn{L}. For numerical stability the weights are
#' capped at their 99th percentile; the cap value and the fraction of weights
#' affected are reported in \code{model_output$censoring_weights}.
#' The RMST regression \eqn{E[Y_i \mid A_i,X_i]} is then fit by weighted least squares,
#' and power is derived from a sandwich variance that ignores uncertainty from
#' estimating \eqn{\hat G}.
#'
#' Note: This implementation models a single censoring process and does not
#' handle competing risks.
#'
#' @param pilot_data A `data.frame` with pilot study data.
#' @param time_var A character string for the time-to-event variable.
#' @param status_var A character string for the event status variable (1=event, 0=censored).
#' @param arm_var A character string for the treatment arm variable (1=treatment, 0=control).
#' @param sample_sizes A numeric vector of sample sizes *per arm* to calculate power for.
#' @param linear_terms Optional character vector of additional covariate names.
#' @param L The numeric value for the RMST truncation time.
#' @param alpha The significance level (Type I error rate).
#' @param verbose Logical; if \code{TRUE}, emit progress messages. Default \code{FALSE}.
#'
#' @return A `list` with:
#' \item{results_data}{A data.frame with \code{N_per_Arm} and \code{Power}.}
#' \item{results_plot}{A ggplot object of the power curve.}
#' \item{results_summary}{A data.frame summarizing the pilot arm effect.}
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
                                alpha = 0.05,
                                verbose = FALSE) {

   # --- 1. Estimate model parameters and sandwich variance from pilot data ---
   est <- .estimate_dc_params(pilot_data, time_var, status_var, arm_var,
                              linear_terms, L)

   # --- 2. Power across N ---
   z_alpha <- stats::qnorm(1 - alpha / 2)
   power_values <- sapply(sample_sizes, function(n_per_arm) {
      .rmst_wald_power(est$beta_effect, est$se_beta_n1, 2 * n_per_arm, z_alpha)
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

   results_summary <- data.frame(
      Statistic = "Assumed RMST Difference (from pilot)",
      Value = est$beta_effect
   )

   list(results_data = results_df, results_plot = p,
        results_summary = results_summary,
        model_output = .dc_model_output(est, arm_var))
}



# Sample Size Search ------------------------------------------------------

#' @title Find Sample Size for RMST with Covariate-Dependent Censoring (Analytic)
#' @description Iteratively finds required per-arm sample size for a target power,
#'   using the same IPCW-based analytic variance as \code{DC.power.analytical}.
#'
#' @details
#' Uses a single censoring Cox model \code{Surv(time, status==0) ~ linear_terms} to form
#' IPCW weights \eqn{\Delta_i^Y/\hat G(Y_i\mid X_i)} and fits a weighted RMST regression
#' among subjects whose truncated RMST outcome is observed. Treatment is excluded from
#' the censoring model by default. Competing risks are not modeled. Variance ignores
#' uncertainty in \eqn{\hat G}. As in \code{DC.power.analytical}, the IPCW
#' weights are capped at their 99th percentile for numerical stability.
#'
#' Note: This implementation models a single censoring process and does not
#' handle competing risks.
#'
#' @param pilot_data A `data.frame` containing pilot study data.
#' @param time_var A character string for the time-to-event variable.
#' @param status_var A character string for the event status variable (1=event, 0=censored).
#' @param arm_var A character string for the treatment arm variable (1=treatment, 0=control).
#' @param target_power A single numeric value for the desired power.
#' @param linear_terms Optional character vector of additional covariate names.
#' @param L The numeric value for the RMST truncation time.
#' @param alpha The significance level (Type I error rate).
#' @param n_start The starting sample size *per arm* for the search.
#' @param n_step The increment in sample size at each step of the search.
#' @param max_n_per_arm The maximum sample size *per arm* to search.
#' @param verbose Logical; if \code{TRUE}, emit progress messages. Default \code{FALSE}.
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
                             max_n_per_arm = 2000,
                             verbose = FALSE) {

   # --- One-time estimation identical to DC.power.analytical ---
   est <- .estimate_dc_params(pilot_data, time_var, status_var, arm_var,
                              linear_terms, L)

   # --- Iterative search ---
   search <- .rmst_analytic_ss_search(est$beta_effect, est$se_beta_n1, 2,
                                      target_power, alpha, n_start, n_step,
                                      max_n_per_arm, "/arm", verbose)
   final_n <- search$final_n

   results_df <- data.frame(Target_Power = target_power, Required_N_per_Arm = final_n)
   search_path_df <- search$search_path_df
   names(search_path_df) <- c("N_per_Arm", "Power")

   p <- ggplot2::ggplot(stats::na.omit(search_path_df), ggplot2::aes(x = N_per_Arm, y = Power)) +
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

   .rmst_verbose_message(verbose, "Calculation summary available in returned results_data.")

   results_summary_ss <- data.frame(
      Statistic = "Assumed RMST Difference (from pilot)", Value = est$beta_effect)

   list(results_data = results_df, results_plot = p,
        results_summary = results_summary_ss,
        model_output = .dc_model_output(est, arm_var))
}

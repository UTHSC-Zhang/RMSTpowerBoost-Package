# Power Calculation -------------------------------------------------------------------


#' @title Analyze Power for a Linear RMST Model (Analytic)
#' @description Performs power analysis using a direct formula based on the
#'   asymptotic variance estimator for the linear RMST model.
#'
#' @details
#' This function implements the analytic power calculation for the direct linear
#' regression model of the Restricted Mean Survival Time (RMST) proposed by Tian et al. (2014).
#' The core of the method is a weighted linear model of the form
#' \deqn{E[Y_i \mid A_i, \mathbf{Z}_i] = \alpha + \tau A_i + \mathbf{Z}_i^\top \boldsymbol{\gamma}}
#' where \eqn{Y_i = \min(T_i, L)} is the event time truncated at \eqn{L}, \eqn{A_i} is the
#' treatment indicator, and \eqn{\tau} is the treatment effect of interest.
#'
#' To handle right-censoring, the method uses Inverse Probability of Censoring
#' Weighting (IPCW). Let \eqn{Y_i = \min(T_i, L)} and
#' \eqn{\Delta_i^Y = 1} if the event occurs before \eqn{L} or follow-up reaches
#' \eqn{L}. The weight is \eqn{w_i = \Delta_i^Y / \hat{G}(Y_i)}, where
#' \eqn{\hat{G}(t) = P(C > t)} is the Kaplan-Meier estimate of the censoring
#' distribution fit on the original time scale. For numerical stability the
#' weights are capped at their 99th percentile; the cap value and the fraction
#' of weights affected are reported in \code{model_output$censoring_weights}.
#'
#' Power is calculated analytically based on the asymptotic properties of the
#' coefficient estimators. The variance of the treatment effect estimator, \eqn{\hat{\tau}}, is derived from a
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
#' @param verbose Logical; if \code{TRUE}, emit progress messages. Default \code{FALSE}.
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
                                    sample_sizes, linear_terms = NULL, L, alpha = 0.05,
                                    verbose = FALSE) {

   # --- 1. Estimate model parameters and sandwich variance from pilot data ---
   .rmst_verbose_message(verbose, "--- Estimating parameters from pilot data for analytic calculation... ---")
   est <- .estimate_linear_params(pilot_data, time_var, status_var, arm_var,
                                  linear_terms, L, verbose)
   .rmst_verbose_message(verbose, "--- Calculating asymptotic variance... ---")

   # --- 2. Calculate Power for Each Sample Size ---
   .rmst_verbose_message(verbose, "--- Calculating power for specified sample sizes... ---")
   z_alpha <- stats::qnorm(1 - alpha / 2)
   power_values <- sapply(sample_sizes, function(n_per_arm) {
      .rmst_wald_power(est$beta_effect, est$se_beta_n1, n_per_arm * 2, z_alpha)
   })

   results_df <- data.frame(N_per_Arm = sample_sizes, Power = power_values)

   results_summary <- data.frame(
      Statistic = "Assumed RMST Difference (from pilot)",
      Value = est$beta_effect
   )

   # --- 3. Create Plot and Return ---
   p <- .rmst_power_curve_plot(
      results_df, "N_per_Arm", "#D55E00",
      title = "Analytic Power Curve: Linear IPCW RMST Model",
      subtitle = "Based on the asymptotic variance from Tian et al. (2014).",
      xlab = "Sample Size Per Arm")

   return(list(results_data = results_df, results_plot = p,
               results_summary = results_summary,
               model_output = .linear_model_output(est, arm_var)))
}

# Sample Size Search ------------------------------------------------------


#' @title Find Sample Size for a Linear RMST Model (Analytic)
#' @description Calculates the required sample size for a target power using an
#'   analytic formula based on the methods of Tian et al. (2014).
#'
#' @details
#' This function performs an iterative search to find the sample size needed to
#' achieve a specified `target_power`. It uses the same underlying theory as
#' `linear.power.analytical`, including the 99th-percentile IPCW weight cap
#' described there. First, it estimates the treatment effect size and its
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
#' @param verbose Logical; if \code{TRUE}, emit progress messages. Default \code{FALSE}.
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
                                 n_start = 50, n_step = 25, max_n_per_arm = 2000,
                                 verbose = FALSE) {

   # --- 1. Estimate Parameters and Variance from Pilot Data (One Time) ---
   .rmst_verbose_message(verbose, "--- Estimating parameters from pilot data for analytic search... ---")
   est <- .estimate_linear_params(pilot_data, time_var, status_var, arm_var,
                                  linear_terms, L, verbose)

   # --- 2. Iterative Search for Sample Size using Analytic Formula ---
   .rmst_verbose_message(verbose, "--- Searching for Sample Size (Method: Analytic) ---")
   search <- .rmst_analytic_ss_search(est$beta_effect, est$se_beta_n1, 2,
                                      target_power, alpha, n_start, n_step,
                                      max_n_per_arm, "/arm", verbose)
   final_n <- search$final_n

   # --- 3. Finalize and Return Results ---
   results_summary <- data.frame(
      Statistic = "Assumed RMST Difference (from pilot)",
      Value = est$beta_effect
   )
   results_df <- data.frame(Target_Power = target_power, Required_N_per_Arm = final_n)
   search_path_df <- search$search_path_df
   names(search_path_df) <- c("N_per_Arm", "Power")

   p <- .rmst_ss_search_plot(
      search_path_df, "N_per_Arm", final_n, target_power,
      title = "Analytic Sample Size Search: Linear IPCW RMST Model",
      subtitle = "Power calculated from formula at each step.",
      xlab = "Sample Size Per Arm")

   .rmst_verbose_message(verbose, "Calculation summary available in returned results_data.")

   return(list(results_data = results_df, results_plot = p,
               results_summary = results_summary,
               model_output = .linear_model_output(est, arm_var)))
}

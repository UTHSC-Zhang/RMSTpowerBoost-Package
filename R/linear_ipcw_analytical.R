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
#' \eqn{\hat{G}(t)} estimates censoring survival conditional on treatment and
#' \code{linear_terms}, using a Cox proportional hazards model with Breslow
#' ties fit on the original time scale and evaluated at each subject's
#' truncated outcome. The weights are used as fitted, without capping or
#' truncation, so that the estimator and its variance reproduce the published
#' formulas exactly; their distribution is reported in
#' \code{model_output$censoring_weights$raw_summary}.
#'
#' Power is calculated analytically based on the asymptotic properties of the
#' coefficient estimators. The variance of the treatment effect estimator, \eqn{\hat{\tau}}, is derived from a
#' robust sandwich variance estimator of the form \eqn{A^{-1}B(A^{-1})'}. This is the
#' simplified sandwich variance of Zhang and Schaubel (2024, Equations 11-12), which
#' treats the fitted IPCW weights as fixed: `A` is the weighted information matrix
#' \eqn{(X'WX)/n}, with the weights entering linearly, and `B` is the empirical second
#' moment of the weighted score contributions,
#' \eqn{(\sum \epsilon_i \epsilon_i')/n} with \eqn{\epsilon_i = w_i x_i (Y_i - x_i'\hat{\beta})}.
#' The intercept-plus-slopes parameterization used here is, for a single stratum,
#' algebraically equivalent to the weighted-centered-covariate form without intercept
#' used in the paper. As in the paper's simulations and data analysis, the additional
#' influence-function terms for estimation of the Cox censoring model given in their
#' Theorem 3.1 are not included.
#'
#' The `test` argument selects which hypothesis test the reported power refers to, and
#' this matters because the two candidates do not use the same variance. The sandwich
#' above is consistent for the sampling variability of the estimator. The ordinary
#' weighted-least-squares standard error reported by `summary(lm(...))`, which is the
#' test that `linear.power.boot()` simulates, is not, because IPCW weights are sampling
#' weights rather than inverse-variance weights. In a representative pilot the WLS
#' standard error runs about 20 percent above the true sampling standard deviation of
#' the estimate, which pushes the rejection threshold out and costs power.
#'
#' With `test = "wls"` (the default) the power refers to the weighted least-squares
#' t-test. The sampling distribution of the estimate is taken from the sandwich
#' variance \eqn{\sigma_{SW}} and the rejection threshold from the WLS standard error
#' \eqn{\sigma_{WLS}}, giving
#' \deqn{\Phi\left(\frac{|\tau| - z_{1-\alpha/2}\,\sigma_{WLS}}{\sigma_{SW}}\right) +
#'       \Phi\left(\frac{-|\tau| - z_{1-\alpha/2}\,\sigma_{WLS}}{\sigma_{SW}}\right).}
#' This is the analysis that `linear.power.boot()` performs, so the analytic and
#' bootstrap engines then describe the same test.
#'
#' With `test = "sandwich"` the power refers to the sandwich Wald test, in which one
#' variance sets both the sampling distribution and the threshold, giving the usual
#' \eqn{\Phi(|\tau|/\sigma_{SW} - z_{1-\alpha/2})}. That is the more efficient test and
#' the one matching the paper's inference, but it is not what the bootstrap simulates.
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
#' @param test Which test the reported power refers to. `"wls"` (default) is the
#'   weighted least-squares t-test, matching `linear.power.boot`. `"sandwich"` is
#'   the sandwich Wald test of Zhang and Schaubel (2024). See Details.
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
                                    verbose = FALSE, test = c("wls", "sandwich")) {
   test <- match.arg(test)

   # --- 1. Estimate model parameters and sandwich variance from pilot data ---
   .rmst_verbose_message(verbose, "--- Estimating parameters from pilot data for analytic calculation... ---")
   est <- .estimate_linear_params(pilot_data, time_var, status_var, arm_var,
                                  linear_terms, L, verbose)
   .rmst_verbose_message(verbose, "--- Calculating asymptotic variance... ---")

   # --- 2. Calculate Power for Each Sample Size ---
   .rmst_verbose_message(verbose, "--- Calculating power for specified sample sizes... ---")
   z_alpha <- stats::qnorm(1 - alpha / 2)
   se_test_n1 <- if (test == "wls") est$se_beta_model_n1 else NULL
   power_values <- sapply(sample_sizes, function(n_per_arm) {
      .rmst_wald_power(est$beta_effect, est$se_beta_n1, n_per_arm * 2, z_alpha,
                       se_test_n1)
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
      subtitle = if (test == "wls")
         "Weighted least-squares t-test; sandwich sampling variance."
      else "Simplified sandwich test (Zhang & Schaubel, 2024).",
      xlab = "Sample Size Per Arm")

   return(list(results_data = results_df, results_plot = p,
               results_summary = results_summary,
               model_output = .linear_model_output(est, arm_var, test)))
}

# Sample Size Search ------------------------------------------------------


#' @title Find Sample Size for a Linear RMST Model (Analytic)
#' @description Calculates the required sample size for a target power using an
#'   analytic formula based on the methods of Tian et al. (2014).
#'
#' @details
#' This function performs an iterative search to find the sample size needed to
#' achieve a specified `target_power`. It uses the same underlying theory as
#' `linear.power.analytical`, including the uncapped IPCW weights, the
#' simplified sandwich variance, and the `test` argument described there. First, it estimates the treatment effect size and its
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
#' @param test Which test the required sample size refers to. `"wls"` (default)
#'   is the weighted least-squares t-test, matching `linear.ss.boot`.
#'   `"sandwich"` is the sandwich Wald test. See `linear.power.analytical`.
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
                                 verbose = FALSE, test = c("wls", "sandwich")) {
   test <- match.arg(test)

   # --- 1. Estimate Parameters and Variance from Pilot Data (One Time) ---
   .rmst_verbose_message(verbose, "--- Estimating parameters from pilot data for analytic search... ---")
   est <- .estimate_linear_params(pilot_data, time_var, status_var, arm_var,
                                  linear_terms, L, verbose)

   # --- 2. Iterative Search for Sample Size using Analytic Formula ---
   .rmst_verbose_message(verbose, "--- Searching for Sample Size (Method: Analytic) ---")
   search <- .rmst_analytic_ss_search(est$beta_effect, est$se_beta_n1, 2,
                                      target_power, alpha, n_start, n_step,
                                      max_n_per_arm, "/arm", verbose,
                                      if (test == "wls") est$se_beta_model_n1 else NULL)
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
               model_output = .linear_model_output(est, arm_var, test)))
}

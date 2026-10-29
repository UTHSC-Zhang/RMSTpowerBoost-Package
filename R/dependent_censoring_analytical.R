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
#' \eqn{L} or follow-up reaches \eqn{L}. The censoring model uses Breslow ties, is fit
#' on the original follow-up, and is evaluated at each subject's truncated outcome,
#' with the uncentered baseline hazard paired with an uncentered linear predictor.
#' The weights are used as fitted, without capping, truncation, or a floor on
#' \eqn{\hat G}, so that the estimator and its variance reproduce the published
#' formulas exactly; their distribution is reported in
#' \code{model_output$censoring_weights$raw_summary}. With no observed censoring the
#' weights are one and no censoring model is fit.
#' The RMST regression \eqn{E[Y_i \mid A_i,X_i]} is then fit by weighted least squares,
#' and the sampling variance of the treatment effect is the sandwich
#' \eqn{A^{-1}B(A^{-1})'} with the fitted weights treated as fixed, so it ignores
#' uncertainty from estimating \eqn{\hat G}.
#'
#' The `test` argument selects which hypothesis test the reported power refers to,
#' and this matters because the two candidates do not use the same variance. The
#' sandwich above is consistent for the sampling variability of the estimator. The
#' ordinary weighted-least-squares standard error reported by `summary(lm(...))` is
#' not, because IPCW weights are sampling weights rather than inverse-variance
#' weights, so it typically runs above the true sampling standard deviation of the
#' estimate, which pushes the rejection threshold out and costs power.
#'
#' With `test = "wls"` (the default) the power refers to the weighted least-squares
#' t-test that a standard `lm(..., weights = w)` analysis reports, as in
#' `linear.power.analytical()`. The sampling distribution of the estimate is taken
#' from the sandwich variance \eqn{\sigma_{SW}} and the rejection threshold from the
#' WLS standard error \eqn{\sigma_{WLS}}, giving
#' \deqn{\Phi\left(\frac{|\tau| - z_{1-\alpha/2}\,\sigma_{WLS}}{\sigma_{SW}}\right) +
#'       \Phi\left(\frac{-|\tau| - z_{1-\alpha/2}\,\sigma_{WLS}}{\sigma_{SW}}\right),}
#' and \code{model_output$coefficient_table} is the t-table that `summary()`
#' reports for the weighted fit.
#'
#' With `test = "sandwich"` the power refers to the sandwich Wald test, in which one
#' variance sets both the sampling distribution and the threshold, giving the usual
#' \eqn{\Phi(|\tau|/\sigma_{SW} - z_{1-\alpha/2})}. That is the more efficient test,
#' and \code{model_output$coefficient_table} then uses the robust variance. Both
#' standard errors are always returned in \code{model_output$variance_components}.
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
#' @param test Which test the reported power refers to. `"wls"` (default) is the
#'   weighted least-squares t-test. `"sandwich"` is the robust sandwich Wald test.
#'   See Details.
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
                                verbose = FALSE,
                                test = c("wls", "sandwich")) {
   test <- match.arg(test)

   # --- 1. Estimate model parameters and sandwich variance from pilot data ---
   est <- .estimate_dc_params(pilot_data, time_var, status_var, arm_var,
                              linear_terms, L)

   # --- 2. Power across N ---
   z_alpha <- stats::qnorm(1 - alpha / 2)
   se_test_n1 <- if (test == "wls") est$se_beta_model_n1 else NULL
   power_values <- sapply(sample_sizes, function(n_per_arm) {
      .rmst_wald_power(est$beta_effect, est$se_beta_n1, 2 * n_per_arm, z_alpha,
                       se_test_n1)
   })

   results_df <- data.frame(N_per_Arm = sample_sizes, Power = power_values)

   p <- ggplot2::ggplot(results_df, ggplot2::aes(x = N_per_Arm, y = Power)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 3) +
      ggplot2::labs(
         title = "Analytic Power: Covariate-Dependent Censoring (IPCW)",
         subtitle = if (test == "wls")
            "Weighted least-squares t-test; sandwich sampling variance. Single censoring mechanism."
         else "Sandwich Wald test. Single censoring mechanism; no competing risks.",
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
        model_output = .dc_model_output(est, arm_var, test))
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
#' uncertainty in \eqn{\hat G}. It uses the same uncapped IPCW weights, sandwich
#' sampling variance, and `test` argument as \code{DC.power.analytical}.
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
#' @param test Which test the required sample size refers to. `"wls"` (default)
#'   is the weighted least-squares t-test; `"sandwich"` is the robust sandwich Wald
#'   test. See `DC.power.analytical`.
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
                             verbose = FALSE,
                             test = c("wls", "sandwich")) {
   test <- match.arg(test)

   # --- One-time estimation identical to DC.power.analytical ---
   est <- .estimate_dc_params(pilot_data, time_var, status_var, arm_var,
                              linear_terms, L)

   # --- Iterative search ---
   search <- .rmst_analytic_ss_search(est$beta_effect, est$se_beta_n1, 2,
                                      target_power, alpha, n_start, n_step,
                                      max_n_per_arm, "/arm", verbose,
                                      if (test == "wls") est$se_beta_model_n1 else NULL)
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
        model_output = .dc_model_output(est, arm_var, test))
}

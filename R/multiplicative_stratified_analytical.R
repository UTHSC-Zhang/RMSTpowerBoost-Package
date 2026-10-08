# Power Calculation -------------------------------------------------------

#' @title Analyze Power for a Multiplicative Stratified RMST Model (Analytic)
#' @description Performs power analysis for a multiplicative, stratified RMST model using an
#'   analytic method based on the work of Wang et al. (2019).
#'
#' @details
#' This function is based on the method for evaluating center (stratum) effects
#' using a multiplicative model for RMST: \eqn{\mu_{ij} = \mu_{0j} \exp\{\beta'Z_i\}}.
#' The method uses IPCW with a stratified Cox model for the censoring distribution
#' fit on the original time scale and weights of the form
#' \eqn{\hat{W}_{ij}\Delta_i^Y}, where \eqn{\Delta_i^Y = 1} if the event occurs
#' before \eqn{L} or follow-up reaches \eqn{L}. The stratified Cox model uses
#' Breslow ties and is evaluated at each subject's truncated outcome. The
#' weights are used as fitted, without capping or truncation, so that the
#' estimator and its variance reproduce the published formulas exactly; their
#' distribution is reported in
#' \code{model_output$censoring_weights$raw_summary}.
#'
#' Formal estimation of \eqn{\beta} requires an iterative solver for the estimating
#' equation given in Equation (8) of Wang et al. (2019). Because this is computationally
#' intensive, this implementation uses a log-linear working model fitted by weighted
#' least squares to pseudo-observations (\code{lm(log(Y_rmst) ~ ...)}) as a tractable
#' approximation. The approximation is consistent when the log-linear mean structure
#' is well-specified but may differ from the formal estimator under strong misspecification.
#'
#' The power calculation relies on the asymptotic variance of the log-RMST ratio
#' estimator, \eqn{\hat{\tau}}. The variance is the robust sandwich estimator
#' \eqn{A_n^{-1} B_n (A_n^{-1})'} of the weighted least-squares estimating
#' equations, with \eqn{A_n = X'WX/n} and \eqn{B_n} the empirical second moment
#' of the weighted residuals, mirroring the sandwich form described in Theorem 1
#' of Wang et al. (2019).
#'
#' The `test` argument selects which hypothesis test the reported power refers to,
#' and this matters because the two candidates do not use the same variance. The
#' sandwich above is consistent for the sampling variability of the estimator. The
#' ordinary weighted-least-squares standard error reported by `summary(lm(...))` is
#' not, because IPCW weights are sampling weights rather than inverse-variance
#' weights, so it typically runs above the true sampling standard deviation of the
#' estimate, which pushes the rejection threshold out and costs power.
#'
#' With `test = "sandwich"` (the default) the power refers to the sandwich Wald test,
#' in which one variance sets both the sampling distribution and the threshold,
#' giving the usual \eqn{\Phi(|\tau|/\sigma_{SW} - z_{1-\alpha/2})}. This is the
#' inference of Wang et al. (2019), and coefficient standard errors in
#' \code{model_output$coefficient_table} then use the robust variance.
#'
#' With `test = "wls"` the power refers to the weighted least-squares t-test: the
#' sampling distribution of the estimate is taken from the sandwich variance
#' \eqn{\sigma_{SW}} and the rejection threshold from the WLS standard error
#' \eqn{\sigma_{WLS}}, giving
#' \deqn{\Phi\left(\frac{|\tau| - z_{1-\alpha/2}\,\sigma_{WLS}}{\sigma_{SW}}\right) +
#'       \Phi\left(\frac{-|\tau| - z_{1-\alpha/2}\,\sigma_{WLS}}{\sigma_{SW}}\right),}
#' and \code{model_output$coefficient_table} becomes the t-table that
#' `summary()` reports for the weighted fit. Both standard errors are always
#' returned in \code{model_output$variance_components}.
#'
#' Note that the default differs from `linear.power.analytical()`, which defaults to
#' `test = "wls"` because `linear.power.boot()` simulates exactly that test. No
#' bootstrap engine in this package simulates the IPCW-weighted t-test for a
#' stratified model: `MS.power.boot()` fits an unweighted model to jackknife
#' pseudo-observations. `"sandwich"` is therefore the default here.
#'
#' @param pilot_data A `data.frame` with pilot study data.
#' @param time_var A character string for the time-to-event variable.
#' @param status_var A character string for the event status variable (1=event, 0=censored).
#' @param arm_var A character string for the treatment arm variable (1=treatment, 0=control).
#' @param strata_var A character string for the stratification variable.
#' @param sample_sizes A numeric vector of sample sizes *per stratum* to calculate power for.
#' @param linear_terms An optional character vector of other covariate names.
#' @param L The numeric value for the RMST truncation time.
#' @param alpha The significance level (Type I error rate).
#' @param verbose Logical; if \code{TRUE}, emit progress messages. Default \code{FALSE}.
#' @param test Which test the reported power refers to. `"sandwich"` (default) is
#'   the robust sandwich Wald test of Wang et al. (2019). `"wls"` is the weighted
#'   least-squares t-test. See Details.
#'
#' @return A `list` containing:
#' \item{results_data}{A `data.frame` with sample sizes and corresponding powers.}
#' \item{results_plot}{A `ggplot` object visualizing the power curve.}
#' @export
#' @importFrom survival Surv coxph basehaz
#' @importFrom stats lm as.formula complete.cases na.omit sd quantile pnorm qnorm coef vcov predict model.matrix residuals
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
#'  L = 10, alpha = 0.05
#' )
#' print(power_results$results_data)
MS.power.analytical <- function(pilot_data, time_var, status_var, arm_var, strata_var,
                                sample_sizes, linear_terms = NULL, L, alpha = 0.05,
                                verbose = FALSE, test = c("sandwich", "wls")) {
   test <- match.arg(test)

   # --- 1. Estimate Parameters from Pilot Data ---
   .rmst_verbose_message(verbose, "--- Estimating parameters from pilot data (log-linear approximation)... ---")
   est <- .estimate_ms_params(pilot_data, time_var, status_var, arm_var,
                              strata_var, linear_terms, L)

   # --- 2. Calculate Power ---
   .rmst_verbose_message(verbose, "--- Calculating power for specified sample sizes... ---")
   z_alpha <- stats::qnorm(1 - alpha / 2)
   se_test_n1 <- if (test == "wls") est$se_beta_model_n1 else NULL
   power_values <- sapply(sample_sizes, function(n_per_stratum) {
      .rmst_wald_power(est$beta_effect, est$se_beta_n1,
                       n_per_stratum * est$n_strata, z_alpha, se_test_n1)
   })

   results_df <- data.frame(N_per_Stratum = sample_sizes, Power = power_values)

   # --- 3. Plot and Return ---
   p <- .rmst_power_curve_plot(
      results_df, "N_per_Stratum", "#E69F00",
      title = "Analytic Power Curve: Multiplicative Stratified RMST Model",
      subtitle = if (test == "wls")
         "Weighted least-squares t-test; sandwich sampling variance."
      else "Log-linear approximation with robust sandwich variance (Wang et al. 2019).",
      xlab = "Sample Size Per Stratum")

   return(list(results_data = results_df, results_plot = p,
               results_summary = NULL,
               model_output = .ms_model_output(est, arm_var, test)))
}

# Sample Size Search ------------------------------------------------------

#' @title Find Sample Size for a Multiplicative Stratified RMST Model (Analytic)
#' @description Calculates the required sample size for a target power using the analytic
#'   (approximate) method from Wang et al. (2019).
#'
#' @details
#' This function estimates the log-RMST ratio and its asymptotic variance once
#' from the pilot data, then increases the per-stratum sample size until the
#' target power is reached or the search limit is hit. It uses the same
#' log-linear approximation and robust sandwich variance as
#' `MS.power.analytical`, including the uncapped IPCW weights and the `test`
#' argument described there.
#'
#' @param pilot_data A `data.frame` containing pilot study data.
#' @param time_var A character string for the time-to-event variable.
#' @param status_var A character string for the event status variable (1=event, 0=censored).
#' @param arm_var A character string for the treatment arm variable (1=treatment, 0=control).
#' @param strata_var A character string for the stratification variable.
#' @param target_power A single numeric value for the desired power.
#' @param linear_terms An optional character vector of other covariate names.
#' @param L The numeric value for the RMST truncation time.
#' @param alpha The significance level (Type I error rate).
#' @param n_start The starting sample size *per stratum* for the search.
#' @param n_step The increment in sample size at each step of the search.
#' @param max_n_per_arm The maximum sample size *per stratum* to search up to.
#' @param verbose Logical; if \code{TRUE}, emit progress messages. Default \code{FALSE}.
#' @param test Which test the required sample size refers to. `"sandwich"` (default)
#'   is the robust sandwich Wald test; `"wls"` is the weighted least-squares t-test.
#'   See `MS.power.analytical`.
#'
#' @return A `list` containing:
#' \item{results_data}{A `data.frame` with the target power and required sample size.}
#' \item{results_plot}{A `ggplot` object visualizing the search path.}
#' \item{results_summary}{A `data.frame` summarizing the estimated log(RMST Ratio).}
#' @export
#' @importFrom survival Surv coxph basehaz
#' @importFrom stats lm as.formula complete.cases na.omit sd quantile pnorm qnorm coef vcov predict model.matrix residuals
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
#'  target_power = 0.80, L = 10,
#'  n_start = 100, n_step = 50
#' )
#' print(ss_results$results_data)
MS.ss.analytical <- function(pilot_data, time_var, status_var, arm_var, strata_var,
                             target_power, linear_terms = NULL, L, alpha = 0.05,
                             n_start = 50, n_step = 25, max_n_per_arm = 2000,
                             verbose = FALSE, test = c("sandwich", "wls")) {
   test <- match.arg(test)

   # --- 1. Estimate Parameters and Variance from Pilot Data (One Time) ---
   .rmst_verbose_message(verbose, "--- Estimating parameters from pilot data (log-linear approximation)... ---")
   est <- .estimate_ms_params(pilot_data, time_var, status_var, arm_var,
                              strata_var, linear_terms, L)

   # --- 2. Iterative Search for Sample Size ---
   .rmst_verbose_message(verbose, "--- Searching for Sample Size (Method: Analytic/Approximation) ---")
   search <- .rmst_analytic_ss_search(est$beta_effect, est$se_beta_n1, est$n_strata,
                                      target_power, alpha, n_start, n_step,
                                      max_n_per_arm, "/stratum", verbose,
                                      if (test == "wls") est$se_beta_model_n1 else NULL)
   final_n <- search$final_n

   # --- 3. Finalize and Return Results ---
   results_summary <- data.frame(Statistic = "Assumed log(RMST Ratio) (from pilot)", Value = est$beta_effect)
   results_df <- data.frame(Target_Power = target_power, Required_N_per_Stratum = final_n)
   search_path_df <- search$search_path_df
   names(search_path_df) <- c("N_per_Stratum", "Power")

   p <- .rmst_ss_search_plot(
      search_path_df, "N_per_Stratum", final_n, target_power,
      title = "Analytic Sample Size Search: Multiplicative Stratified RMST Model",
      subtitle = "Power calculated from formula at each step (approximate method).",
      xlab = "Sample Size Per Stratum")

   .rmst_verbose_message(verbose, "Calculation summary available in returned results_data.")

   return(list(results_data = results_df, results_plot = p,
               results_summary = results_summary,
               model_output = .ms_model_output(est, arm_var, test)))
}

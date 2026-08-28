#' @title Analyze Power for a Stratified Additive RMST Model (Analytic)
#' @description Performs power analysis for a stratified, additive RMST model using the
#'   analytic variance estimator based on the method of Zhang & Schaubel (2024).
#'
#' @details
#' This function computes power for the semiparametric additive RMST model
#' \eqn{\mu_{ij} = \mu_{0j} + \beta'Z_i}, where `i` indexes subjects and `j`
#' indexes strata.
#'
#' The method uses Inverse Probability of Censoring Weighting (IPCW), with
#' weights \eqn{\hat{W}_{ij}\Delta_i^Y} derived from a stratified Cox model for
#' the censoring times fit on the original time scale. Here
#' \eqn{\Delta_i^Y = 1} if the event occurs before \eqn{L} or follow-up reaches
#' \eqn{L}. The regression coefficient \eqn{\hat{\beta}} is estimated by
#' centering the covariates and RMST values within each stratum and then solving
#' the resulting estimating equations in closed form.
#'
#' For numerical stability, the estimated IPCW weights are capped at their
#' 99th percentile before estimation. This is a practical stabilization step
#' that is not part of the Zhang & Schaubel (2024) theory; the cap value and
#' the fraction of weights affected are reported in
#' \code{model_output$censoring_weights}.
#'
#' Power is obtained from the asymptotic sandwich variance of \eqn{\hat{\beta}}.
#' This implementation uses the robust variance estimator
#' \eqn{A_n^{-1} B_n (A_n^{-1})'}, where \eqn{A_n} and \eqn{B_n} are empirical
#' estimates of the variance components.
#'
#' @references
#' Zhang, Y. and Schaubel, D. E. (2024). Semiparametric Additive Modeling of
#' the Restricted Mean Survival Time. \emph{Biometrical Journal}, 66:e202200371.
#' \doi{10.1002/bimj.202200371}
#'
#' @param pilot_data A `data.frame` containing pilot study data.
#' @param time_var A character string for the time-to-event variable.
#' @param status_var A character string for the event status variable (1=event, 0=censored).
#' @param arm_var A character string for the treatment arm variable (1=treatment, 0=control).
#' @param strata_var A character string for the stratification variable.
#' @param sample_sizes A numeric vector of sample sizes *per stratum* to calculate power for.
#' @param linear_terms An optional character vector of other covariate names.
#' @param L The numeric value for the RMST truncation time.
#' @param alpha The significance level (Type I error rate).
#' @param verbose Logical; if \code{TRUE}, emit progress messages. Default \code{FALSE}.
#'
#' @return A `list` containing:
#' \item{results_data}{A `data.frame` with specified sample sizes and corresponding powers.}
#' \item{results_plot}{A `ggplot` object visualizing the power curve.}
#' @export
#' @importFrom stats as.formula complete.cases lm na.omit pnorm qnorm quantile sd vcov weighted.mean predict
#' @importFrom survival Surv coxph basehaz
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_minimal ylim
#' @importFrom dplyr group_by summarise across all_of mutate select left_join
#' @importFrom magrittr %>%
#'
#' @examples
#' set.seed(123)
#' pilot_df_strat <- data.frame(
#'  time = rexp(150, 0.1),
#'  status = rbinom(150, 1, 0.8),
#'  arm = rep(0:1, each = 75),
#'  region = factor(rep(c("A", "B", "C"), each = 50)),
#'  age = rnorm(150, 60, 10)
#' )
#' # Introduce an additive treatment effect
#' pilot_df_strat$time[pilot_df_strat$arm == 1] <-
#'   pilot_df_strat$time[pilot_df_strat$arm == 1] + 1.5
#'
#' power_results <- additive.power.analytical(
#'   pilot_data = pilot_df_strat,
#'   time_var = "time", status_var = "status", arm_var = "arm", strata_var = "region",
#'   sample_sizes = c(100, 150, 200),
#'   linear_terms = "age",
#'   L = 12
#' )
#' print(power_results$results_data)
#' print(power_results$results_plot)
#'
additive.power.analytical <- function(pilot_data, time_var, status_var, arm_var, strata_var,
                                      sample_sizes, linear_terms = NULL, L, alpha = 0.05,
                                      verbose = FALSE) {

   # --- 1. Estimate model parameters and sandwich variance from pilot data ---
   .rmst_verbose_message(verbose, "--- Estimating parameters from pilot data... ---")
   .rmst_verbose_message(verbose, "--- Estimating additive effect via stratum-centering... ---")
   est <- .estimate_additive_params(pilot_data, time_var, status_var, arm_var,
                                    strata_var, linear_terms, L)
   .rmst_verbose_message(verbose, "--- Calculating asymptotic variance... ---")

   # --- 2. Calculate Power ---
   .rmst_verbose_message(verbose, "--- Calculating power for specified sample sizes... ---")
   z_alpha <- stats::qnorm(1 - alpha / 2)
   power_values <- sapply(sample_sizes, function(n_per_stratum) {
      .rmst_wald_power(est$beta_effect, est$se_beta_n1,
                       n_per_stratum * est$n_strata, z_alpha)
   })

   results_df <- data.frame(N_per_Stratum = sample_sizes, Power = power_values)

   # --- 3. Plot and Return ---
   p <- .rmst_power_curve_plot(
      results_df, "N_per_Stratum", "#0072B2",
      title = "Analytic Power Curve: Additive Stratified RMST Model",
      subtitle = "Based on stratum-centered estimating equations.",
      xlab = "Sample Size Per Stratum")

   return(list(results_data = results_df, results_plot = p,
               results_summary = NULL,
               model_output = .additive_model_output(est, arm_var, strata_var)))
}


#' @title Find Sample Size for a Stratified Additive RMST Model (Analytic)
#' @description Calculates the required sample size for a target power using the analytic
#'   method for a stratified, additive RMST model.
#'
#' @details
#' This function estimates the additive treatment effect and its asymptotic
#' variance once from the pilot data, then increases the per-stratum sample
#' size until the target power is reached or the search limit is hit. It uses
#' the same stratum-centering framework as `additive.power.analytical`,
#' including the 99th-percentile IPCW weight cap described there.
#'
#' @references
#' Zhang, Y. and Schaubel, D. E. (2024). Semiparametric Additive Modeling of
#' the Restricted Mean Survival Time. \emph{Biometrical Journal}, 66:e202200371.
#' \doi{10.1002/bimj.202200371}
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
#'
#' @return A `list` containing:
#' \item{results_data}{A `data.frame` with the target power and required sample size.}
#' \item{results_plot}{A `ggplot` object visualizing the search path.}
#' \item{results_summary}{A `data.frame` summarizing the estimated treatment effect.}
#' @export
#' @importFrom stats as.formula complete.cases lm na.omit pnorm qnorm quantile sd vcov weighted.mean predict
#' @importFrom survival Surv coxph basehaz
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline geom_vline labs theme_minimal
#' @importFrom knitr kable
#' @importFrom dplyr group_by summarise across all_of mutate select left_join
#' @importFrom magrittr %>%
#' @examples
#' set.seed(123)
#' pilot_df_strat <- data.frame(
#'  time = rexp(150, 0.1),
#'  status = rbinom(150, 1, 0.8),
#'  arm = rep(0:1, each = 75),
#'  region = factor(rep(c("A", "B", "C"), each = 50)),
#'  age = rnorm(150, 60, 10)
#' )
#' # Introduce an additive treatment effect
#' pilot_df_strat$time[pilot_df_strat$arm == 1] <-
#'   pilot_df_strat$time[pilot_df_strat$arm == 1] + 1.5
#'
#'   # Find the required sample size per stratum for 80% power
#'   ss_results <- additive.ss.analytical(
#'     pilot_data = pilot_df_strat,
#'     time_var = "time", status_var = "status",
#'     arm_var = "arm", strata_var = "region",
#'     target_power = 0.50,
#'     L = 18, #
#'     n_start = 200,
#'     n_step = 50,
#'     max_n_per_arm = 1000
#'   )
#'   print(ss_results$results_data)
#'   print(ss_results$results_plot)
additive.ss.analytical <- function(pilot_data, time_var, status_var, arm_var, strata_var,
                                   target_power, linear_terms = NULL, L, alpha = 0.05,
                                   n_start = 50, n_step = 25, max_n_per_arm = 2000,
                                   verbose = FALSE) {

   # --- 1. Estimate Parameters and Variance from Pilot Data (One Time) ---
   .rmst_verbose_message(verbose, "--- Estimating parameters from pilot data for analytic search... ---")
   est <- .estimate_additive_params(pilot_data, time_var, status_var, arm_var,
                                    strata_var, linear_terms, L)

   # --- 2. Iterative Search for Sample Size ---
   .rmst_verbose_message(verbose, "--- Searching for Sample Size (Method: Additive Analytic) ---")
   search <- .rmst_analytic_ss_search(est$beta_effect, est$se_beta_n1, est$n_strata,
                                      target_power, alpha, n_start, n_step,
                                      max_n_per_arm, "/stratum", verbose)
   final_n <- search$final_n

   # --- 3. Finalize and Return Results ---
   results_summary <- data.frame(Statistic = "Assumed RMST Difference (from pilot)", Value = est$beta_effect)
   results_df <- data.frame(Target_Power = target_power, Required_N_per_Stratum = final_n)
   search_path_df <- search$search_path_df
   names(search_path_df) <- c("N_per_Stratum", "Power")

   p <- .rmst_ss_search_plot(
      search_path_df, "N_per_Stratum", final_n, target_power,
      title = "Analytic Sample Size Search: Additive Stratified RMST Model",
      subtitle = "Power calculated from formula at each step.",
      xlab = "Sample Size Per Stratum")

   .rmst_verbose_message(verbose, "Calculation summary available in returned results_data.")

   return(list(results_data = results_df, results_plot = p,
               results_summary = results_summary,
               model_output = .additive_model_output(est, arm_var, strata_var)))
}

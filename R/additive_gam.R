#' @title Calculate Power for a Semiparametric Additive RMST Model
#' @description Performs a power analysis for given sample sizes using a flexible,
#'   semiparametric additive model for the RMST based on pseudo-observations.
#'
#' @note `status_var` should be `1` for an event, `0` for censored. `arm_var`
#'   should be `1` for treatment, `0` for control.
#'
#' @param pilot_data A data.frame with pilot study data.
#' @param time_var,status_var,arm_var Strings for column names.
#' @param sample_sizes A numeric vector of sample sizes per arm/stratum.
#' @param strata_var An optional string for a stratification variable.
#' @param linear_terms Optional character vector of covariates with a linear effect.
#' @param smooth_terms Optional character vector of covariates with a non-linear effect.
#' @param tau The truncation time for RMST.
#' @param n_sim Number of bootstrap simulations.
#' @param alpha The significance level.
#' @param parallel.cores Number of cores for parallel processing. Default is 1 (sequential).
#'
#' @return A list containing results data, a ggplot object, and a summary data.frame.
#' @references This function uses a pseudo-observation approach. For an alternative
#'   IPCW-based method, see Zhang, Y., & Schaubel, D. E. (2024), Biometrical Journal.
#'
#' @importFrom survival survfit Surv
#' @importFrom stats stepfun integrate lm as.formula complete.cases na.omit sd quantile proc.time
#' @importFrom mgcv gam summary.gam
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline labs theme_minimal ylim
#' @importFrom knitr kable
#' @importFrom foreach foreach %dopar%
#'
#' @examples
#' # Create a small pilot dataset for the example
#' pilot_data_gam <- data.frame(
#'   surv_time = stats::rexp(100, rate = 0.1),
#'   censor_status = stats::rbinom(100, 1, 0.7),
#'   arm = rep(0:1, each = 50),
#'   age = stats::rnorm(100, 50, 10),
#'   biomarker = stats::rnorm(100, 10, 3)
#' )
#'
#' # --- Power Calculation Example ---
#' results_power <- design_rmst_gam_power(
#'   pilot_data = pilot_data_gam,
#'   time_var = "surv_time",
#'   status_var = "censor_status",
#'   arm_var = "arm",
#'   sample_sizes = c(100, 150),
#'   linear_terms = c("age"),
#'   smooth_terms = c("biomarker"),
#'   tau = 10,
#'   n_sim = 10 # Low n_sim for example speed
#' )
#' @export
design_rmst_gam_power <- function(pilot_data, time_var, status_var, arm_var, strata_var = NULL,
                                  sample_sizes, linear_terms = NULL, smooth_terms = NULL,
                                  tau, n_sim = 1000, alpha = 0.05,
                                  parallel.cores = 1) {
   start_time <- proc.time()
   if (is.null(sample_sizes)) stop("You must provide a numeric vector for 'sample_sizes'.")
   if (parallel.cores > 1) {
      if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
         stop("Packages 'future' and 'future.apply' are required for parallel processing.")
      }
   }

   is_stratified <- !is.null(strata_var)

   # --- Prepare data ---
   core_vars <- c(time_var, status_var, arm_var, strata_var)
   if (is.null(linear_terms) && is.null(smooth_terms)) {
      linear_terms <- setdiff(names(pilot_data), core_vars)
      if(length(linear_terms) > 0) message(paste("Using unspecified columns as linear terms:", paste(linear_terms, collapse=", ")))
   }
   all_vars <- c(core_vars, linear_terms, smooth_terms)
   pilot_data <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]

   # --- Model Formula Construction ---
   smooth_part <- if (!is.null(smooth_terms)) paste0("s(", smooth_terms, ")") else NULL
   if (is_stratified) {
      pilot_groups <- split(pilot_data, pilot_data[[strata_var]])
      all_terms <- c(strata_var, paste0(arm_var, ":", strata_var), smooth_part, linear_terms)
      test_term_pattern <- paste0(arm_var, "1:", strata_var)
   } else {
      pilot_groups <- list(pilot_data)
      all_terms <- c(arm_var, smooth_part, linear_terms)
      test_term_pattern <- arm_var
   }
   model_rhs <- paste(all_terms[!sapply(all_terms, is.null)], collapse = " + ")
   model_formula <- stats::as.formula(paste("pseudo_obs ~", model_rhs))

   # --- Helper to calculate jackknife pseudo-observations ---
   get_pseudo_obs <- function(time, status, tau) {
      n <- length(time)
      if (n == 0) return(numeric(0))
      km_fit_full <- survival::survfit(survival::Surv(time, status) ~ 1)
      km_step_full <- stats::stepfun(km_fit_full$time, c(1, km_fit_full$surv))
      rmst_full <- tryCatch(stats::integrate(km_step_full, 0, tau, subdivisions=2000, stop.on.error = FALSE)$value, error = function(e) 0)
      rmst_loo <- sapply(1:n, function(i) {
         if(n > 1) {
            km_fit_loo <- survival::survfit(survival::Surv(time[-i], status[-i]) ~ 1)
            km_step_loo <- stats::stepfun(km_fit_loo$time, c(1, km_fit_loo$surv))
            tryCatch(stats::integrate(km_step_loo, 0, tau, subdivisions=2000, stop.on.error = FALSE)$value, error = function(e) 0)
         } else { 0 }
      })
      return(n * rmst_full - (n - 1) * rmst_loo)
   }

   group_label <- if(is_stratified) "/stratum" else "/arm"
   cat("--- Calculating Power (Method: Additive GAM for RMST) ---\n")
   message("Model: pseudo_obs ~ ", model_rhs)

   all_sim_outputs <- vector("list", length(sample_sizes))

   # Set up parallel plan
   if (parallel.cores > 1) {
      future::plan(future::multisession, workers = parallel.cores)
   } else {
      future::plan(future::sequential)
   }
   on.exit(future::plan(future::sequential), add = TRUE) # Ensure plan is reset

   for (i in seq_along(sample_sizes)) {
      n_per_group <- sample_sizes[i]
      cat("Simulating for n =", n_per_group, group_label, "...\n")

      sim_results_list <- future.apply::future_lapply(1:n_sim, function(j) {
         p_val <- NA_real_; est_val <- NA_real_; se_val <- NA_real_
         if (is_stratified) {
            boot_list <- lapply(pilot_groups, function(group_df) group_df[sample(1:nrow(group_df), size = n_per_group, replace = TRUE), ])
         } else {
            group_arms <- split(pilot_groups[[1]], pilot_groups[[1]][[arm_var]])
            boot_list <- lapply(group_arms, function(arm_df) arm_df[sample(1:nrow(arm_df), size = n_per_group, replace = TRUE), ])
         }
         boot_data <- do.call(rbind, boot_list)
         pseudo_obs_list <- by(boot_data, boot_data[[if(is_stratified) strata_var else arm_var]], function(sub_data) {
            sub_data$pseudo_obs <- get_pseudo_obs(sub_data[[time_var]], sub_data[[status_var]], tau)
            sub_data
         })
         boot_data <- do.call(rbind, pseudo_obs_list)
         fit <- tryCatch(mgcv::gam(model_formula, data = boot_data), error = function(e) NULL)
         if (!is.null(fit)) {
            sfit <- mgcv::summary.gam(fit)
            p_table <- sfit$p.table
            matching_rows <- grep(test_term_pattern, rownames(p_table), fixed = TRUE)
            if (length(matching_rows) > 0) {
               p_val <- min(p_table[matching_rows, "Pr(>|t|)"], na.rm = TRUE)
               est_val <- mean(p_table[matching_rows, "Estimate"], na.rm = TRUE)
               se_val <- mean(p_table[matching_rows, "Std. Error"], na.rm = TRUE)
            }
         }
         return(list(p_value = p_val, estimate = est_val, std_error = se_val))
      }, future.seed = TRUE)

      p_values <- vapply(sim_results_list, `[[`, "p_value", FUN.VALUE = numeric(1))
      estimates <- vapply(sim_results_list, `[[`, "estimate", FUN.VALUE = numeric(1))
      std_errors <- vapply(sim_results_list, `[[`, "std_error", FUN.VALUE = numeric(1))

      all_sim_outputs[[i]] <- list(power = mean(p_values < alpha, na.rm = TRUE),
                                   estimates = estimates,
                                   std_errors = std_errors)
   }

   end_time <- proc.time()



   power_values <- sapply(all_sim_outputs, `[[`, "power")
   results_df <- data.frame(N_per_Group = sample_sizes, Power = power_values)
   best <- which.max(sample_sizes)
   est <- stats::na.omit(all_sim_outputs[[best]]$estimates)
   se  <- stats::na.omit(all_sim_outputs[[best]]$std_errors)
   results_summary <- NULL
   if (length(est) > 1) {
      results_summary <- data.frame(
         Statistic = c("Mean RMST Difference", "Mean Standard Error", "95% CI Lower", "95% CI Upper"),
         Value = c(mean(est), mean(se, na.rm=TRUE), mean(est) - 1.96 * stats::sd(est), mean(est) + 1.96 * stats::sd(est))
      )
   }

   p <- ggplot2::ggplot(results_df, ggplot2::aes(x = N_per_Group, y = Power)) +
      ggplot2::geom_line(color = "#0072B2", linewidth = 1) +
      ggplot2::geom_point(color = "#0072B2", size = 3) +
      ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
      ggplot2::labs(title = "Power Curve: Additive GAM RMST Model",
                    x = paste0("Sample Size Per Arm", if(is_stratified) "/Stratum" else ""), y = "Estimated Power") +
      ggplot2::ylim(0, 1) + ggplot2::theme_minimal()

   elapsed_time <- round((end_time - start_time)["elapsed"], 2)
   message(paste("Total simulation time:", elapsed_time, "seconds"))
   cat("\n--- Simulation Summary ---\n")
   if (!is.null(results_summary)) {
      print(knitr::kable(results_summary, caption = "Estimated Treatment Effect (RMST Difference)"))
   } else {
      cat("No valid estimates were generated to create a summary.\n")
   }
   print(p)

   return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}



#' @title Find Sample Size for a Semiparametric Additive RMST Model
#' @description Performs an iterative sample size search to achieve a target power
#'   using a flexible, semiparametric additive model for the RMST.
#'
#' @note `status_var` should be `1` for an event, `0` for censored. `arm_var`
#'   should be `1` for treatment, `0` for control.
#'
#' @param pilot_data A data.frame with pilot study data.
#' @param time_var,status_var,arm_var Strings for column names.
#' @param target_power A single numeric value for the target power.
#' @param strata_var An optional string for a stratification variable.
#' @param linear_terms Optional character vector of covariates with a linear effect.
#' @param smooth_terms Optional character vector of covariates with a non-linear effect.
#' @param tau The truncation time for RMST.
#' @param n_sim Number of bootstrap simulations per search step.
#' @param alpha The significance level.
#' @param parallel.cores Number of cores for parallel processing. Default is 1 (sequential).
#' @param patience Number of consecutive non-improving steps in the iterative
#'   search before terminating. Default is 5.
#' @param n_start The starting sample size per arm/stratum for the search.
#' @param n_step The increment in sample size at each step of the search.
#' @param max_n_per_arm The maximum sample size per arm/stratum to search up to.
#'
#' @return A list containing results data, a ggplot object of the search path,
#'   and a summary data.frame.
#' @references This function uses a pseudo-observation approach. For an alternative
#'   IPCW-based method, see Zhang, Y., & Schaubel, D. E. (2024), Biometrical Journal.
#'
#' @importFrom survival survfit Surv
#' @importFrom stats stepfun integrate lm as.formula complete.cases na.omit sd quantile proc.time
#' @importFrom mgcv gam summary.gam
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline geom_vline labs theme_minimal
#' @importFrom knitr kable
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#'
#' @examples
#' # Create a small pilot dataset for the example
#' pilot_data_gam <- data.frame(
#'   surv_time = c(stats::rexp(50, 0.1), stats::rexp(50, 0.05)), #<-- effect
#'   censor_status = stats::rbinom(100, 1, 0.8),
#'   arm = rep(0:1, each = 50),
#'   age = stats::rnorm(100, 50, 10)
#' )
#'
#' # --- Sample Size Calculation Example ---
#' results_n <- design_rmst_gam_ss(
#'   pilot_data = pilot_data_gam,
#'   time_var = "surv_time",
#'   status_var = "censor_status",
#'   arm_var = "arm",
#'   target_power = 0.80,
#'   linear_terms = c("age"),
#'   tau = 10,
#'   n_sim = 20, # Low n_sim for example speed
#'   n_step = 100,
#'   patience = 2
#' )
#' @export
design_rmst_gam_ss <- function(pilot_data, time_var, status_var, arm_var, strata_var = NULL,
                               target_power,
                               linear_terms = NULL, smooth_terms = NULL,
                               tau, n_sim = 1000, alpha = 0.05,
                               parallel.cores = 1, patience = 5,
                               n_start = 50, n_step = 25, max_n_per_arm = 2000) {

   start_time <- proc.time()
   if (is.null(target_power) || length(target_power) != 1 || !is.numeric(target_power)) {
      stop("You must provide a single numeric value for 'target_power'.")
   }
   if (parallel.cores > 1) {
      if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
         stop("Packages 'future' and 'future.apply' are required for parallel processing.")
      }
   }

   is_stratified <- !is.null(strata_var)

   # --- Prepare data ---
   core_vars <- c(time_var, status_var, arm_var, strata_var)
   if (is.null(linear_terms) && is.null(smooth_terms)) {
      linear_terms <- setdiff(names(pilot_data), core_vars)
      if(length(linear_terms) > 0) message(paste("Using unspecified columns as linear terms:", paste(linear_terms, collapse=", ")))
   }
   all_vars <- c(core_vars, linear_terms, smooth_terms)
   pilot_data <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]

   # --- Model Formula Construction ---
   smooth_part <- if (!is.null(smooth_terms)) paste0("s(", smooth_terms, ")") else NULL
   if (is_stratified) {
      pilot_groups <- split(pilot_data, pilot_data[[strata_var]])
      all_terms <- c(strata_var, paste0(arm_var, ":", strata_var), smooth_part, linear_terms)
      test_term_pattern <- paste0(arm_var, "1:", strata_var)
   } else {
      pilot_groups <- list(pilot_data)
      all_terms <- c(arm_var, smooth_part, linear_terms)
      test_term_pattern <- arm_var
   }
   model_rhs <- paste(all_terms[!sapply(all_terms, is.null)], collapse = " + ")
   model_formula <- stats::as.formula(paste("pseudo_obs ~", model_rhs))

   # --- Helper to calculate jackknife pseudo-observations ---
   get_pseudo_obs <- function(time, status, tau) {
      n <- length(time)
      if (n == 0) return(numeric(0))
      km_fit_full <- survival::survfit(survival::Surv(time, status) ~ 1)
      km_step_full <- stats::stepfun(km_fit_full$time, c(1, km_fit_full$surv))
      rmst_full <- tryCatch(stats::integrate(km_step_full, 0, tau, subdivisions=2000, stop.on.error = FALSE)$value, error = function(e) 0)
      rmst_loo <- sapply(1:n, function(i) {
         if(n > 1) {
            km_fit_loo <- survival::survfit(survival::Surv(time[-i], status[-i]) ~ 1)
            km_step_loo <- stats::stepfun(km_fit_loo$time, c(1, km_fit_loo$surv))
            tryCatch(stats::integrate(km_step_loo, 0, tau, subdivisions=2000, stop.on.error = FALSE)$value, error = function(e) 0)
         } else { 0 }
      })
      return(n * rmst_full - (n - 1) * rmst_loo)
   }

   group_label <- if(is_stratified) "/stratum" else "/arm"
   cat("--- Searching for Sample Size (Method: Additive GAM for RMST) ---\n")
   message("Model: pseudo_obs ~ ", model_rhs)

   # Set up parallel plan
   if (parallel.cores > 1) {
      future::plan(future::multisession, workers = parallel.cores)
   } else {
      future::plan(future::sequential)
   }
   on.exit(future::plan(future::sequential), add = TRUE) # Ensure plan is reset

   # --- Main Search Loop ---
   current_n <- n_start
   max_power_so_far <- -1
   stagnation_counter <- 0
   n_at_max_power <- n_start
   best_sim_output <- NULL
   search_results <- list()
   final_n <- NA

   while (current_n <= max_n_per_arm) {
      cat(paste0("  N = ", current_n, group_label, ", Calculating Power..."))

      sim_results_list <- future.apply::future_lapply(1:n_sim, function(j) {
         p_val <- NA_real_; est_val <- NA_real_; se_val <- NA_real_
         if (is_stratified) {
            boot_list <- lapply(pilot_groups, function(group_df) group_df[sample(1:nrow(group_df), size = current_n, replace = TRUE), ])
         } else {
            group_arms <- split(pilot_groups[[1]], pilot_groups[[1]][[arm_var]])
            boot_list <- lapply(group_arms, function(arm_df) arm_df[sample(1:nrow(arm_df), size = current_n, replace = TRUE), ])
         }
         boot_data <- do.call(rbind, boot_list)
         pseudo_obs_list <- by(boot_data, boot_data[[if(is_stratified) strata_var else arm_var]], function(sub_data) {
            sub_data$pseudo_obs <- get_pseudo_obs(sub_data[[time_var]], sub_data[[status_var]], tau)
            sub_data
         })
         boot_data <- do.call(rbind, pseudo_obs_list)
         fit <- tryCatch(mgcv::gam(model_formula, data = boot_data), error = function(e) NULL)
         if (!is.null(fit)) {
            sfit <- mgcv::summary.gam(fit)
            p_table <- sfit$p.table
            matching_rows <- grep(test_term_pattern, rownames(p_table), fixed = TRUE)
            if (length(matching_rows) > 0) {
               p_val <- min(p_table[matching_rows, "Pr(>|t|)"], na.rm = TRUE)
               est_val <- mean(p_table[matching_rows, "Estimate"], na.rm = TRUE)
               se_val <- mean(p_table[matching_rows, "Std. Error"], na.rm = TRUE)
            }
         }
         return(list(p_value = p_val, estimate = est_val, std_error = se_val))
      }, future.seed = TRUE)

      p_values <- vapply(sim_results_list, `[[`, "p_value", FUN.VALUE = numeric(1))
      estimates <- vapply(sim_results_list, `[[`, "estimate", FUN.VALUE = numeric(1))
      std_errors <- vapply(sim_results_list, `[[`, "std_error", FUN.VALUE = numeric(1))
      sim_output <- list(power = mean(p_values < alpha, na.rm = TRUE),
                         estimates = estimates, std_errors = std_errors)

      calculated_power <- if(is.finite(sim_output$power)) sim_output$power else 0
      search_results[[as.character(current_n)]] <- calculated_power
      cat(paste0(" Power = ", round(calculated_power, 3), "\n"))

      if (calculated_power >= target_power) {
         message("Success: Target power reached at N = ", current_n, group_label, ".")
         best_sim_output <- sim_output
         final_n <- current_n
         break
      }
      if (calculated_power > max_power_so_far) {
         max_power_so_far <- calculated_power
         n_at_max_power <- current_n
         best_sim_output <- sim_output
         stagnation_counter <- 0
      } else {
         stagnation_counter <- stagnation_counter + 1
      }
      if (stagnation_counter >= patience) {
         warning(paste0("Search terminated due to stagnation. Best N = ", n_at_max_power,
                        " with power = ", round(max_power_so_far, 3)), call. = FALSE)
         final_n <- n_at_max_power
         break
      }
      current_n <- current_n + n_step
   }

   if (is.na(final_n) && current_n > max_n_per_arm) {
      warning(paste("Target power", target_power, "not achieved by max N of", max_n_per_arm))
      final_n <- n_at_max_power
   }

   # --- Finalize Summary, Plot, and Results ---
   results_summary <- NULL
   if (!is.null(best_sim_output)) {
      est <- stats::na.omit(best_sim_output$estimates)
      se  <- stats::na.omit(best_sim_output$std_errors)
      if(length(est) > 1) {
         results_summary <- data.frame(
            Statistic = c("Mean RMST Difference", "Mean Standard Error", "95% CI Lower", "95% CI Upper"),
            Value = c(mean(est), mean(se, na.rm=TRUE), mean(est) - 1.96 * stats::sd(est), mean(est) + 1.96 * stats::sd(est))
         )
      }
   }

   results_df <- data.frame(Target_Power = target_power, Required_N_per_Group = final_n)
   search_path_df <- data.frame(N_per_Group = as.integer(names(search_results)),
                                Power = unlist(search_results))

   p <- ggplot2::ggplot(na.omit(search_path_df), ggplot2::aes(x = N_per_Group, y = Power)) +
      ggplot2::geom_line(color = "#009E73", linewidth = 1) +
      ggplot2::geom_point(color = "#009E73", size = 3) +
      ggplot2::geom_hline(yintercept = target_power, linetype = "dashed", color = "red") +
      ggplot2::labs(title = "Sample Size Search Path: Additive GAM RMST Model",
                    x = paste0("Sample Size Per Arm", if(is_stratified) "/Stratum" else ""),
                    y = "Calculated Power") +
      ggplot2::theme_minimal()

   end_time <- proc.time()
   elapsed_time <- round((end_time - start_time)["elapsed"], 2)
   message(paste("Total simulation time:", elapsed_time, "seconds"))

   cat("\n--- Simulation Summary ---\n")
   if (!is.null(results_summary)) {
      print(knitr::kable(results_summary, caption = "Estimated Treatment Effect (RMST Difference)"))
   } else {
      cat("No valid estimates were generated to create a summary.\n")
   }

   print(p)

   return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}

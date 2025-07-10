# Power Calculation -------------------------------------------------------

#' @title Analyze Power for a Linear RMST Model via Simulation
#' @description Performs a power analysis for given sample sizes based on the direct
#'   linear regression model for RMST, using a bootstrap simulation approach.
#'
#' @details
#' This function estimates power by generating a number of bootstrap
#' samples (`n_sim`) from the provided pilot data by resampling with replacement.
#' For each bootstrap sample, it performs the following steps:
#' 1.  Estimates the censoring distribution using the Kaplan-Meier method (`survival::survfit`).
#' 2.  Calculates Inverse Probability of Censoring Weights (IPCW) for each observation.
#' 3.  Fits a weighted linear model (`stats::lm`) to the RMST of the uncensored subjects.
#' 4.  Extracts the p-value for the treatment `arm_var` coefficient.
#'
#' The final power for a given sample size is the proportion of the `n_sim` simulations
#' where this p-value is less than the significance level `alpha`. This simulation-based
#' approach is robust but can be computationally intensive.
#'
#' @note `status_var` should be `1` for an event, `0` for censored. `arm_var`
#'   should be `1` for treatment, `0` for control.
#'
#' @param pilot_data A `data.frame` with pilot study data.
#' @param time_var A character string for the time-to-event variable.
#' @param status_var A character string for the event status variable.
#' @param arm_var A character string for the treatment arm variable.
#' @param sample_sizes A numeric vector of sample sizes *per arm* to calculate power for.
#' @param linear_terms Optional character vector of other covariates for the linear model.
#' @param L The numeric truncation time for RMST.
#' @param n_sim The number of bootstrap simulations to run for each sample size.
#' @param alpha The significance level (Type I error rate).
#'
#' @return A `list` containing:
#' \item{results_data}{A `data.frame` of sample sizes and corresponding estimated power.}
#' \item{results_plot}{A `ggplot` object visualizing the power curve.}
#' \item{results_summary}{A `data.frame` with summary statistics for the estimated treatment effect from the largest sample size simulation.}
#'
#' @importFrom survival survfit Surv
#' @importFrom stats lm as.formula complete.cases na.omit sd quantile
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_minimal ylim
#' @importFrom knitr kable
#' @export
#' @examples
#' \dontrun{
#' pilot_df <- data.frame(
#'   time = rexp(100, 0.1),
#'   status = rbinom(100, 1, 0.7),
#'   arm = rep(0:1, each = 50),
#'   age = rnorm(100, 60, 8)
#' )
#' # Introduce a treatment effect for a more interesting example
#' pilot_df$time[pilot_df$arm == 1] <- pilot_df$time[pilot_df$arm == 1] * 1.5
#'
#' power_results <- linear.power.boot(
#'   pilot_data = pilot_df,
#'   time_var = "time",
#'   status_var = "status",
#'   arm_var = "arm",
#'   linear_terms = "age",
#'   sample_sizes = c(100, 150, 200),
#'   L = 10,
#'   n_sim = 200 # Use more simulations in practice (e.g., 1000)
#' )
#' print(power_results$results_data)
#' print(power_results$results_plot)
#' }
linear.power.boot <- function(pilot_data, time_var, status_var, arm_var,
                              sample_sizes,linear_terms = NULL, L, n_sim = 1000, alpha = 0.05)
{

   start_time <- proc.time()
   if (is.null(sample_sizes)) stop("You must provide a numeric vector for 'sample_sizes'.")

   # --- Prepare data ---
   core_vars <- c(time_var, status_var, arm_var)
   all_vars <- c(core_vars, linear_terms)
   pilot_data <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]
   pilot_groups <- split(pilot_data, pilot_data[[arm_var]])

   model_rhs <- paste(c(arm_var, linear_terms), collapse = " + ")
   model_formula <- as.formula(paste("Y_rmst ~", model_rhs))

   # --- Run simulation ---
   cat("--- Calculating Power (Method: Linear RMST with IPCW) ---\n")
   message("Model: Y_rmst ~ ", model_rhs)
   results_summary <- NULL

   all_sim_outputs <- vector("list", length(sample_sizes))

   for (i in seq_along(sample_sizes)) {
      n_per_arm <- sample_sizes[i]
      cat("Simulating for n =", n_per_arm, "per arm...\n")

      p_values <- rep(NA_real_, n_sim)
      estimates <- rep(NA_real_, n_sim)
      std_errors <- rep(NA_real_, n_sim)

      for (j in seq_len(n_sim)) {
         boot_list <- lapply(pilot_groups, function(df) df[sample(seq_len(nrow(df)), size = n_per_arm, replace = TRUE), ])
         boot_data <- do.call(rbind, boot_list)
         boot_data[[arm_var]] <- factor(boot_data[[arm_var]], levels = c(0, 1))

         is_censored <- boot_data[[status_var]] == 0
         cens_fit <- tryCatch(survival::survfit(Surv(boot_data[[time_var]], is_censored) ~ 1), error = function(e) NULL)
         if (is.null(cens_fit)) next
         surv_summary <- tryCatch(summary(cens_fit, times = pmin(boot_data[[time_var]], L), extend = TRUE), error = function(e) NULL)
         if (is.null(surv_summary)) next

         weights <- 1 / surv_summary$surv
         finite_weights <- weights[is.finite(weights)]
         if (length(finite_weights) > 0) {
            weight_cap <- stats::quantile(finite_weights, probs = 0.99, na.rm = TRUE)
            weights[weights > weight_cap] <- weight_cap
         }
         weights[!is.finite(weights)] <- NA

         boot_data$Y_rmst <- pmin(boot_data[[time_var]], L)
         fit_data <- boot_data[boot_data[[status_var]] == 1 & is.finite(weights), ]
         fit_weights <- weights[boot_data[[status_var]] == 1 & is.finite(weights)]

         if (nrow(fit_data) > (length(all_vars) + 1)) {
            fit <- tryCatch(lm(model_formula, data = fit_data, weights = fit_weights), error = function(e) NULL)
            if (!is.null(fit)) {
               sfit <- summary(fit)
               # Robustly find the coefficient name, e.g., 'arm1'
               test_term_pattern <- paste0("^", arm_var, "1$")
               test_term <- grep(test_term_pattern, rownames(sfit$coefficients), value = TRUE)[1]
               if (!is.na(test_term)) {
                  p_values[j] <- tryCatch(sfit$coefficients[test_term, "Pr(>|t|)"], error = function(e) NA)
                  estimates[j] <- tryCatch(sfit$coefficients[test_term, "Estimate"], error = function(e) NA)
                  std_errors[j] <- tryCatch(sfit$coefficients[test_term, "Std. Error"], error = function(e) NA)
               }
            }
         }
      }
      all_sim_outputs[[i]] <- list(power = mean(p_values < alpha, na.rm = TRUE),
                                   estimates = estimates,
                                   std_errors = std_errors)
   }

   power_values <- sapply(all_sim_outputs, `[[`, "power")
   results_df <- data.frame(N_per_Arm = sample_sizes, Power = power_values)
   best <- which.max(sample_sizes)
   est <- na.omit(all_sim_outputs[[best]]$estimates)
   se  <- na.omit(all_sim_outputs[[best]]$std_errors)

   if (length(est) > 1) {
      results_summary <- data.frame(
         Statistic = c("Mean RMST Difference", "Mean Standard Error", "95% CI Lower", "95% CI Upper"),
         Value = c(mean(est), mean(se, na.rm=TRUE), mean(est) - 1.96 * sd(est), mean(est) + 1.96 * sd(est))
      )
   }

   p <- ggplot2::ggplot(results_df, ggplot2::aes(x = N_per_Arm, y = Power)) +
      ggplot2::geom_line(color = "#D55E00", linewidth = 1) +
      ggplot2::geom_point(color = "#D55E00", size = 3) +
      ggplot2::labs(title = "Power Curve: Linear IPCW RMST Model",
                    x = "Sample Size Per Arm", y = "Estimated Power") +
      ggplot2::ylim(0, 1) + ggplot2::theme_minimal()

   end_time <- proc.time()
   elapsed_time <- round((end_time - start_time)["elapsed"], 2)
   message(paste("Total simulation time:", elapsed_time, "seconds"))

   cat("\n--- Simulation Summary ---\n")
   if (!is.null(results_summary)) {
      print(knitr::kable(results_summary, caption = "Estimated Treatment Effect (RMST Difference)"))
   } else {
      cat("No valid estimates were generated to create a summary.\n")
   }

   return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}

# Sample_Size_Search ------------------------------------------------------


#' @title Find Sample Size for a Linear RMST Model via Simulation
#' @description Performs an iterative sample size search to achieve a target power
#'   based on the direct linear regression model for RMST, using bootstrap simulation.
#'
#' @details
#' This function iteratively searches for the required sample size to achieve
#' a specified `target_power`. At each step of the search, it runs a full bootstrap
#' simulation (`n_sim` iterations), as described in `linear.power.boot`, to
#' estimate the power for the current sample size. The search stops when the
#' target power is achieved or other stopping criteria (e.g., `patience`) are met.
#' Due to the nested simulation structure, this function can be very time-consuming.
#'
#' @note `status_var` should be `1` for an event, `0` for censored. `arm_var`
#'   should be `1` for treatment, `0` for control.
#'
#' @param pilot_data A `data.frame` with pilot study data.
#' @param time_var A character string for the time-to-event variable.
#' @param status_var A character string for the event status variable.
#' @param arm_var A character string for the treatment arm variable.
#' @param target_power A single numeric value for the target power (e.g., 0.80).
#' @param linear_terms Optional character vector of other covariates for the linear model.
#' @param L The numeric truncation time for RMST.
#' @param n_sim The number of bootstrap simulations per search step.
#' @param alpha The significance level.
#' @param patience The number of consecutive non-improving steps in the search before terminating.
#' @param n_start The starting sample size *per arm* for the search.
#' @param n_step The increment in sample size at each step of the search.
#' @param max_n_per_arm The maximum sample size *per arm* to search up to.
#'
#' @return A `list` containing:
#' \item{results_data}{A `data.frame` with the target power and the final required sample size per arm.}
#' \item{results_plot}{A `ggplot` object showing the search path.}
#' \item{results_summary}{A `data.frame` with summary statistics for the estimated treatment effect from the final simulation.}
#'
#' @importFrom survival survfit Surv
#' @importFrom stats lm as.formula complete.cases na.omit sd quantile
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline geom_vline labs theme_minimal
#' @importFrom knitr kable
#' @export
#' @examples
#' \dontrun{
#' pilot_df_effect <- data.frame(
#'   time = c(rexp(50, 0.1), rexp(50, 0.05)), # Effect present
#'   status = rbinom(100, 1, 0.8),
#'   arm = rep(0:1, each = 50)
#' )
#' ss_results <- linear.ss.boot(
#'   pilot_data = pilot_df_effect,
#'   time_var = "time",
#'   status_var = "status",
#'   arm_var = "arm",
#'   target_power = 0.80,
#'   L = 10,
#'   n_sim = 200, # Low n_sim for example
#'   patience = 2,
#'   n_start = 100,
#'   n_step = 50,
#'   max_n_per_arm = 500
#' )
#' print(ss_results$results_data)
#' print(ss_results$results_plot)
#' }
linear.ss.boot <- function(pilot_data, time_var, status_var, arm_var,
                           target_power,
                           linear_terms = NULL, L, n_sim = 1000, alpha = 0.05,
                           patience = 5,
                           n_start = 50, n_step = 25, max_n_per_arm = 2000) {

   start_time <- proc.time()
   if (is.null(target_power) || length(target_power) != 1 || !is.numeric(target_power)) {
      stop("You must provide a single numeric value for 'target_power'.")
   }

   # --- Prepare data ---
   core_vars <- c(time_var, status_var, arm_var)
   all_vars <- c(core_vars, linear_terms)
   pilot_data <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]
   pilot_groups <- split(pilot_data, pilot_data[[arm_var]])
   model_rhs <- paste(c(arm_var, linear_terms), collapse = " + ")
   model_formula <- as.formula(paste("Y_rmst ~", model_rhs))

   # --- Run simulation ---
   cat("--- Searching for Sample Size (Method: Linear RMST with IPCW) ---\n")
   message("Model: Y_rmst ~ ", model_rhs)

   # Initialize variables for the search
   cat(paste0("\n--- Searching for N for ", target_power * 100, "% Power ---\n"))
   current_n <- n_start
   max_power_so_far <- -1
   stagnation_counter <- 0
   n_at_max_power <- n_start
   best_sim_output <- NULL
   search_results <- list()
   final_n <- NA

   # --- Iterative Search Loop ---
   while (current_n <= max_n_per_arm) {
      p_values <- rep(NA_real_, n_sim)
      estimates <- rep(NA_real_, n_sim)
      std_errors <- rep(NA_real_, n_sim)

      for (j in seq_len(n_sim)) {
         boot_list <- lapply(pilot_groups, function(df) df[sample(seq_len(nrow(df)), size = current_n, replace = TRUE), ])
         boot_data <- do.call(rbind, boot_list)
         boot_data[[arm_var]] <- factor(boot_data[[arm_var]], levels = c(0, 1))
         is_censored <- boot_data[[status_var]] == 0
         cens_fit <- tryCatch(survival::survfit(Surv(boot_data[[time_var]], is_censored) ~ 1), error = function(e) NULL)
         if (is.null(cens_fit)) next
         surv_summary <- tryCatch(summary(cens_fit, times = pmin(boot_data[[time_var]], L), extend = TRUE), error = function(e) NULL)
         if (is.null(surv_summary)) next
         weights <- 1 / surv_summary$surv
         finite_weights <- weights[is.finite(weights)]
         if (length(finite_weights) > 0) {
            weight_cap <- stats::quantile(finite_weights, probs = 0.99, na.rm = TRUE)
            weights[weights > weight_cap] <- weight_cap
         }
         weights[!is.finite(weights)] <- NA
         boot_data$Y_rmst <- pmin(boot_data[[time_var]], L)
         fit_data <- boot_data[boot_data[[status_var]] == 1 & is.finite(weights), ]
         fit_weights <- weights[boot_data[[status_var]] == 1 & is.finite(weights)]

         if (nrow(fit_data) > (length(all_vars) + 1)) {
            fit <- tryCatch(lm(model_formula, data = fit_data, weights = fit_weights), error = function(e) NULL)
            if (!is.null(fit)) {
               sfit <- summary(fit)
               test_term_pattern <- paste0("^", arm_var, "1$")
               test_term <- grep(test_term_pattern, rownames(sfit$coefficients), value = TRUE)[1]
               if (!is.na(test_term)) {
                  p_values[j] <- tryCatch(sfit$coefficients[test_term, "Pr(>|t|)"], error = function(e) NA)
                  estimates[j] <- tryCatch(sfit$coefficients[test_term, "Estimate"], error = function(e) NA)
                  std_errors[j] <- tryCatch(sfit$coefficients[test_term, "Std. Error"], error = function(e) NA)
               }
            }
         }
      }

      sim_output <- list(power = mean(p_values < alpha, na.rm = TRUE),
                         estimates = estimates, std_errors = std_errors)

      calculated_power <- if(is.finite(sim_output$power)) sim_output$power else 0
      search_results[[as.character(current_n)]] <- calculated_power
      cat(paste0("  N = ", current_n, "/arm, Calculated Power = ", round(calculated_power, 3), "\n"))

      if (calculated_power >= target_power) {
         message("Success: Target power reached at N = ", current_n, "/arm.")
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
         final_n <- n_at_max_power
         warning(paste0("Search terminated due to stagnation. Returning best N found: ", final_n,
                        " (Power = ", round(max_power_so_far, 3), ")"), call. = FALSE)
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
      est <- na.omit(best_sim_output$estimates)
      se  <- na.omit(best_sim_output$std_errors)
      if(length(est) > 1) {
         results_summary <- data.frame(
            Statistic = c("Mean RMST Difference", "Mean Standard Error", "95% CI Lower", "95% CI Upper"),
            Value = c(mean(est), mean(se, na.rm=TRUE), mean(est) - 1.96 * sd(est), mean(est) + 1.96 * sd(est))
         )
      }
   }

   results_df <- data.frame(Target_Power = target_power, Required_N_per_Arm = final_n)
   search_path_df <- data.frame(N_per_Arm = as.integer(names(search_results)),
                                Power = unlist(search_results))

   p <- ggplot2::ggplot(na.omit(search_path_df), ggplot2::aes(x = N_per_Arm, y = Power)) +
      ggplot2::geom_line(color = "#009E73", linewidth = 1) +
      ggplot2::geom_point(color = "#009E73", size = 3) +
      ggplot2::geom_hline(yintercept = target_power, linetype = "dashed", color = "red") +
      ggplot2::geom_vline(xintercept = final_n, linetype = "dotted", color = "blue") +
      ggplot2::labs(title = "Sample Size Search Path: Linear IPCW RMST Model",
                    x = "Sample Size Per Arm", y = "Calculated Power") +
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

   return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}


#' @title Analyze Power for a Simple Linear RMST Model
#' @description Performs a power analysis for given sample sizes based on the direct
#'   linear regression model for RMST (Tian et al., 2014).
#'
#' @note `status_var` should be `1` for an event, `0` for censored. `arm_var`
#'   should be `1` for treatment, `0` for control.
#'
#' @param pilot_data A data.frame with pilot study data.
#' @param time_var,status_var,arm_var Strings for column names.
#' @param sample_sizes A numeric vector of sample sizes per arm to calculate power for.
#' @param linear_terms Optional vector of covariates for the linear model. If NULL,
#'   defaults to all other columns in `pilot_data`.
#' @param tau The truncation time for RMST.
#' @param n_sim Number of bootstrap simulations.
#' @param alpha The significance level.
#'
#' @return A list containing `results_data` (a data.frame of sample sizes and corresponding power),
#'   `results_plot` (a ggplot object), and `results_summary` (a data.frame with the estimated
#'   treatment effect statistics from the largest sample size simulation).
#'
#' @references Tian, L., et al. (2014). *Biostatistics*, 15(2), 222-233.
#'
#' @importFrom survival survfit Surv
#' @importFrom stats lm as.formula complete.cases na.omit sd quantile proc.time
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline labs theme_minimal ylim
#' @importFrom knitr kable
#'
#' @examples
#' # Create a small pilot dataset for the example
#' pilot_data_linear <- data.frame(
#'   event_time = rexp(100, rate = 0.1),
#'   event_status = rbinom(100, 1, 0.7),
#'   treatment_arm = rep(0:1, each = 50),
#'   age = rnorm(100, 50, 10)
#' )
#'
#' # --- Power Calculation Example ---
#' results_power <- design_rmst_ipcw_power(
#'   pilot_data = pilot_data_linear,
#'   time_var = "event_time",
#'   status_var = "event_status",
#'   arm_var = "treatment_arm",
#'   sample_sizes = c(100, 150),
#'   tau = 10,
#'   n_sim = 10 # Low n_sim for example speed
#' )
#'
#' @export
design_rmst_ipcw_power <- function(pilot_data, time_var, status_var, arm_var,
                                   sample_sizes,linear_terms = NULL, tau, n_sim = 1000, alpha = 0.05)
{

   start_time <- proc.time()
   if (is.null(sample_sizes)) stop("You must provide a numeric vector for 'sample_sizes'.")

   # --- Prepare data ---
   core_vars <- c(time_var, status_var, arm_var)
   if (is.null(linear_terms)) {
      linear_terms <- setdiff(names(pilot_data), core_vars)
      if (length(linear_terms) > 0) message("Using unspecified columns as linear terms: ", paste(linear_terms, collapse = ", "))
   }

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
         surv_summary <- tryCatch(summary(cens_fit, times = pmin(boot_data[[time_var]], tau), extend = TRUE), error = function(e) NULL)
         if (is.null(surv_summary)) next

         weights <- 1 / surv_summary$surv
         finite_weights <- weights[is.finite(weights)]
         if (length(finite_weights) > 0) {
            weight_cap <- stats::quantile(finite_weights, probs = 0.99, na.rm = TRUE)
            weights[weights > weight_cap] <- weight_cap
         }
         weights[!is.finite(weights)] <- NA

         boot_data$Y_rmst <- pmin(boot_data[[time_var]], tau)
         fit_data <- boot_data[boot_data[[status_var]] == 1 & is.finite(weights), ]
         fit_weights <- weights[boot_data[[status_var]] == 1 & is.finite(weights)]

         if (nrow(fit_data) > (length(all.vars(model_formula)) - 1)) {
            fit <- tryCatch(lm(model_formula, data = fit_data, weights = fit_weights), error = function(e) NULL)
            if (!is.null(fit)) {
               sfit <- summary(fit)
               test_term <- grep(paste0("^", arm_var), rownames(sfit$coefficients), value = TRUE)[1]
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
      # ggplot2::geom_hline(yintercept = 0.75, linetype = "dashed", color = "blue") +
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
   print(p)

   return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}
#' @title Find Sample Size for a Linear RMST Model
#' @description Performs an iterative sample size search to achieve a target power
#'   based on the direct linear regression model for RMST.
#'
#' @note `status_var` should be `1` for an event, `0` for censored. `arm_var`
#'   should be `1` for treatment, `0` for control.
#'
#' @param pilot_data A data.frame with pilot study data.
#' @param time_var,status_var,arm_var Strings for column names.
#' @param target_power A single numeric value for the target power.
#' @param linear_terms Optional vector of covariates for the linear model.
#' @param tau The truncation time for RMST.
#' @param n_sim Number of bootstrap simulations per search step.
#' @param alpha The significance level.
#' @param patience Number of consecutive non-improving steps in the iterative
#'   search before terminating. Default is 5.
#' @param n_start The starting sample size per arm for the search.
#' @param n_step The increment in sample size at each step of the search.
#' @param max_n_per_arm The maximum sample size per arm to search up to.
#'
#' @return A list containing `results_data` (the final required N), a ggplot object
#'   of the search path, and a `results_summary` data.frame.
#'
#' @importFrom survival survfit Surv
#' @importFrom stats lm as.formula complete.cases na.omit sd quantile proc.time
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline geom_vline labs theme_minimal
#' @importFrom knitr kable
#'
#' @examples
#' # Use the 'veteran' dataset for a stable example
#' data(veteran)
#' veteran_pilot <- veteran %>%
#'   dplyr::mutate(
#'     status = ifelse(status == 1, 1, 0),
#'     trt = ifelse(trt == 2, 1, 0)
#'   ) %>%
#'   tidyr::drop_na(time, status, trt, karno, age)
#'
#' results_n <- design_rmst_ipcw_ss(
#'   pilot_data = veteran_pilot,
#'   time_var = "time",
#'   status_var = "status",
#'   arm_var = "trt",
#'   target_power = 0.80,
#'   linear_terms = c("age", "karno"),
#'   tau = 365,
#'   n_sim = 50, # Low n_sim for example speed
#'   patience = 3,
#'   n_step = 100
#' )
#' @export
design_rmst_ipcw_ss <- function(pilot_data, time_var, status_var, arm_var,
                                target_power,
                                linear_terms = NULL, tau, n_sim = 1000, alpha = 0.05,
                                patience = 5,
                                n_start = 50, n_step = 25, max_n_per_arm = 2000) {

   start_time <- proc.time()
   if (is.null(target_power) || length(target_power) != 1 || !is.numeric(target_power)) {
      stop("You must provide a single numeric value for 'target_power'.")
   }

   # --- Prepare data ---
   core_vars <- c(time_var, status_var, arm_var)
   if (is.null(linear_terms)) {
      linear_terms <- setdiff(names(pilot_data), core_vars)
      if (length(linear_terms) > 0) message("Using unspecified columns as linear terms: ", paste(linear_terms, collapse = ", "))
   }
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
         surv_summary <- tryCatch(summary(cens_fit, times = pmin(boot_data[[time_var]], tau), extend = TRUE), error = function(e) NULL)
         if (is.null(surv_summary)) next
         weights <- 1 / surv_summary$surv
         finite_weights <- weights[is.finite(weights)]
         if (length(finite_weights) > 0) {
            weight_cap <- stats::quantile(finite_weights, probs = 0.99, na.rm = TRUE)
            weights[weights > weight_cap] <- weight_cap
         }
         weights[!is.finite(weights)] <- NA
         boot_data$Y_rmst <- pmin(boot_data[[time_var]], tau)
         fit_data <- boot_data[boot_data[[status_var]] == 1 & is.finite(weights), ]
         fit_weights <- weights[boot_data[[status_var]] == 1 & is.finite(weights)]
         if (nrow(fit_data) > (length(all.vars(model_formula)) - 1)) {
            fit <- tryCatch(lm(model_formula, data = fit_data, weights = fit_weights), error = function(e) NULL)
            if (!is.null(fit)) {
               sfit <- summary(fit)
               test_term <- grep(paste0("^", arm_var), rownames(sfit$coefficients), value = TRUE)[1]
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

   print(p)

   return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}

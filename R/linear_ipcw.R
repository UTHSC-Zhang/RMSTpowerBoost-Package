#' @title Analyze Power/Sample Size using a Simple Linear RMST Model
#' @description Performs power or sample size analysis based on the direct
#'   linear regression model for RMST (Tian et al., 2014). This method is
#'   best suited for studies where censoring can be assumed to be independent
#'   of covariates.
#'
#' @note `status_var` should be `1` for an event, `0` for censored. `arm_var`
#'   should be `1` for treatment, `0` for control.
#'
#' @param pilot_data A data.frame with pilot study data.
#' @param time_var,status_var,arm_var Strings for column names.
#' @param sample_sizes Optional vector of sample sizes per arm.
#' @param target_powers Optional vector of target powers.
#' @param linear_terms Optional vector of covariates for the linear model. If NULL,
#'   defaults to all other columns in `pilot_data`.
#' @param tau The truncation time for RMST.
#' @param n_sim Number of bootstrap simulations.
#' @param alpha The significance level.
#' @param parallel.cores Number of cores for parallel processing. Default is 1 (no parallel).
#' @param patience Number of consecutive non-improving steps in the sample size
#'   search before terminating that search. Default is 5.
#' @param ... Additional arguments for sample size search (n_start, n_step, etc.).
#'
#' @return A list containing `results_data`, `results_plot`, and `results_summary`.
#'
#' @references Tian, L., et al. (2014). *Biostatistics*, 15(2), 222-233.
#'
#' @importFrom survival survfit Surv
#' @importFrom stats lm as.formula complete.cases na.omit
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline labs theme_minimal ylim
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
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
#' # Fast-running example for CRAN checks
#' results_power <- design_rmst_linear_ipcw(
#'   pilot_data = pilot_data_linear,
#'   time_var = "event_time",
#'   status_var = "event_status",
#'   arm_var = "treatment_arm",
#'   sample_sizes = c(50, 100),
#'   tau = 10,
#'   n_sim = 10 # Low n_sim for example speed
#' )
#'
#' \donttest{
#' # --- Sample Size Calculation with Parallel Processing ---
#' # This example is more intensive and demonstrates parallel computing.
#' # It will not be run on CRAN.
#' results_n <- design_rmst_linear_ipcw(
#'   pilot_data = pilot_data_linear,
#'   time_var = "event_time",
#'   status_var = "event_status",
#'   arm_var = "treatment_arm",
#'   target_powers = c(0.80),
#'   tau = 10,
#'   n_sim = 100,
#'   parallel.cores = 2 # Use 2 cores
#' )
#' }
#'
#' @export
design_rmst_linear_ipcw <- function(pilot_data, time_var, status_var, arm_var,
                                    sample_sizes = NULL, target_powers = NULL,
                                    linear_terms = NULL, tau, n_sim = 1000, alpha = 0.05,
                                    parallel.cores = 1, patience = 5,
                                    n_start = 50, n_step = 25, max_n_per_arm = 2000) {
   start_time <- Sys.time()
   if (is.null(sample_sizes) && is.null(target_powers)) stop("Must provide 'sample_sizes' or 'target_powers'.")

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

   # --- Simulation function ---
   run_power_sim <- function(n_per_arm) {
      for (i in seq_len(n_sim)) {
         # Bootstrap
         boot_data <- do.call(rbind, lapply(pilot_groups, function(df) {
            df[sample(seq_len(nrow(df)), size = n_per_arm, replace = TRUE), ]
         }))

         # Make sure arm_var is a factor with levels 0 and 1
         boot_data[[arm_var]] <- factor(boot_data[[arm_var]], levels = c(0, 1))

         # IPCW weights
         is_censored <- boot_data[[status_var]] == 0
         surv_obj <- survival::Surv(boot_data[[time_var]], is_censored)
         cens_fit <- tryCatch(survival::survfit(surv_obj ~ 1), error = function(e) NULL)
         if (is.null(cens_fit)) next

         surv_summary <- tryCatch(summary(cens_fit, times = pmin(boot_data[[time_var]], tau), extend = TRUE), error = function(e) NULL)
         if (is.null(surv_summary)) next

         weights <- 1 / surv_summary$surv
         weights[!is.finite(weights) | weights > 1e4] <- NA
         boot_data$Y_rmst <- pmin(boot_data[[time_var]], tau)

         # Use only uncensored observations
         fit_data <- boot_data[boot_data[[status_var]] == 1, ]
         fit_weights <- weights[boot_data[[status_var]] == 1]

         valid_idx <- is.finite(fit_weights)
         fit_data <- fit_data[valid_idx, ]
         fit_weights <- fit_weights[valid_idx]

         # if (nrow(fit_data) < length(linear_terms) + 2) next

         fit <- tryCatch(lm(model_formula, data = fit_data, weights = fit_weights), error = function(e) NULL)
         if (is.null(fit)) next

         sfit <- summary(fit)
         term_names <- rownames(sfit$coefficients)

         # Find test term robustly
         test_term <- grep(paste0("^", arm_var), term_names, value = TRUE)[1]
         if (is.na(test_term)) {
            cat("Skipped sim", i, ": test_term not found in model coefficients\n")
            next
         }

         # Debug output
         cat("Sim", i, ": Coefs =", paste(term_names, collapse = ", "), "\n")
         cat("Using test term:", test_term, "\n")

         p_val <- tryCatch(sfit$coefficients[test_term, "Pr(>|t|)"], error = function(e) NA)
         estimate <- tryCatch(sfit$coefficients[test_term, "Estimate"], error = function(e) NA)
         std_error <- tryCatch(sfit$coefficients[test_term, "Std. Error"], error = function(e) NA)

         p_values[i] <- p_val
         estimates[i] <- estimate
         std_errors[i] <- std_error
      }

      return(list(
         power = mean(p_values < alpha, na.rm = TRUE),
         estimates = estimates,
         std_errors = std_errors
      ))
   }

   # --- Run simulation ---
   cat("--- Calculating Power/Sample Size (Method: Linear RMST with IPCW) ---\n")
   message("Model: Y_rmst ~ ", model_rhs)
   results_summary <- NULL

   if (!is.null(sample_sizes)) {
      results_list <- vector("list", length(sample_sizes))

      for (i in seq_along(sample_sizes)) {
         n <- sample_sizes[i]
         cat("Simulating for n =", n, "per arm...\n")
         results_list[[i]] <- run_power_sim(n)
      }

      power_values <- sapply(results_list, function(x) x$power)
      results_df <- data.frame(N_per_Arm = sample_sizes, Power = power_values)

      # Summary from the largest n
      best <- which.max(sample_sizes)
      est <- na.omit(results_list[[best]]$estimates)
      se  <- na.omit(results_list[[best]]$std_errors)

      if (length(est) > 1) {
         results_summary <- data.frame(
            Statistic = c("Mean RMST Difference", "Mean Standard Error", "95% CI Lower", "95% CI Upper"),
            Value = c(mean(est), mean(se), mean(est) - 1.96 * sd(est), mean(est) + 1.96 * sd(est))
         )
      }

      p <- ggplot2::ggplot(results_df, ggplot2::aes(x = N_per_Arm, y = Power)) +
         ggplot2::geom_line(color = "#D55E00", linewidth = 1) +
         ggplot2::geom_point(color = "#D55E00", size = 3) +
         ggplot2::labs(title = "Power Curve: Linear IPCW RMST Model",
                       x = "Sample Size Per Arm", y = "Estimated Power") +
         ggplot2::ylim(0, 1) + ggplot2::theme_minimal()

   }
   else {
      required_n_values <- sapply(sort(target_powers), function(target_power) {
         cat(paste0("\n--- Searching for N for ", target_power * 100, "% Power ---\n"))
         current_n <- n_start
         max_power_so_far <- -1
         stagnation_counter <- 0
         n_at_max_power <- n_start
         best_sim_output <- NULL # <<-- NEW: Store the best simulation output

         while (current_n <= max_n_per_arm) {
            sim_output <- run_power_sim(current_n)
            calculated_power <- sim_output$power
            if (!is.finite(calculated_power)) {
               calculated_power <- 0
            }
            cat(paste0("  N = ", current_n, "/arm, Calculated Power = ", round(calculated_power, 3), "\n"))

            if (calculated_power >= target_power) {
               est <- na.omit(sim_output$estimates)
               se <- na.omit(sim_output$std_errors)
               if(length(est) > 1) {
                  results_summary <<- data.frame(
                     Statistic = c("Mean RMST Difference", "Mean Standard Error", "95% CI Lower", "95% CI Upper"),
                     Value = c(mean(est), mean(se, na.rm=TRUE), mean(est) - 1.96 * sd(est), mean(est) + 1.96 * sd(est))
                  )
               }
               return(current_n)
            }

            if (calculated_power > max_power_so_far) {
               max_power_so_far <- calculated_power
               n_at_max_power <- current_n
               best_sim_output <- sim_output # <<-- NEW: Save the best result
               stagnation_counter <- 0
            } else {
               stagnation_counter <- stagnation_counter + 1
            }

            if (stagnation_counter >= patience) {
               warning(paste0("Search for target power ", target_power, " terminated due to stagnation. ",
                              "Returning N=", n_at_max_power, " which achieved the highest power of ", round(max_power_so_far, 3), "."),
                       call. = FALSE)

               # --- NEW: Generate summary from the best result found ---
               if (!is.null(best_sim_output)) {
                  est <- na.omit(best_sim_output$estimates)
                  se <- na.omit(best_sim_output$std_errors)
                  if(length(est) > 1) {
                     results_summary <<- data.frame(
                        Statistic = c("Mean RMST Difference", "Mean Standard Error", "95% CI Lower", "95% CI Upper"),
                        Value = c(mean(est), mean(se, na.rm=TRUE), mean(est) - 1.96 * sd(est), mean(est) + 1.96 * sd(est))
                     )
                  }
               }
               # --- End NEW --
               return(n_at_max_power)
            }
            current_n <- current_n + n_step
         }
         warning(paste("Target power", target_power, "not achieved by max N of", max_n_per_arm))
         return(NA)
      })

      p <- ggplot2::ggplot(na.omit(results_df), ggplot2::aes(x = Target_Power, y = Required_N_per_Arm)) +
         ggplot2::geom_line(color = "#009E73", linewidth = 1) +
         ggplot2::geom_point(color = "#009E73", size = 3) +
         ggplot2::labs(title = "Sample Size Curve: Linear IPCW RMST Model",
                       x = "Target Power", y = "Required Sample Size Per Arm") +
         ggplot2::theme_minimal()
   }

   cat("\n--- Simulation Summary ---\n")
   if (!is.null(results_summary)) {
      print(knitr::kable(results_summary, caption = "Estimated Treatment Effect (RMST Difference)"))
   }
   else {
      cat("No valid estimates were generated to create a summary.\n")
   }
   print(p)
   end_time <- Sys.time()
   return(list(results_data = results_df, results_plot = p, results_summary = results_summary,
               execution_time = end_time - start_time))
}

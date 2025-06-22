#' @title Analyze Power/Sample Size for RMST with Dependent Censoring
#' @description Performs power analysis for a two-arm study with a competing risk,
#'   implementing the double IPCW method of Wang & Schaubel (2018).
#'
#' @note The primary event in `status_var` should be `1`/`0`. The dependent censoring
#'   event in `dep_cens_status_var` should be `1`/`0`. `arm_var` should be `1`/`0`.
#'
#' @param pilot_data A data.frame with pilot study data.
#' @param time_var,status_var,arm_var,dep_cens_status_var Strings for column names.
#' @param sample_sizes Optional vector of sample sizes per arm.
#' @param target_powers Optional vector of target powers.
#' @param linear_terms Optional vector of covariates for the model. If NULL,
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
#' @references Wang, X., & Schaubel, D. E. (2018). *Lifetime Data Analysis*, 24, 176-199.
#'
#' @importFrom survival coxph Surv basehaz
#' @importFrom stats lm as.formula complete.cases predict na.omit sd
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline labs theme_minimal ylim
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#'
#' @examples
#' # Create a small pilot dataset with a competing risk
#' set.seed(123) # for reproducibility
#' pilot_data_cr <- data.frame(
#'   time = rexp(200, rate = 0.1),
#'   arm = rep(0:1, each = 100),
#'   age = rnorm(200, 60, 8)
#' )
#' event_probs <- runif(200)
#' pilot_data_cr$status_primary <- ifelse(event_probs < 0.35, 1, 0)
#' pilot_data_cr$status_competing <- ifelse(event_probs >= 0.35 & event_probs < 0.60, 1, 0)
#' pilot_data_cr$status_primary[pilot_data_cr$status_competing == 1] <- 0
#'
#' # --- Power Calculation Example ---
#' results_power_cr <- design_rmst_dependent_censoring(
#'   pilot_data = pilot_data_cr,
#'   time_var = "time",
#'   status_var = "status_primary",
#'   dep_cens_status_var = "status_competing",
#'   arm_var = "arm",
#'   sample_sizes = c(250),
#'   tau = 10,
#'   n_sim = 20 # Low n_sim for example speed
#' )
#'
#' \donttest{
#' # --- Sample Size Calculation with Parallel Processing ---
#' results_n_cr <- design_rmst_dependent_censoring(
#'   pilot_data = pilot_data_cr,
#'   time_var = "time",
#'   status_var = "status_primary",
#'   dep_cens_status_var = "status_competing",
#'   arm_var = "arm",
#'   target_powers = c(0.80),
#'   tau = 10,
#'   n_sim = 100,
#'   parallel.cores = 2
#' )
#' }
#'
#' @export
design_rmst_dependent_censoring <- function(pilot_data, time_var, status_var, arm_var, dep_cens_status_var,
                                            sample_sizes = NULL, target_powers = NULL,
                                            linear_terms = NULL, tau, n_sim = 1000, alpha = 0.05,
                                            parallel.cores = 1, patience = 5,
                                            n_start = 50, n_step = 25, max_n_per_arm = 2000) {

   # --- Input Validation & Package Checks ---
   if (parallel.cores > 1) {
      if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
         stop("Packages 'future' and 'future.apply' are required for parallel processing.")
      }
   }
   if (is.null(sample_sizes) && is.null(target_powers)) stop("Must provide 'sample_sizes' or 'target_powers'.")

   # --- Data Preparation & Model Formula ---
   core_vars <- c(time_var, status_var, arm_var, dep_cens_status_var)
   if (is.null(linear_terms)) {
      linear_terms <- setdiff(names(pilot_data), core_vars)
      if(length(linear_terms) > 0) message(paste("Using unspecified columns as linear terms:", paste(linear_terms, collapse=", ")))
   }
   all_vars <- c(core_vars, linear_terms)
   cleaned_pilot_data <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]
   pilot_groups <- split(cleaned_pilot_data, cleaned_pilot_data[[arm_var]])

   cox_rhs <- paste(c(arm_var, linear_terms), collapse=" + ")
   cens_formula <- as.formula(paste("Surv(", time_var, ", cens_ind) ~", cox_rhs))
   dep_cens_formula <- as.formula(paste("Surv(", time_var, ", dep_cens_ind) ~", cox_rhs))
   model_formula <- as.formula(paste("Y_rmst ~", cox_rhs))
   test_term <- arm_var

   # --- Internal Simulation Function ---
   run_power_sim <- function(n_per_arm) {

      single_iteration <- function() {
         boot_list <- lapply(pilot_groups, function(group_df) group_df[sample(1:nrow(group_df), size = n_per_arm, replace = TRUE), ])
         boot_data <- do.call(rbind, boot_list)

         boot_data$cens_ind <- boot_data[[status_var]] == 0 & boot_data[[dep_cens_status_var]] == 0
         fit_cens <- tryCatch(survival::coxph(cens_formula, data = boot_data, ties="breslow"), error=function(e) NULL)

         boot_data$dep_cens_ind <- boot_data[[dep_cens_status_var]] == 1
         fit_dep_cens <- tryCatch(survival::coxph(dep_cens_formula, data = boot_data, ties="breslow"), error=function(e) NULL)

         if(is.null(fit_cens) || is.null(fit_dep_cens)) {
            return(list(p_value = NA, estimate = NA, std_error = NA))
         }

         bh_cens <- tryCatch(survival::basehaz(fit_cens, centered = FALSE), error = function(e) NULL)
         if(is.null(bh_cens)) { H_cens <- rep(0, nrow(boot_data)); } else {
            cumhaz_func_cens <- stepfun(bh_cens$time, c(0, bh_cens$hazard))
            H_cens <- cumhaz_func_cens(pmin(boot_data[[time_var]], tau)) * exp(predict(fit_cens, type="lp"))
         }

         bh_dep_cens <- tryCatch(survival::basehaz(fit_dep_cens, centered = FALSE), error=function(e) NULL)
         if(is.null(bh_dep_cens)) { H_dep_cens <- rep(0, nrow(boot_data)); } else {
            cumhaz_func_dep <- stepfun(bh_dep_cens$time, c(0, bh_dep_cens$hazard))
            H_dep_cens <- cumhaz_func_dep(pmin(boot_data[[time_var]], tau)) * exp(predict(fit_dep_cens, type="lp"))
         }

         weights <- exp(H_cens + H_dep_cens)
         weights[!is.finite(weights)] <- max(weights[is.finite(weights)], 0, na.rm = TRUE)
         boot_data$Y_rmst <- pmin(boot_data[[time_var]], tau)

         # The regression must ONLY be fit on subjects who had the primary event.
         fit_data <- boot_data[boot_data[[status_var]] == 1, ]
         fit_weights <- weights[boot_data[[status_var]] == 1]

         if (nrow(fit_data) > (length(linear_terms) + 2)) {
            fit <- tryCatch(lm(model_formula, data=fit_data, weights=fit_weights), error = function(e) NULL)
            if (!is.null(fit)) {
               sfit <- summary(fit)
               p_val <- tryCatch(sfit$coefficients[test_term, "Pr(>|t|)"], error = function(e) NA)
               estimate <- tryCatch(sfit$coefficients[test_term, "Estimate"], error = function(e) NA)
               std_error <- tryCatch(sfit$coefficients[test_term, "Std. Error"], error = function(e) NA)

            }
         }
         return(list(p_value = p_val, estimate = estimate, std_error = std_error))
      }

      if (parallel.cores > 1) {
         old_plan <- future::plan(future::multisession, workers = parallel.cores)
         on.exit(future::plan(old_plan), add = TRUE)
         all_sim_outputs <- future.apply::future_lapply(1:n_sim, function(i) single_iteration(), future.seed = TRUE)
      } else {
         all_sim_outputs <- lapply(1:n_sim, function(i) single_iteration())
      }

      p_values <- vapply(all_sim_outputs, `[[`, "p_value", FUN.VALUE = numeric(1))
      estimates <- vapply(all_sim_outputs, `[[`, "estimate", FUN.VALUE = numeric(1))
      std_errors <- vapply(all_sim_outputs, `[[`, "std_error", FUN.VALUE = numeric(1))
      power <- mean(p_values < alpha, na.rm = TRUE)
      return(list(power = power, estimates = estimates, std_errors = std_errors))
   }

   # --- Main Execution Logic ---
   cat("--- Calculating Power/Sample Size (Method: Dependent Censoring RMST Model) ---\n")
   message("Model: Y_rmst ~ ", cox_rhs)

   results_summary <- NULL

   if (!is.null(sample_sizes)) {
      # Power Calculation Mode
      sim_outputs <- lapply(sample_sizes, function(n) {
         cat(paste0("Simulating for n = ", n, "/arm...\n"))
         run_power_sim(n)
      })
      power_values <- sapply(sim_outputs, `[[`, "power")
      results_df <- data.frame(N_per_Arm = sample_sizes, Power = power_values)

      est <- na.omit(sim_outputs[[1]]$estimates)
      se <- na.omit(sim_outputs[[1]]$std_errors)
      if(length(est) > 1){
         results_summary <- data.frame(
            Statistic = c("Mean RMST Difference", "Mean Standard Error", "95% CI Lower", "95% CI Upper"),
            Value = c(mean(est), mean(se, na.rm=TRUE), mean(est) - 1.96 * sd(est), mean(est) + 1.96 * sd(est))
         )
      }

      p <- ggplot2::ggplot(results_df, ggplot2::aes(x = N_per_Arm, y = Power)) +
         ggplot2::geom_line(color = "#56B4E9", linewidth = 1) +
         ggplot2::geom_point(color = "#56B4E9", size = 3) +
         ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
         ggplot2::labs(title = "Power Curve: Dependent Censoring RMST Model",
                       x = "Sample Size Per Arm", y = "Estimated Power") +
         ggplot2::ylim(0, 1) + ggplot2::theme_minimal()

   } else {
      # Sample Size Calculation Mode
      required_n_values <- sapply(sort(target_powers), function(target_power) {
         cat(paste0("\n--- Searching for N for ", target_power * 100, "% Power ---\n"))
         current_n <- n_start
         max_power_so_far <- -1
         stagnation_counter <- 0
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
               stagnation_counter <- 0
            } else {
               stagnation_counter <- stagnation_counter + 1
            }

            if (stagnation_counter >= patience) {
               warning(paste("Search for target power", target_power, "terminated: power has not improved for", patience, "steps."), call. = FALSE)
               break # Exit the while loop for this target_power
            }


            current_n <- current_n + n_step
         }
         warning(paste("Target power", target_power, "not achieved by max N of", max_n_per_arm))
         return(NA)
      })
      results_df <- data.frame(Target_Power = sort(target_powers), Required_N_per_Arm = required_n_values)

      p <- ggplot2::ggplot(na.omit(results_df), ggplot2::aes(x = Target_Power, y = Required_N_per_Arm)) +
         ggplot2::geom_line(color = "#009E73", linewidth = 1) +
         ggplot2::geom_point(color = "#009E73", size = 3) +
         ggplot2::geom_vline(xintercept = 0.8, linetype = "dashed", color = "red") +
         ggplot2::labs(title = "Sample Size Curve: Dependent Censoring RMST Model",
                       x = "Target Power", y = "Required Sample Size Per Arm") +
         ggplot2::theme_minimal()
   }

   cat("\n--- Simulation Summary ---\n")
   if(!is.null(results_summary)) {
      print(knitr::kable(results_summary, caption = "Estimated Treatment Effect (RMST Difference)"))
   } else {
      cat("No valid estimates were generated to create a summary.\n")
   }
   cat("\n")
   print(p)

   return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}

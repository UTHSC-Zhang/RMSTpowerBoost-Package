#' @title Analyze Power/Sample Size for an Additive RMST Model
#' @description Performs power or sample size analysis for a two-arm study using
#'   a flexible, semiparametric additive model for the RMST. This method allows
#'   for non-linear covariate effects via splines and supports optional stratification.
#'
#' @details This function uses a bootstrap-based approach on a pilot dataset.
#'   The core estimation is performed by fitting a Generalized Additive Model
#'   (GAM) from the `mgcv` package to jackknife pseudo-observations of the RMST.
#'   The user can specify which covariates have linear effects and which have
#'   smooth, non-linear effects.
#'
#' @note This function assumes the data is pre-processed.
#'   `status_var` should be coded as `1` for an event and `0` for censored.
#'   `arm_var` should be coded as `1` for the treatment/experimental group and `0` for the control group.
#'
#' @param pilot_data A data.frame containing pilot study data.
#' @param time_var,status_var,arm_var Strings specifying column names.
#' @param strata_var An optional string for a stratification variable.
#' @param sample_sizes Optional numeric vector of sample sizes per arm/stratum.
#' @param target_powers Optional numeric vector of target power levels.
#' @param linear_terms Optional character vector of covariates with a linear effect.
#' @param smooth_terms Optional character vector of covariates with a non-linear (spline) effect.
#' @param tau The truncation time for RMST.
#' @param n_sim Number of bootstrap simulations.
#' @param alpha The significance level.
#' @param parallel.cores Number of cores for parallel processing. Default is 1 (no parallel).
#' @param patience Number of consecutive non-improving steps in the sample size
#'   search before terminating that search. Default is 5.
#' @param n_start,n_step,max_n_per_arm Search parameters for sample size mode.
#'
#' @return A list containing `results_data`, `results_plot`, and `results_summary`.
#'   `results_summary` provides the mean, standard error, and confidence interval for the
#'   estimated RMST difference between treatment arms.
#'
#' @references Zhang, Y., & Schaubel, D. E. (2024). *Biometrical Journal*.
#'
#' @importFrom survival survfit Surv
#' @importFrom mgcv gam summary.gam
#' @importFrom stats lm as.formula complete.cases integrate na.omit stepfun predict sd
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_vline geom_hline labs theme_minimal ylim
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#'
#' @examples
#' # Create a small pilot dataset for the example
#' pilot_data_gam <- data.frame(
#'   surv_time = rexp(100, rate = 0.1),
#'   censor_status = rbinom(100, 1, 0.7),
#'   arm = rep(0:1, each = 50),
#'   age = rnorm(100, 50, 10),
#'   biomarker = rnorm(100, 10, 3)
#' )
#'
#' # --- Power Calculation with Smooth Term ---
#' # Fast-running example for CRAN checks
#' results_power_gam <- design_rmst_additive_gam(
#'   pilot_data = pilot_data_gam,
#'   time_var = "surv_time",
#'   status_var = "censor_status",
#'   arm_var = "arm",
#'   linear_terms = "age",
#'   smooth_terms = "biomarker",
#'   sample_sizes = c(100),
#'   tau = 10,
#'   n_sim = 10 # Low n_sim for example speed
#' )
#'
#' \donttest{
#' # --- Sample Size Calculation with Parallel Processing ---
#' # This example is more intensive and demonstrates parallel computing.
#' results_n_gam <- design_rmst_additive_gam(
#'   pilot_data = pilot_data_gam,
#'   time_var = "surv_time",
#'   status_var = "censor_status",
#'   arm_var = "arm",
#'   target_powers = c(0.80),
#'   tau = 10,
#'   n_sim = 100,
#'   parallel.cores = 2
#' )
#' }
#'
#' @export
design_rmst_additive_gam <- function(pilot_data, time_var, status_var, arm_var, strata_var = NULL,
                                     sample_sizes = NULL, target_powers = NULL,
                                     linear_terms = NULL, smooth_terms = NULL,
                                     tau, n_sim = 1000, alpha = 0.05,
                                     parallel.cores = 1, patience = 5,
                                     n_start = 50, n_step = 25, max_n_per_arm = 2000) {

   # --- Input Validation & Package Checks ---
   if (is.null(sample_sizes) && is.null(target_powers)) stop("Must provide 'sample_sizes' or 'target_powers'.")
   if (!is.null(sample_sizes) && !is.null(target_powers)) stop("Cannot provide both 'sample_sizes' and 'target_powers'.")
   if (parallel.cores > 1) {
      if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
         stop("Packages 'future' and 'future.apply' are required for parallel processing.")
      }
   }
   if (!requireNamespace("survival", quietly = TRUE)) stop("'survival' package is required.")
   if (!requireNamespace("ggplot2", quietly = TRUE)) stop("'ggplot2' package is required.")
   if (!requireNamespace("mgcv", quietly = TRUE)) stop("'mgcv' package is required.")

   is_stratified <- !is.null(strata_var)

   # --- Data Preparation ---
   core_vars <- c(time_var, status_var, arm_var, strata_var)
   if (is.null(linear_terms) && is.null(smooth_terms)) {
      linear_terms <- setdiff(names(pilot_data), core_vars)
      if(length(linear_terms) > 0) message(paste("Using unspecified columns as linear terms:", paste(linear_terms, collapse=", ")))
   }
   all_vars <- c(core_vars, linear_terms, smooth_terms)
   missing_cols <- setdiff(all_vars, names(pilot_data))
   if (length(missing_cols) > 0) {
      stop(paste("The following required columns are not in pilot_data:", paste(missing_cols, collapse = ", ")))
   }
   if (!all(pilot_data[[status_var]] %in% c(0, 1))) stop("'status_var' must be coded as 0 (censored) and 1 (event).")
   if (!all(pilot_data[[arm_var]] %in% c(0, 1))) stop("'arm_var' must be coded as 0 (control) and 1 (treatment).")

   cleaned_pilot_data <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]

   # --- Helper to calculate jackknife pseudo-observations ---
   get_pseudo_obs <- function(time, status, tau) {
      n <- length(time)
      if (n == 0) return(numeric(0))
      km_fit_full <- survival::survfit(survival::Surv(time, status) ~ 1)
      km_step_full <- stepfun(km_fit_full$time, c(1, km_fit_full$surv))
      rmst_full <- tryCatch(integrate(km_step_full, 0, tau, subdivisions=2000, stop.on.error = FALSE)$value, error = function(e) 0)
      rmst_loo <- sapply(1:n, function(i) {
         if(n > 1) {
            km_fit_loo <- survival::survfit(survival::Surv(time[-i], status[-i]) ~ 1)
            km_step_loo <- stepfun(km_fit_loo$time, c(1, km_fit_loo$surv))
            tryCatch(integrate(km_step_loo, 0, tau, subdivisions=2000, stop.on.error = FALSE)$value, error = function(e) 0)
         } else { 0 }
      })
      return(n * rmst_full - (n - 1) * rmst_loo)
   }

   # --- Model Formula Construction ---
   smooth_part <- if (!is.null(smooth_terms)) paste0("s(", smooth_terms, ")") else NULL

   if (is_stratified) {
      pilot_groups <- split(cleaned_pilot_data, cleaned_pilot_data[[strata_var]])
      all_terms <- c(strata_var, paste0(arm_var, ":", strata_var), smooth_part, linear_terms)
      test_term_pattern <- paste0(arm_var, "1:", strata_var) # Pattern for interaction term
   } else {
      pilot_groups <- list(cleaned_pilot_data) # Treat as a single stratum
      all_terms <- c(arm_var, smooth_part, linear_terms)
      test_term_pattern <- arm_var
   }
   model_rhs <- paste(all_terms[!sapply(all_terms, is.null)], collapse = " + ")
   model_formula <- as.formula(paste("pseudo_obs ~", model_rhs))

   # --- Internal Simulation Function ---
   run_power_sim <- function(n_per_group) {

      single_iteration <- function() {
         # Stratified or non-stratified sampling
         if (is_stratified) {
            boot_list <- lapply(pilot_groups, function(group_df) {
               group_df[sample(1:nrow(group_df), size = n_per_group, replace = TRUE), ]
            })
         } else {
            group_arms <- split(pilot_groups[[1]], pilot_groups[[1]][[arm_var]])
            boot_list <- lapply(group_arms, function(arm_df) {
               arm_df[sample(1:nrow(arm_df), size = n_per_group, replace = TRUE), ]
            })
         }
         boot_data <- do.call(rbind, boot_list)

         # Calculate pseudo-observations within each bootstrap sample's strata
         pseudo_obs_list <- by(boot_data, boot_data[[if(is_stratified) strata_var else arm_var]], function(sub_data) {
            sub_data$pseudo_obs <- get_pseudo_obs(sub_data[[time_var]], sub_data[[status_var]], tau)
            sub_data
         })
         boot_data <- do.call(rbind, pseudo_obs_list)

         fit <- tryCatch(mgcv::gam(model_formula, data = boot_data), error = function(e) NULL)
         if (!is.null(fit)) {
            sfit <- mgcv::summary.gam(fit)
            p_table <- sfit$p.table
            matching_rows <- grep(test_term_pattern, rownames(p_table), fixed=TRUE)

            if (length(matching_rows) > 0) {
               p_val <- min(p_table[matching_rows, "Pr(>|t|)"], na.rm=TRUE)
               estimate <- mean(p_table[matching_rows, "Estimate"], na.rm=TRUE)
               std_error <- mean(p_table[matching_rows, "Std. Error"], na.rm=TRUE)
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
   group_label <- if(is_stratified) "/stratum" else "/arm"
   cat("--- Calculating Power/Sample Size (Method: Additive GAM for RMST) ---\n")
   message("Model: pseudo_obs ~ ", model_rhs)

   results_summary <- NULL

   if (!is.null(sample_sizes)) {
      # Power Calculation Mode
      sim_outputs <- lapply(sample_sizes, function(n) {
         cat(paste0("Simulating for n = ", n, group_label, "...\n"))
         run_power_sim(n)
      })
      power_values <- sapply(sim_outputs, `[[`, "power")
      results_df <- data.frame(N_per_Group = sample_sizes, Power = power_values)

      est <- na.omit(sim_outputs[[1]]$estimates)
      se <- na.omit(sim_outputs[[1]]$std_errors)

      if(length(est) > 1) {
         results_summary <- data.frame(
            Statistic = c("Mean RMST Difference", "Mean Standard Error", "95% CI Lower", "95% CI Upper"),
            Value = c(mean(est), mean(se, na.rm=TRUE), mean(est) - 1.96 * sd(est), mean(est) + 1.96 * sd(est))
         )
      }

      p <- ggplot2::ggplot(results_df, ggplot2::aes(x = N_per_Group, y = Power)) +
         ggplot2::geom_line(color = "#0072B2", linewidth = 1) +
         ggplot2::geom_point(color = "#0072B2", size = 3) +
         ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
         ggplot2::labs(title = "Power Curve: Additive GAM RMST Model",
                       x = paste0("Sample Size Per Arm", if(is_stratified) "/Stratum" else ""), y = "Estimated Power") +
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

            cat(paste0("  N = ", current_n, group_label, ", Calculated Power = ", round(calculated_power, 3), "\n"))

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
      results_df <- data.frame(Target_Power = sort(target_powers), Required_N_per_Group = required_n_values)

      p <- ggplot2::ggplot(na.omit(results_df), ggplot2::aes(x = Target_Power, y = Required_N_per_Group)) +
         ggplot2::geom_line(color = "#009E73", linewidth = 1) +
         ggplot2::geom_point(color = "#009E73", size = 3) +
         ggplot2::geom_vline(xintercept = 0.8, linetype = "dashed", color = "red") +
         ggplot2::labs(title = "Sample Size Curve: Additive GAM RMST Model",
                       x = "Target Power", y = paste0("Required Sample Size Per Arm", if(is_stratified) "/Stratum" else "")) +
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

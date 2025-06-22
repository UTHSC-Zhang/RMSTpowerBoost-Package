#' @title Analyze Power/Sample Size for a Multiplicative Stratified RMST Model
#' @description Performs power analysis based on the multiplicative model for RMST
#'   (Wang et al., 2019), designed for stratified trials.
#'
#' @note `status_var` should be `1`/`0`. `arm_var` should be `1`/`0`. `strata_var`
#'   is a mandatory argument.
#'
#' @param pilot_data A data.frame with pilot study data.
#' @param time_var,status_var,arm_var,strata_var Strings for column names.
#' @param sample_sizes Optional vector of sample sizes per stratum.
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
#' @references Wang, X., et al. (2019). *Statistics in Medicine*, 38, 5133-5145.
#'
#' @importFrom survival survfit Surv
#' @importFrom stats lm as.formula complete.cases integrate na.omit quantile stepfun sd
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline geom_vline labs theme_minimal ylim
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#'
#' @examples
#' # Create a small pilot dataset with strata
#' pilot_data_strat <- data.frame(
#'   time = rexp(120, rate = 0.15),
#'   status = rbinom(120, 1, 0.6),
#'   arm = rep(0:1, each = 60),
#'   # Create a stratification variable with 3 levels
#'   region = factor(rep(c("A", "B", "C"), each = 40)),
#'   age = rnorm(120, 55, 8)
#' )
#'
#' # --- Power Calculation Example ---
#' results_power_strat <- design_rmst_multiplicative_stratified(
#'   pilot_data = pilot_data_strat,
#'   time_var = "time",
#'   status_var = "status",
#'   arm_var = "arm",
#'   strata_var = "region",
#'   sample_sizes = c(50),
#'   tau = 10,
#'   n_sim = 10 # Low n_sim for example speed
#' )
#'
#' \donttest{
#' # --- Sample Size Calculation with Parallel Processing ---
#' results_n_strat <- design_rmst_multiplicative_stratified(
#'   pilot_data = pilot_data_strat,
#'   time_var = "time",
#'   status_var = "status",
#'   arm_var = "arm",
#'   strata_var = "region",
#'   target_powers = c(0.80),
#'   tau = 10,
#'   n_sim = 100,
#'   parallel.cores = 2
#' )
#' }
#'
#' @export
design_rmst_multiplicative_stratified <- function(pilot_data, time_var, status_var, arm_var, strata_var,
                                                  sample_sizes = NULL, target_powers = NULL,
                                                  linear_terms = NULL, tau, n_sim = 1000, alpha = 0.05,
                                                  parallel.cores = 1, patience = 5,
                                                  n_start = 50, n_step = 25, max_n_per_arm = 2000) {

   # --- Input Validation & Package Checks ---
   if (is.null(sample_sizes) && is.null(target_powers)) stop("Must provide 'sample_sizes' or 'target_powers'.")
   if (parallel.cores > 1) {
      if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
         stop("Packages 'future' and 'future.apply' are required for parallel processing.")
      }
   }

   # --- Data Preparation & Model Formula ---
   core_vars <- c(time_var, status_var, arm_var, strata_var)
   if (is.null(linear_terms)) {
      linear_terms <- setdiff(names(pilot_data), core_vars)
      if(length(linear_terms) > 0) message(paste("Using unspecified columns as linear terms:", paste(linear_terms, collapse=", ")))
   }
   all_vars <- c(core_vars, linear_terms)
   cleaned_pilot_data <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]
   pilot_groups <- split(cleaned_pilot_data, cleaned_pilot_data[[strata_var]])

   # --- Model Formula Construction ---
   all_terms <- c(strata_var, paste0(strata_var, ":", arm_var), linear_terms)
   model_rhs <- paste(all_terms[!sapply(all_terms, is.null)], collapse = " + ")
   model_formula <- as.formula(paste("log_pseudo_obs ~", model_rhs))
   test_term_pattern <- paste0(":", arm_var)

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

   # --- Internal Simulation Function ---
   run_power_sim <- function(n_per_stratum) {

      single_iteration <- function() {
         boot_list <- lapply(pilot_groups, function(group_df) group_df[sample(1:nrow(group_df), size = n_per_stratum, replace = TRUE), ])
         boot_data <- do.call(rbind, boot_list)

         pseudo_obs_list <- by(boot_data, boot_data[[strata_var]], function(sub_data) {
            sub_data$pseudo_obs <- get_pseudo_obs(sub_data[[time_var]], sub_data[[status_var]], tau)
            sub_data
         })
         boot_data <- do.call(rbind, pseudo_obs_list)

         boot_data$log_pseudo_obs <- log(boot_data$pseudo_obs)
         boot_data <- boot_data[is.finite(boot_data$log_pseudo_obs) & boot_data$pseudo_obs > 0,]

         if(nrow(boot_data) > (length(linear_terms) + 2*length(unique(boot_data[[strata_var]])))) {
            fit <- tryCatch(stats::lm(model_formula, data = boot_data), error = function(e) NULL)
            if (!is.null(fit)) {
               sfit <- summary(fit)
               coeffs <- sfit$coefficients
               matching_rows <- grep(test_term_pattern, rownames(coeffs), fixed=TRUE)
               if (length(matching_rows) > 0) {
                  p_vals <- coeffs[matching_rows, "Pr(>|t|)"]
                  p_val <- min(p_vals, na.rm=TRUE)
                  estimate <- mean(coeffs[matching_rows, "Estimate"], na.rm=TRUE)
                  std_error <- mean(coeffs[matching_rows, "Std. Error"], na.rm=TRUE)
               }
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
      valid_estimates <- estimates[is.finite(estimates)]
      return(list(power = power, estimates = exp(valid_estimates), std_errors = std_errors))
   }

   # --- Main Execution Logic ---
   cat("--- Calculating Power/Sample Size (Method: Multiplicative Stratified RMST Model) ---\n")
   message("Model: log(pseudo_obs) ~ ", model_rhs)

   results_summary <- NULL

   if (!is.null(sample_sizes)) {
      # Power Calculation Mode
      sim_outputs <- lapply(sample_sizes, function(n) {
         cat(paste0("Simulating for n = ", n, "/stratum...\n"))
         run_power_sim(n)
      })
      power_values <- sapply(sim_outputs, `[[`, "power")
      results_df <- data.frame(N_per_Stratum = sample_sizes, Power = power_values)

      est <- na.omit(sim_outputs[[1]]$estimates)
      if(length(est) > 1) {
         results_summary <- data.frame(
            Statistic = c("Mean RMST Ratio", "95% CI Lower", "95% CI Upper"),
            Value = c(mean(est), quantile(est, 0.025), quantile(est, 0.975))
         )
      }

      p <- ggplot2::ggplot(results_df, ggplot2::aes(x = N_per_Stratum, y = Power)) +
         ggplot2::geom_line(color = "#E69F00", linewidth = 1) +
         ggplot2::geom_point(color = "#E69F00", size = 3) +
         ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
         ggplot2::labs(title = "Power Curve: Multiplicative Stratified RMST Model",
                       x = "Sample Size Per Stratum", y = "Estimated Power") +
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

            cat(paste0("  N = ", current_n, "/stratum, Calculated Power = ", round(calculated_power, 3), "\n"))

            if (calculated_power >= target_power) {
               est <- na.omit(sim_output$estimates)
               if(length(est) > 1) {
                  results_summary <<- data.frame(
                     Statistic = c("Mean RMST Ratio", "95% CI Lower", "95% CI Upper"),
                     Value = c(mean(est), quantile(est, 0.025), quantile(est, 0.975))
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
      results_df <- data.frame(Target_Power = sort(target_powers), Required_N_per_Stratum = required_n_values)

      p <- ggplot2::ggplot(na.omit(results_df), ggplot2::aes(x = Target_Power, y = Required_N_per_Stratum)) +
         ggplot2::geom_line(color = "#009E73", linewidth = 1) +
         ggplot2::geom_point(color = "#009E73", size = 3) +
         ggplot2::geom_vline(xintercept = 0.8, linetype = "dashed", color = "red") +
         ggplot2::labs(title = "Sample Size Curve: Multiplicative Stratified RMST Model",
                       x = "Target Power", y = "Required Sample Size Per Stratum") +
         ggplot2::theme_minimal()
   }

   cat("\n--- Simulation Summary ---\n")
   if(!is.null(results_summary)) {
      print(knitr::kable(results_summary, caption = "Estimated Treatment Effect (RMST Ratio)"))
   } else {
      cat("No valid estimates were generated to create a summary.\n")
   }
   cat("\n")
   print(p)

   return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}

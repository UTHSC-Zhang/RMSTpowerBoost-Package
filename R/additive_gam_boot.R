# Shared Simulation Runner -------------------------------------------------

#' Internal setup shared by GAM.power.boot and GAM.ss.boot
#'
#' Prepares the pilot data and model formula once and returns the replicate
#' runner, so that the bootstrap procedure is written exactly once for both the
#' power and the sample-size functions.
#'
#' @return A list with the cleaned `pilot_data`, the `model_rhs` string, and
#'   `run_power_sim`, a function of the per-arm (or per-stratum) sample size
#'   that returns the simulated power and the per-replicate draws. It runs its
#'   replicates through `future.apply::future_lapply()`, so the caller sets the
#'   future plan.
#' @noRd
.gam_simulation_runner <- function(pilot_data, time_var, status_var, arm_var,
                                   strata_var, linear_terms, smooth_terms, L,
                                   n_sim, alpha) {
   is_stratified <- !is.null(strata_var)

   # --- Prepare data ---
   core_vars <- c(time_var, status_var, arm_var, strata_var)
   all_vars <- c(core_vars, linear_terms, smooth_terms)
   pilot_data <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]
   # A numeric stratum code must give stratum-specific intercepts, not a slope.
   if (is_stratified) pilot_data[[strata_var]] <- factor(pilot_data[[strata_var]])

   # --- Model Formula Construction ---
   smooth_part <- if (!is.null(smooth_terms)) paste0("s(", smooth_terms, ")") else NULL
   if (is_stratified) {
      pilot_groups <- split(pilot_data, pilot_data[[strata_var]])
      # Stratum-specific intercepts with a common treatment effect, matching the
      # single-effect hypothesis tested by the analytical stratified models.
      # (Testing the minimum p-value across stratum-by-arm interactions would
      # inflate the type-I error and overstate power.)
      all_terms <- c(strata_var, arm_var, smooth_part, linear_terms)
   } else {
      pilot_groups <- list(pilot_data)
      all_terms <- c(arm_var, smooth_part, linear_terms)
   }
   model_rhs <- paste(all_terms[!sapply(all_terms, is.null)], collapse = " + ")
   model_formula <- stats::as.formula(paste("pseudo_obs ~", model_rhs))
   # The treatment coefficient, named exactly. Matching on a prefix would also
   # pick up any term whose name starts with arm_var (e.g. "Age" for arm "A"),
   # and taking the minimum p-value over those rows would test the wrong
   # hypothesis.
   test_term <- paste0(arm_var, "1")

   # --- Jackknife pseudo-observations (shared helper) ---
   get_pseudo_obs <- .rmst_jackknife_pseudo_obs

   run_power_sim <- function(n_per_group) {
      sim_results_list <- future.apply::future_lapply(1:n_sim, function(j) {
         p_val <- NA_real_; est_val <- NA_real_; se_val <- NA_real_
         if (is_stratified) {
            boot_list <- lapply(pilot_groups, function(group_df) group_df[sample(1:nrow(group_df), size = n_per_group, replace = TRUE), ])
         } else {
            group_arms <- split(pilot_groups[[1]], pilot_groups[[1]][[arm_var]])
            boot_list <- lapply(group_arms, function(arm_df) arm_df[sample(1:nrow(arm_df), size = n_per_group, replace = TRUE), ])
         }
         boot_data <- do.call(rbind, boot_list)
         # Explicitly set arm to factor for consistent coefficient naming
         boot_data[[arm_var]] <- factor(boot_data[[arm_var]], levels = c(0, 1))
         grouping_var <- if(is_stratified) strata_var else arm_var
         pseudo_obs_list <- by(boot_data, boot_data[[grouping_var]], function(sub_data) {
            sub_data$pseudo_obs <- get_pseudo_obs(sub_data[[time_var]], sub_data[[status_var]], L)
            sub_data
         })
         boot_data <- do.call(rbind, pseudo_obs_list)
         fit <- tryCatch(mgcv::gam(model_formula, data = boot_data), error = function(e) NULL)
         if (!is.null(fit)) {
            p_table <- mgcv::summary.gam(fit)$p.table
            if (test_term %in% rownames(p_table)) {
               p_val <- p_table[test_term, "Pr(>|t|)"]
               est_val <- p_table[test_term, "Estimate"]
               se_val <- p_table[test_term, "Std. Error"]
            }
         }
         return(list(p_value = p_val, estimate = est_val, std_error = se_val))
      }, future.seed = TRUE)

      p_values <- vapply(sim_results_list, `[[`, "p_value", FUN.VALUE = numeric(1))
      estimates <- vapply(sim_results_list, `[[`, "estimate", FUN.VALUE = numeric(1))
      std_errors <- vapply(sim_results_list, `[[`, "std_error", FUN.VALUE = numeric(1))
      finite_p_values <- p_values[is.finite(p_values)]

      list(power = if (length(finite_p_values) > 0) mean(finite_p_values < alpha) else 0,
           estimates  = estimates,
           std_errors = std_errors,
           p_values   = p_values)
   }

   list(pilot_data = pilot_data, model_rhs = model_rhs, run_power_sim = run_power_sim)
}


# Power Calculations -------------------------------------------------------
#' @title Calculate Power for a Semiparametric Additive RMST Model via Simulation
#' @description Performs a power analysis for given sample sizes using a flexible,
#'   semiparametric additive model for the RMST based on pseudo-observations.
#'
#' @details
#' This function estimates power by bootstrap resampling the pilot data,
#' computing RMST pseudo-observations, fitting a semiparametric additive model,
#' and recording whether the treatment effect is significant. If `strata_var`
#' is supplied, resampling is carried out within strata and the model includes
#' stratum-specific intercepts with a common treatment effect (the same
#' single-effect hypothesis tested by the analytical stratified models). The
#' stratification variable is treated as a factor, so a numeric stratum code
#' still gives one intercept per stratum rather than a single slope. The model
#' formula can include linear terms and non-linear smooth terms (`s()`).
#'
#' The jackknife pseudo-observations are computed from the Kaplan-Meier curve
#' within each arm (within each stratum when `strata_var` is supplied), with the
#' area under the curve up to `L` computed exactly as a sum of rectangles.
#'
#' Power is estimated as the proportion of simulations with a treatment-effect
#' p-value below `alpha`. The treatment effect is the coefficient of `arm_var`
#' coded as a 0/1 factor, selected by its exact name, so a covariate whose name
#' begins with the arm variable's name is never mistaken for it. This
#' pseudo-observation approach models RMST directly without using the IPCW
#' estimating equations used elsewhere in the package.
#'
#' @note `status_var` should be `1` for an event, `0` for censored. `arm_var`
#'   should be `1` for treatment, `0` for control.
#'
#' @param pilot_data A `data.frame` with pilot study data.
#' @param time_var A character string for the time-to-event variable.
#' @param status_var A character string for the event status variable.
#' @param arm_var A character string for the treatment arm variable.
#' @param strata_var An optional character string for a stratification variable.
#' @param sample_sizes A numeric vector of sample sizes per arm/stratum.
#' @param linear_terms Optional character vector of covariates with a linear effect.
#' @param smooth_terms Optional character vector of covariates with a non-linear effect.
#' @param L The numeric truncation time for RMST.
#' @param n_sim Number of bootstrap simulations.
#' @param alpha The significance level.
#' @param parallel.cores Number of cores for parallel processing.
#' @param verbose Logical; if \code{TRUE}, emit progress messages. Default \code{FALSE}.
#'
#' @return A `list` containing:
#' \item{results_data}{A `data.frame` of sample sizes and corresponding estimated power.}
#' \item{results_plot}{A `ggplot` object visualizing the power curve.}
#' \item{results_summary}{A `data.frame` with a summary of the estimated treatment effect.}
#'
#' @importFrom survival survfit Surv
#' @importFrom stats stepfun integrate lm as.formula complete.cases na.omit sd quantile predict
#' @importFrom mgcv gam summary.gam
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline labs theme_minimal ylim
#' @importFrom knitr kable
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @export
#' @examples
#' \donttest{
#' pilot_df <- data.frame(
#'   time = rexp(100, 0.08),
#'   status = rbinom(100, 1, 0.7),
#'   arm = rep(0:1, each = 50),
#'   age = rnorm(100, 60, 10)
#' )
#' # Add a treatment effect
#' pilot_df$time[pilot_df$arm == 1] <- pilot_df$time[pilot_df$arm == 1] * 1.3
#'
#' power_results <- GAM.power.boot(
#'   pilot_data = pilot_df,
#'   time_var = "time",
#'   status_var = "status",
#'   arm_var = "arm",
#'   sample_sizes = c(100, 150),
#'   linear_terms = "age",
#'   L = 15,
#'   n_sim = 100, # Use more sims in practice, e.g., 1000
#'   parallel.cores = 2
#' )
#' print(power_results$results_data)
#' print(power_results$results_plot)
#' }
GAM.power.boot <- function(pilot_data, time_var, status_var, arm_var, strata_var = NULL,
                           sample_sizes, linear_terms = NULL, smooth_terms = NULL,
                           L, n_sim = 1000, alpha = 0.05,
                           parallel.cores = 1, verbose = FALSE) {
   start_time <- proc.time()
   if (is.null(sample_sizes)) stop("You must provide a numeric vector for 'sample_sizes'.")
   if (!is.numeric(sample_sizes) || any(!is.finite(sample_sizes)) ||
       any(sample_sizes < 2) || any(sample_sizes != floor(sample_sizes))) {
      stop("'sample_sizes' must contain integers greater than or equal to 2.")
   }
   if (parallel.cores > 1) {
      if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
         stop("Packages 'future' and 'future.apply' are required for parallel processing.")
      }
   }

   is_stratified <- !is.null(strata_var)
   runner <- .gam_simulation_runner(pilot_data, time_var, status_var, arm_var,
                                    strata_var, linear_terms, smooth_terms, L,
                                    n_sim, alpha)
   pilot_data <- runner$pilot_data

   group_label <- if(is_stratified) "/stratum" else "/arm"
   .rmst_verbose_message(verbose, "--- Calculating Power (Method: Additive GAM for RMST) ---")
   .rmst_verbose_message(verbose, "Model: pseudo_obs ~ ", runner$model_rhs)

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
      .rmst_verbose_message(verbose, "Simulating for n = ", n_per_group, group_label, "...")
      all_sim_outputs[[i]] <- runner$run_power_sim(n_per_group)
   }

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
      ggplot2::labs(title = "Power Curve: Additive GAM RMST Model",
                    x = paste0("Sample Size Per ", if(is_stratified) "Stratum" else "Arm"),
                    y = "Estimated Power") +
      ggplot2::ylim(0, 1) + ggplot2::theme_minimal()

   end_time <- proc.time()
   elapsed_time <- round((end_time - start_time)["elapsed"], 2)
   .rmst_verbose_message(verbose, "Total simulation time: ", elapsed_time, " seconds")
   if (isTRUE(verbose) && is.null(results_summary)) {
      .rmst_verbose_message(verbose, "No valid estimates were generated to create a summary.")
   }

   best_out <- all_sim_outputs[[which.max(sample_sizes)]]
   sim_draws <- data.frame(
      replicate  = seq_along(best_out$estimates),
      estimate   = best_out$estimates,
      std_error  = best_out$std_errors,
      p_value    = best_out$p_values
   )
   model_output <- list(
      coefficient_table   = NULL,
      treatment_effect    = results_summary,
      arm_specific_rmst   = NULL,
      variance_components = NULL,
      censoring_weights   = NULL,
      diagnostics         = list(n_used = nrow(pilot_data), n_events = sum(pilot_data[[status_var]]),
                                 n_sim = n_sim, convergence_ok = TRUE, singular_flag = FALSE),
      simulation_draws    = sim_draws
   )

   return(list(results_data = results_df, results_plot = p,
               results_summary = results_summary, model_output = model_output))
}


# Sample Size Search ------------------------------------------------------

#' @title Find Sample Size for a Semiparametric Additive RMST Model via Simulation
#' @description Performs an iterative sample size search to achieve a target power
#'   using a flexible, semiparametric additive model for the RMST.
#'
#' @details
#' This function increases the sample size step by step until the estimated
#' power reaches `target_power` or a stopping rule is met. At each step it runs
#' the same bootstrap procedure described in `GAM.power.boot`, including the
#' factor treatment of `strata_var`, the exact pseudo-observations, and the
#' exact selection of the treatment coefficient.
#'
#' @note This function's methodology is bootstrap-based.
#'
#' @param pilot_data A `data.frame` with pilot study data.
#' @param time_var A character string for the time-to-event variable.
#' @param status_var A character string for the event status variable.
#' @param arm_var A character string for the treatment arm variable.
#' @param strata_var An optional string for a stratification variable.
#' @param target_power A single numeric value for the target power.
#' @param linear_terms Optional character vector of covariates with a linear effect.
#' @param smooth_terms Optional character vector of covariates with a non-linear effect.
#' @param L The numeric truncation time for RMST.
#' @param n_sim Number of bootstrap simulations per search step.
#' @param alpha The significance level.
#' @param parallel.cores Number of cores for parallel processing.
#' @param patience Number of consecutive non-improving steps in the search before terminating.
#' @param n_start The starting sample size per arm/stratum for the search.
#' @param n_step The increment in sample size at each step of the search.
#' @param max_n_per_arm The maximum sample size per arm/stratum to search up to.
#' @param verbose Logical; if \code{TRUE}, emit progress messages. Default \code{FALSE}.
#'
#' @return A `list` containing:
#' \item{results_data}{A `data.frame` with the target power and required sample size.}
#' \item{results_plot}{A `ggplot` object of the search path.}
#' \item{results_summary}{A `data.frame` summarizing the estimated treatment effect.}
#'
#' @importFrom survival survfit Surv
#' @importFrom stats stepfun integrate lm as.formula complete.cases na.omit sd quantile predict
#' @importFrom mgcv gam summary.gam
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline geom_vline labs theme_minimal
#' @importFrom knitr kable
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @export
#' @examples
#' \donttest{
#' pilot_df_effect <- data.frame(
#'   time = c(stats::rexp(50, 0.1), stats::rexp(50, 0.04)), # Effect
#'   status = stats::rbinom(100, 1, 0.9),
#'   arm = rep(0:1, each = 50)
#' )
#'
#' ss_results <- GAM.ss.boot(
#'   pilot_data = pilot_df_effect,
#'   time_var = "time",
#'   status_var = "status",
#'   arm_var = "arm",
#'   target_power = 0.80,
#'   L = 15,
#'   n_sim = 100,      # Low n_sim for example
#'   n_start = 100,
#'   n_step = 50,
#'   patience = 2,
#'   parallel.cores = 2
#' )
#' print(ss_results$results_data)
#' print(ss_results$results_plot)
#' }
GAM.ss.boot <- function(pilot_data, time_var, status_var, arm_var, strata_var = NULL,
                        target_power,
                        linear_terms = NULL, smooth_terms = NULL,
                        L, n_sim = 1000, alpha = 0.05,
                        parallel.cores = 1, patience = 5,
                        n_start = 50, n_step = 25, max_n_per_arm = 2000,
                        verbose = FALSE) {

   start_time <- proc.time()
   if (is.null(target_power) || length(target_power) != 1 || !is.numeric(target_power)) {
      stop("You must provide a single numeric value for 'target_power'.")
   }
    if (!is.numeric(n_start) || length(n_start) != 1 || !is.finite(n_start) ||
        n_start < 2 || n_start != floor(n_start)) {
       stop("'n_start' must be a single integer greater than or equal to 2.")
    }
    if (!is.numeric(n_step) || length(n_step) != 1 || !is.finite(n_step) ||
        n_step < 1 || n_step != floor(n_step)) {
       stop("'n_step' must be a single positive integer.")
    }
    if (!is.numeric(max_n_per_arm) || length(max_n_per_arm) != 1 || !is.finite(max_n_per_arm) ||
        max_n_per_arm < 2 || max_n_per_arm != floor(max_n_per_arm)) {
       stop("'max_n_per_arm' must be a single integer greater than or equal to 2.")
    }
    if (max_n_per_arm < n_start) {
       stop("'max_n_per_arm' must be greater than or equal to 'n_start'.")
    }
   if (parallel.cores > 1) {
      if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
         stop("Packages 'future' and 'future.apply' are required for parallel processing.")
      }
   }

   is_stratified <- !is.null(strata_var)
   runner <- .gam_simulation_runner(pilot_data, time_var, status_var, arm_var,
                                    strata_var, linear_terms, smooth_terms, L,
                                    n_sim, alpha)
   pilot_data <- runner$pilot_data

   group_label <- if(is_stratified) "/stratum" else "/arm"
   .rmst_verbose_message(verbose, "--- Searching for Sample Size (Method: Additive GAM for RMST) ---")
   .rmst_verbose_message(verbose, "Model: pseudo_obs ~ ", runner$model_rhs)

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
      .rmst_verbose_message(verbose, "  N = ", current_n, group_label, ", calculating power...")
      sim_output <- runner$run_power_sim(current_n)

      calculated_power <- if(is.finite(sim_output$power)) sim_output$power else 0
      search_results[[as.character(current_n)]] <- calculated_power
      .rmst_verbose_message(verbose, "  Power = ", round(calculated_power, 3))

      if (calculated_power >= target_power) {
         .rmst_verbose_message(verbose, "Success: Target power reached at N = ", current_n, group_label, ".")
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
      ggplot2::geom_vline(xintercept = final_n, linetype = "dotted", color = "blue") +
      ggplot2::labs(title = "Sample Size Search Path: Additive GAM RMST Model",
                    x = paste0("Sample Size Per ", if(is_stratified) "Stratum" else "Arm"),
                    y = "Calculated Power") +
      ggplot2::theme_minimal()

   end_time <- proc.time()
   elapsed_time <- round((end_time - start_time)["elapsed"], 2)
   .rmst_verbose_message(verbose, "Total simulation time: ", elapsed_time, " seconds")
   if (isTRUE(verbose) && is.null(results_summary)) {
      .rmst_verbose_message(verbose, "No valid estimates were generated to create a summary.")
   }

   sim_draws <- if (!is.null(best_sim_output)) {
      data.frame(replicate = seq_along(best_sim_output$estimates),
                 estimate  = best_sim_output$estimates,
                 std_error = best_sim_output$std_errors,
                 p_value   = best_sim_output$p_values)
   } else NULL
   model_output <- list(
      coefficient_table   = NULL,
      treatment_effect    = results_summary,
      arm_specific_rmst   = NULL,
      variance_components = NULL,
      censoring_weights   = NULL,
      diagnostics         = list(n_used = nrow(pilot_data), n_events = sum(pilot_data[[status_var]]),
                                 n_sim = n_sim, convergence_ok = TRUE, singular_flag = FALSE),
      simulation_draws    = sim_draws
   )

   return(list(results_data = results_df, results_plot = p,
               results_summary = results_summary, model_output = model_output))
}


# estimation-helpers.R
# Internal shared machinery for the analytical power/sample-size engines.
# Each engine has a single .estimate_*_params() function (called by both its
# *.power.* and *.ss.* wrapper) and a .*_model_output() builder, so the
# statistical procedure is written exactly once per engine.

# --- Common data preparation ------------------------------------------------

#' Truncate follow-up at L and flag events / complete observations
#' @noRd
.rmst_prepare_df <- function(pilot_data, all_vars, time_var, status_var, L) {
   df <- pilot_data[stats::complete.cases(pilot_data[, all_vars, drop = FALSE]), ]
   df$Y_rmst <- pmin(df[[time_var]], L)
   df$is_event <- df[[status_var]] == 1
   df$is_complete <- df$is_event | df[[time_var]] >= L
   df
}

#' Wald power at a given total sample size
#' @noRd
.rmst_wald_power <- function(beta_effect, se_beta_n1, total_n, z_alpha) {
   se_final <- se_beta_n1 / sqrt(total_n)
   stats::pnorm((abs(beta_effect) / se_final) - z_alpha)
}

#' Iterative sample-size search on the analytic power formula
#' @noRd
.rmst_analytic_ss_search <- function(beta_effect, se_beta_n1, group_multiplier,
                                     target_power, alpha, n_start, n_step,
                                     max_n, group_label, verbose) {
   z_alpha <- stats::qnorm(1 - alpha / 2)
   current_n <- n_start
   search_path <- list()
   final_n <- NA_integer_
   while (current_n <= max_n) {
      calculated_power <- .rmst_wald_power(beta_effect, se_beta_n1,
                                           current_n * group_multiplier, z_alpha)
      if (!is.finite(calculated_power)) calculated_power <- 0
      search_path[[as.character(current_n)]] <- calculated_power
      .rmst_verbose_message(verbose, "  N = ", current_n, group_label,
                            ", calculated power = ", round(calculated_power, 3))
      if (calculated_power >= target_power) {
         final_n <- current_n
         break
      }
      current_n <- current_n + n_step
   }
   if (is.na(final_n)) {
      warning(paste("Target power", target_power, "not achieved by max N of", max_n),
              call. = FALSE)
      final_n <- max_n
   }
   list(final_n = final_n,
        search_path_df = data.frame(N = as.integer(names(search_path)),
                                    Power = as.numeric(unlist(search_path))))
}

# --- Shared plotting --------------------------------------------------------

#' Power-curve plot shared by the analytical engines
#' @noRd
.rmst_power_curve_plot <- function(results_df, x_var, color, title, subtitle, xlab) {
   ggplot2::ggplot(results_df, ggplot2::aes(x = .data[[x_var]], y = Power)) +
      ggplot2::geom_line(color = color, linewidth = 1) +
      ggplot2::geom_point(color = color, size = 3) +
      ggplot2::labs(title = title, subtitle = subtitle,
                    x = xlab, y = "Estimated Power") +
      ggplot2::ylim(0, 1) + ggplot2::theme_minimal()
}

#' Sample-size search-path plot shared by the analytical engines
#' @noRd
.rmst_ss_search_plot <- function(search_path_df, x_var, final_n, target_power,
                                 title, subtitle, xlab) {
   ggplot2::ggplot(stats::na.omit(search_path_df),
                   ggplot2::aes(x = .data[[x_var]], y = Power)) +
      ggplot2::geom_line(color = "#009E73", linewidth = 1) +
      ggplot2::geom_point(color = "#009E73", size = 3) +
      ggplot2::geom_hline(yintercept = target_power, linetype = "dashed", color = "red") +
      ggplot2::geom_vline(xintercept = final_n, linetype = "dotted", color = "blue") +
      ggplot2::labs(title = title, subtitle = subtitle,
                    x = xlab, y = "Calculated Power") +
      ggplot2::theme_minimal()
}

# --- Stratified-Cox IPCW ----------------------------------------------------

#' Fit a censoring Cox model, turning non-convergence into a clear error
#' @noRd
.fit_ms_censoring_model <- function(cens_formula, df) {
   tryCatch(
      withCallingHandlers(
         survival::coxph(cens_formula, data = df, ties = "breslow"),
         warning = function(w) {
            msg <- conditionMessage(w)
            if (grepl("did not converge|ran out of iterations|infinite", msg, ignore.case = TRUE)) {
               stop("Not enough data points to fit the censoring model: the Cox model for censoring did not converge.", call. = FALSE)
            }
         }
      ),
      error = function(e) {
         stop("Not enough data points to fit the censoring model: ", conditionMessage(e), call. = FALSE)
      }
   )
}

#' Per-subject baseline cumulative censoring hazard at Y_rmst, by stratum
#'
#' basehaz() labels strata either as "level" or "strata_var=level" depending
#' on the survival package version; accept both and fail loudly if neither
#' matches, since silently missing strata would zero out the censoring hazard.
#' @noRd
.strat_cum_hazard <- function(df, bh_cens, strata_var) {
   H_cens <- numeric(nrow(df))
   bh_labels <- if (is.null(bh_cens$strata)) NULL else as.character(bh_cens$strata)
   for (st in unique(df[[strata_var]])) {
      is_stratum <- df[[strata_var]] == st
      is_bh_stratum <- if (is.null(bh_labels)) rep(TRUE, nrow(bh_cens)) else
         bh_labels == as.character(st) | bh_labels == paste0(strata_var, "=", st)
      if (sum(is_bh_stratum) == 0) {
         stop("Could not locate the baseline censoring hazard for stratum '", st,
              "'; IPCW weights cannot be computed.", call. = FALSE)
      }
      H_cens[is_stratum] <- stats::stepfun(
         bh_cens$time[is_bh_stratum],
         c(0, bh_cens$hazard[is_bh_stratum]))(df$Y_rmst[is_stratum])
   }
   H_cens
}

#' Raw IPCW weights W = exp(Lambda_C(Y)) from a stratified Cox censoring model
#' @noRd
.ipcw_stratified_cox_weights <- function(df, time_var, status_var, covariates,
                                         strata_var, strict = FALSE) {
   cens_formula <- stats::as.formula(paste0(
      "survival::Surv(", time_var, ", ", status_var, " == 0) ~ ",
      paste(covariates, collapse = " + "),
      " + survival::strata(", strata_var, ")"))
   fit_cens <- if (strict) .fit_ms_censoring_model(cens_formula, df)
               else survival::coxph(cens_formula, data = df, ties = "breslow")
   bh_cens <- survival::basehaz(fit_cens, centered = FALSE)
   H_cens <- .strat_cum_hazard(df, bh_cens, strata_var)
   exp(H_cens * exp(stats::predict(fit_cens, newdata = df, type = "lp", reference = "zero")))
}

#' Cap weights at their 99th percentile and zero out non-finite values
#' @noRd
.cap_weights <- function(weights, positive_only = TRUE) {
   weight_cap <- NA_real_
   keep <- if (positive_only) is.finite(weights) & weights > 0 else is.finite(weights)
   finite_weights <- weights[keep]
   if (length(finite_weights) > 0) {
      weight_cap <- stats::quantile(finite_weights, probs = 0.99, na.rm = TRUE)
      weights[weights > weight_cap] <- weight_cap
   }
   weights[!is.finite(weights)] <- 0
   list(weights = weights, cap = weight_cap)
}

# --- Additive stratified engine (Zhang & Schaubel 2024) ----------------------

#' Estimate the stratified additive RMST model from pilot data
#' @noRd
.estimate_additive_params <- function(pilot_data, time_var, status_var, arm_var,
                                      strata_var, linear_terms, L) {
   covariates <- c(arm_var, linear_terms)
   all_vars <- c(time_var, status_var, strata_var, covariates)
   df <- .rmst_prepare_df(pilot_data, all_vars, time_var, status_var, L)
   n_pilot <- nrow(df)

   df$weights <- .ipcw_stratified_cox_weights(df, time_var, status_var,
                                              covariates, strata_var)
   df$weights[!df$is_complete] <- 0
   capped <- .cap_weights(df$weights, positive_only = TRUE)
   df$weights <- capped$weights
   weight_cap <- capped$cap

   # Stratum-centering: weighted (W * Delta_Y) means of the outcome and covariates
   vars_to_center <- c("Y_rmst", covariates)
   stratum_means <- df %>%
      dplyr::group_by(.data[[strata_var]]) %>%
      dplyr::summarise(
         dplyr::across(
            dplyr::all_of(vars_to_center),
            ~ weighted.mean(.x, w = weights, na.rm = TRUE)
         ),
         .groups = 'drop'
      )
   names(stratum_means) <- c(strata_var, paste0(vars_to_center, "_mean"))

   df_centered <- df %>% dplyr::left_join(stratum_means, by = strata_var)
   for (cov in vars_to_center) {
      df_centered[[paste0(cov, "_tilde")]] <- df_centered[[cov]] - df_centered[[paste0(cov, "_mean")]]
   }

   Z_tilde <- as.matrix(df_centered[, paste0(covariates, "_tilde")])
   Y_tilde <- df_centered[["Y_rmst_tilde"]]
   W <- df_centered$weights

   A_hat_num <- crossprod(Z_tilde * sqrt(W))
   dimnames(A_hat_num) <- list(covariates, covariates)
   A_hat <- A_hat_num / n_pilot

   A_hat_inv <- tryCatch({
      solve(A_hat)
   }, error = function(e) {
      stop("The covariate matrix (A_hat) is singular and cannot be inverted.\nThis may be caused by a lack of variation in the covariates among subjects with an event within one or more strata.\nPlease inspect the pilot data for issues like perfect separation.", call. = FALSE)
   })

   beta_hat <- (A_hat_inv / n_pilot) %*% (t(Z_tilde * W) %*% Y_tilde)
   rownames(beta_hat) <- covariates
   beta_effect <- beta_hat[arm_var, 1]

   mu0_hats <- stratum_means %>%
      dplyr::mutate(
         Z_matrix = as.matrix(dplyr::select(., dplyr::all_of(paste0(covariates, "_mean")))),
         mu0_hat = .data[["Y_rmst_mean"]] - Z_matrix %*% beta_hat
      ) %>%
      dplyr::select(dplyr::all_of(strata_var), mu0_hat)

   df_final <- df_centered %>% dplyr::left_join(mu0_hats, by = strata_var)
   Z_matrix <- as.matrix(df_final[, covariates])
   df_final$residuals <- df_final$Y_rmst - (df_final$mu0_hat + as.vector(Z_matrix %*% beta_hat))

   epsilon <- apply(Z_tilde, 2, function(z_col) z_col * W * df_final$residuals)
   B_hat <- crossprod(epsilon) / n_pilot
   dimnames(B_hat) <- list(covariates, covariates)

   V_hat_n <- A_hat_inv %*% B_hat %*% t(A_hat_inv)

   list(df = df, covariates = covariates, n_pilot = n_pilot,
        n_strata = length(unique(df[[strata_var]])),
        beta_hat = beta_hat, beta_effect = beta_effect,
        mu0_hats = mu0_hats,
        A_hat = A_hat, B_hat = B_hat, V_hat_n = V_hat_n,
        se_beta_n1 = sqrt(V_hat_n[arm_var, arm_var]),
        weight_cap = weight_cap)
}

#' Model-output summary for the stratified additive engine
#' @noRd
.additive_model_output <- function(est, arm_var, strata_var) {
   se_all  <- sqrt(diag(est$V_hat_n) / est$n_pilot)
   z_vals  <- as.numeric(est$beta_hat[, 1]) / se_all
   p_vals  <- 2 * stats::pnorm(-abs(z_vals))
   coef_tbl <- data.frame(
      term      = est$covariates,
      estimate  = as.numeric(est$beta_hat[, 1]),
      std_error = se_all,
      ci_lower  = as.numeric(est$beta_hat[, 1]) - 1.96 * se_all,
      ci_upper  = as.numeric(est$beta_hat[, 1]) + 1.96 * se_all,
      test_stat = z_vals,
      p_value   = p_vals,
      row.names = NULL,
      stringsAsFactors = FALSE
   )
   se_arm_pilot <- est$se_beta_n1 / sqrt(est$n_pilot)
   trt_eff <- data.frame(
      estimand  = "RMST Difference (additive)",
      estimate  = est$beta_effect,
      std_error = se_arm_pilot,
      ci_lower  = est$beta_effect - 1.96 * se_arm_pilot,
      ci_upper  = est$beta_effect + 1.96 * se_arm_pilot,
      stringsAsFactors = FALSE
   )
   strat_col <- est$mu0_hats[[strata_var]]
   mu0_vec   <- est$mu0_hats$mu0_hat
   arm_rmst <- do.call(rbind, c(
      lapply(seq_along(strat_col), function(i)
         data.frame(arm = 0, stratum = as.character(strat_col[i]),
                    rmst_estimate = mu0_vec[i],
                    std_error = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                    scale = "original", stringsAsFactors = FALSE)),
      lapply(seq_along(strat_col), function(i)
         data.frame(arm = 1, stratum = as.character(strat_col[i]),
                    rmst_estimate = mu0_vec[i] + est$beta_effect,
                    std_error = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                    scale = "original", stringsAsFactors = FALSE))
   ))
   capped_frac <- if (is.finite(est$weight_cap))
      mean(est$df$weights[est$df$is_complete] >= est$weight_cap, na.rm = TRUE) else NA_real_
   list(
      coefficient_table   = coef_tbl,
      treatment_effect    = trt_eff,
      arm_specific_rmst   = arm_rmst,
      variance_components = list(A_hat = est$A_hat, B_hat = est$B_hat,
                                 V_hat_n = est$V_hat_n, se_effect_n1 = est$se_beta_n1),
      censoring_weights   = list(
         raw_summary     = stats::quantile(est$df$weights, c(0, .25, .5, .75, .99, 1), na.rm = TRUE),
         cap_value       = est$weight_cap,
         capped_fraction = capped_frac),
      diagnostics         = list(n_used = sum(est$df$is_complete), n_events = sum(est$df$is_event),
                                 n_complete = sum(est$df$is_complete),
                                 convergence_ok = TRUE, singular_flag = FALSE),
      simulation_draws    = NULL
   )
}

# --- Multiplicative stratified engine (Wang et al. 2019, log-linear WLS) -----

#' Estimate the stratified multiplicative RMST working model from pilot data
#' @noRd
.estimate_ms_params <- function(pilot_data, time_var, status_var, arm_var,
                                strata_var, linear_terms, L) {
   covariates <- c(arm_var, linear_terms)
   all_vars <- c(time_var, status_var, arm_var, strata_var, linear_terms)
   df <- .rmst_prepare_df(pilot_data, all_vars, time_var, status_var, L)
   n_pilot <- nrow(df)

   df$weights <- .ipcw_stratified_cox_weights(df, time_var, status_var,
                                              covariates, strata_var, strict = TRUE)
   capped <- .cap_weights(df$weights, positive_only = FALSE)
   df$weights <- capped$weights
   weight_cap <- capped$cap
   df$w_delta <- df$weights * df$is_complete
   df$w_delta[is.na(df$w_delta)] <- 0

   # Log-linear working model fitted by weighted least squares
   fit_data <- df[df$w_delta > 0 & df$Y_rmst > 0, ]
   fit_weights <- fit_data$w_delta
   log_model_formula <- stats::as.formula(
      paste("log(Y_rmst) ~", paste(covariates, collapse = " + "), "+", strata_var))
   if (nrow(fit_data) < (length(covariates) + length(unique(fit_data[[strata_var]])))) {
      stop("Not enough data points to fit the approximation model after filtering.", call. = FALSE)
   }
   fit_log_lm <- stats::lm(log_model_formula, data = fit_data, weights = fit_weights)

   beta_all <- stats::coef(fit_log_lm)
   beta_effect <- beta_all[[arm_var]]

   # Robust sandwich variance A^-1 B A^-1 of the WLS estimating equations,
   # consistent with the other analytical engines (columns with NA coefficients
   # from a rank-deficient fit are dropped).
   X <- stats::model.matrix(fit_log_lm)
   ok <- !is.na(beta_all)
   X <- X[, ok, drop = FALSE]
   r <- stats::residuals(fit_log_lm)

   A_hat <- crossprod(X * sqrt(fit_weights)) / n_pilot
   A_hat_inv <- tryCatch({
      solve(A_hat)
   }, error = function(e) {
      stop("The design matrix (A_hat) is singular and cannot be inverted. This may be caused by a lack of variation in the covariates within one or more strata.", call. = FALSE)
   })
   epsilon <- X * (fit_weights * r)
   B_hat <- crossprod(epsilon) / n_pilot
   V_hat_n <- A_hat_inv %*% B_hat %*% t(A_hat_inv)

   list(df = df, fit_data = fit_data, fit_log_lm = fit_log_lm,
        covariates = covariates, n_pilot = n_pilot,
        n_strata = length(unique(df[[strata_var]])),
        beta_all = beta_all[ok], beta_effect = beta_effect,
        A_hat = A_hat, B_hat = B_hat, V_hat_n = V_hat_n,
        se_beta_n1 = sqrt(V_hat_n[arm_var, arm_var]),
        weight_cap = weight_cap)
}

#' Model-output summary for the multiplicative engine (robust SEs)
#' @noRd
.ms_model_output <- function(est, arm_var) {
   se_all <- sqrt(diag(est$V_hat_n) / est$n_pilot)
   z_vals <- as.numeric(est$beta_all) / se_all
   p_vals <- 2 * stats::pnorm(-abs(z_vals))
   coef_tbl <- data.frame(
      term      = names(est$beta_all),
      estimate  = as.numeric(est$beta_all),
      std_error = se_all,
      ci_lower  = as.numeric(est$beta_all) - 1.96 * se_all,
      ci_upper  = as.numeric(est$beta_all) + 1.96 * se_all,
      test_stat = z_vals,
      p_value   = p_vals,
      row.names = NULL,
      stringsAsFactors = FALSE
   )
   se_arm_pilot <- se_all[names(est$beta_all) == arm_var]
   trt_eff <- data.frame(
      estimand  = c("log(RMST Ratio)", "RMST Ratio"),
      estimate  = c(est$beta_effect, exp(est$beta_effect)),
      std_error = c(se_arm_pilot, NA_real_),
      ci_lower  = c(est$beta_effect - 1.96 * se_arm_pilot,
                    exp(est$beta_effect - 1.96 * se_arm_pilot)),
      ci_upper  = c(est$beta_effect + 1.96 * se_arm_pilot,
                    exp(est$beta_effect + 1.96 * se_arm_pilot)),
      stringsAsFactors = FALSE
   )
   pred_log <- stats::predict(est$fit_log_lm, newdata = est$fit_data)
   arms <- sort(unique(est$fit_data[[arm_var]]))
   arm_rmst <- do.call(rbind, lapply(arms, function(a) {
      idx <- est$fit_data[[arm_var]] == a
      mu_log  <- mean(pred_log[idx], na.rm = TRUE)
      mu_orig <- mean(exp(pred_log[idx]), na.rm = TRUE)
      rbind(
         data.frame(arm = a, rmst_estimate = mu_log,  std_error = NA_real_,
                    ci_lower = NA_real_, ci_upper = NA_real_,
                    scale = "log",      stringsAsFactors = FALSE),
         data.frame(arm = a, rmst_estimate = mu_orig, std_error = NA_real_,
                    ci_lower = NA_real_, ci_upper = NA_real_,
                    scale = "original", stringsAsFactors = FALSE)
      )
   }))
   capped_frac <- if (is.finite(est$weight_cap))
      mean(est$df$weights[est$df$is_complete] >= est$weight_cap, na.rm = TRUE) else NA_real_
   full_rank <- isTRUE(est$fit_log_lm$rank == length(stats::coef(est$fit_log_lm)))
   list(
      coefficient_table   = coef_tbl,
      treatment_effect    = trt_eff,
      arm_specific_rmst   = arm_rmst,
      variance_components = list(A_hat = est$A_hat, B_hat = est$B_hat,
                                 V_hat_n = est$V_hat_n, se_effect_n1 = est$se_beta_n1),
      censoring_weights   = list(
         raw_summary     = stats::quantile(est$df$weights, c(0, .25, .5, .75, .99, 1), na.rm = TRUE),
         cap_value       = est$weight_cap,
         capped_fraction = capped_frac),
      diagnostics         = list(n_used = nrow(est$fit_data), n_events = sum(est$df$is_event),
                                 n_complete = sum(est$df$is_complete),
                                 convergence_ok = full_rank, singular_flag = !full_rank),
      simulation_draws    = NULL
   )
}

# --- Linear IPCW engine (Tian et al. 2014) -----------------------------------

#' Estimate the linear IPCW RMST model from pilot data
#' @noRd
.estimate_linear_params <- function(pilot_data, time_var, status_var, arm_var,
                                    linear_terms, L, verbose = FALSE) {
   core_vars <- c(time_var, status_var, arm_var)
   all_vars <- c(core_vars, linear_terms)
   df <- .rmst_prepare_df(pilot_data, all_vars, time_var, status_var, L)
   n_pilot <- nrow(df)

   factor_arm_str <- paste0("factor(", arm_var, ")")
   model_rhs <- paste(c(factor_arm_str, linear_terms), collapse = " + ")
   model_formula <- stats::as.formula(paste("Y_rmst ~", model_rhs))
   .rmst_verbose_message(verbose, "Model: Y_rmst ~ ", model_rhs)

   df$is_censored <- df[[status_var]] == 0

   # Censoring model: marginal Kaplan-Meier for G(t)
   cens_formula <- stats::as.formula(paste0("survival::Surv(", time_var, ", is_censored) ~ 1"))
   cens_fit <- survival::survfit(cens_formula, data = df)
   cens_surv_prob <- stats::stepfun(cens_fit$time, c(1, cens_fit$surv))(df$Y_rmst)
   df$weights <- df$is_complete / cens_surv_prob

   capped <- .cap_weights(df$weights, positive_only = TRUE)
   df$weights <- capped$weights
   weight_cap <- capped$cap

   fit_data <- df[df$weights > 0, ]
   fit_weights <- fit_data$weights
   if (length(unique(fit_data[[arm_var]])) < 2) {
      stop("Pilot data contains events in only one arm after filtering.", call. = FALSE)
   }

   fit_lm <- stats::lm(model_formula, data = fit_data, weights = fit_weights)
   beta_hat <- stats::coef(fit_lm)
   arm_pattern <- paste0("^factor\\(", arm_var, "\\)1$")
   arm_coeff_name <- names(beta_hat)[grep(arm_pattern, names(beta_hat))]
   if (length(arm_coeff_name) == 0) stop("Could not find treatment effect coefficient.")
   beta_effect <- beta_hat[arm_coeff_name]

   X <- stats::model.matrix(model_formula, data = df)
   A_hat <- crossprod(X * sqrt(df$weights), X * sqrt(df$weights)) / n_pilot
   A_hat_inv <- tryCatch({
      solve(A_hat)
   }, error = function(e) {
      stop("The design matrix (A_hat) is singular. This can happen with small pilot datasets or perfect separation.", call. = FALSE)
   })

   df$predicted_rmst <- stats::predict(fit_lm, newdata = df)
   residuals <- df$Y_rmst - df$predicted_rmst
   epsilon <- X * residuals * df$weights
   epsilon[is.na(epsilon)] <- 0
   B_hat <- crossprod(epsilon) / n_pilot

   V_hat_n <- A_hat_inv %*% B_hat %*% t(A_hat_inv)

   list(df = df, fit_data = fit_data, fit_weights = fit_weights, fit_lm = fit_lm,
        n_pilot = n_pilot, arm_coeff_name = arm_coeff_name, beta_effect = beta_effect,
        A_hat = A_hat, B_hat = B_hat, V_hat_n = V_hat_n,
        se_beta_n1 = sqrt(V_hat_n[arm_coeff_name, arm_coeff_name]),
        weight_cap = weight_cap)
}

#' Model-output summary for the linear IPCW engine
#' @noRd
.linear_model_output <- function(est, arm_var) {
   cs <- summary(est$fit_lm)$coefficients
   coef_tbl <- data.frame(
      term      = rownames(cs),
      estimate  = cs[, 1L],
      std_error = cs[, 2L],
      ci_lower  = cs[, 1L] - 1.96 * cs[, 2L],
      ci_upper  = cs[, 1L] + 1.96 * cs[, 2L],
      test_stat = cs[, 3L],
      p_value   = cs[, 4L],
      row.names = NULL,
      stringsAsFactors = FALSE
   )
   se_arm_pilot <- est$se_beta_n1 / sqrt(est$n_pilot)
   trt_eff <- data.frame(
      estimand  = "RMST Difference",
      estimate  = est$beta_effect,
      std_error = se_arm_pilot,
      ci_lower  = est$beta_effect - 1.96 * se_arm_pilot,
      ci_upper  = est$beta_effect + 1.96 * se_arm_pilot,
      stringsAsFactors = FALSE
   )
   arms <- sort(unique(est$fit_data[[arm_var]]))
   arm_rmst <- do.call(rbind, lapply(arms, function(a) {
      idx   <- est$fit_data[[arm_var]] == a
      sub_w <- est$fit_weights[idx]
      sub_y <- est$fit_data$Y_rmst[idx]
      mu    <- if (sum(sub_w) > 0) stats::weighted.mean(sub_y, sub_w) else mean(sub_y)
      data.frame(arm = a, rmst_estimate = mu, std_error = NA_real_,
                 ci_lower = NA_real_, ci_upper = NA_real_,
                 scale = "original", stringsAsFactors = FALSE)
   }))
   capped_frac <- if (is.finite(est$weight_cap))
      mean(est$df$weights[est$df$is_complete] >= est$weight_cap, na.rm = TRUE) else NA_real_
   list(
      coefficient_table   = coef_tbl,
      treatment_effect    = trt_eff,
      arm_specific_rmst   = arm_rmst,
      variance_components = list(A_hat = est$A_hat, B_hat = est$B_hat,
                                 V_hat_n = est$V_hat_n, se_effect_n1 = est$se_beta_n1),
      censoring_weights   = list(
         raw_summary     = stats::quantile(est$df$weights, c(0, .25, .5, .75, .99, 1), na.rm = TRUE),
         cap_value       = est$weight_cap,
         capped_fraction = capped_frac),
      diagnostics         = list(n_used = nrow(est$fit_data), n_events = sum(est$df$is_event),
                                 n_complete = sum(est$df$is_complete),
                                 convergence_ok = TRUE, singular_flag = FALSE),
      simulation_draws    = NULL
   )
}

# --- Covariate-dependent censoring engine (Wang & Schaubel 2018) --------------

#' Estimate the DC-IPCW RMST model from pilot data
#' @noRd
.estimate_dc_params <- function(pilot_data, time_var, status_var, arm_var,
                                linear_terms, L) {
   covariates <- c(arm_var, linear_terms)
   all_vars <- unique(c(time_var, status_var, covariates))
   df <- .rmst_prepare_df(pilot_data, all_vars, time_var, status_var, L)
   n_pilot <- nrow(df)
   if (n_pilot < 10) stop("Too few complete cases after filtering.", call. = FALSE)

   # Censoring model: covariate-dependent only (no arm, no competing risks)
   cens_rhs <- if (is.null(linear_terms) || length(linear_terms) == 0) "1" else paste(linear_terms, collapse = " + ")
   cens_formula <- stats::as.formula(
      paste0("survival::Surv(", time_var, ", ", status_var, "==0) ~ ", cens_rhs)
   )
   fit_cens <- survival::coxph(cens_formula, data = df, ties = "breslow")

   bh <- survival::basehaz(fit_cens, centered = FALSE)
   H0_step <- stats::stepfun(bh$time, c(0, bh$hazard))
   lp <- if (cens_rhs == "1") 0 else stats::predict(fit_cens, newdata = df, type = "lp")
   Hc <- H0_step(df$Y_rmst) * exp(lp)
   Ghat <- exp(-Hc)

   # IPCW for subjects whose truncated RMST outcome is observed (stabilize)
   eps <- 1e-6
   w <- df$is_complete / pmax(Ghat, eps)
   w[!is.finite(w)] <- 0
   capped <- .cap_weights(w, positive_only = TRUE)
   df$w <- capped$weights
   cap <- capped$cap

   model_rhs <- paste(c(arm_var, linear_terms), collapse = " + ")
   model_formula <- stats::as.formula(paste("Y_rmst ~", model_rhs))
   fit_data <- df[df$w > 0, ]
   fit_wls <- stats::lm(model_formula, data = fit_data, weights = fit_data$w)
   beta <- stats::coef(fit_wls)
   if (!(arm_var %in% names(beta))) stop("Arm coefficient not found in model.", call. = FALSE)
   beta_arm <- beta[arm_var]

   X <- stats::model.matrix(model_formula, data = df)
   r <- df$Y_rmst - as.numeric(X %*% beta)
   W <- df$w

   A_hat <- crossprod(X, X * W) / n_pilot
   meat <- X * (W * r)
   B_hat <- crossprod(meat) / n_pilot

   Ainv <- solve(A_hat)
   Vn <- Ainv %*% B_hat %*% Ainv
   var_beta_n1 <- as.numeric(Vn[colnames(X) == arm_var, colnames(X) == arm_var])

   list(df = df, fit_data = fit_data, fit_wls = fit_wls,
        n_pilot = n_pilot, beta_effect = beta_arm,
        A_hat = A_hat, B_hat = B_hat, V_hat_n = Vn,
        se_beta_n1 = sqrt(var_beta_n1),
        weight_cap = cap)
}

#' Model-output summary for the DC-IPCW engine
#' @noRd
.dc_model_output <- function(est, arm_var) {
   cs <- summary(est$fit_wls)$coefficients
   coef_tbl <- data.frame(
      term      = rownames(cs),
      estimate  = cs[, 1L],
      std_error = cs[, 2L],
      ci_lower  = cs[, 1L] - 1.96 * cs[, 2L],
      ci_upper  = cs[, 1L] + 1.96 * cs[, 2L],
      test_stat = cs[, 3L],
      p_value   = cs[, 4L],
      row.names = NULL,
      stringsAsFactors = FALSE
   )
   se_arm_pilot <- est$se_beta_n1 / sqrt(est$n_pilot)
   trt_eff <- data.frame(
      estimand  = "RMST Difference (DC-IPCW)",
      estimate  = est$beta_effect,
      std_error = se_arm_pilot,
      ci_lower  = est$beta_effect - 1.96 * se_arm_pilot,
      ci_upper  = est$beta_effect + 1.96 * se_arm_pilot,
      stringsAsFactors = FALSE
   )
   arms <- sort(unique(est$df[[arm_var]]))
   arm_rmst <- do.call(rbind, lapply(arms, function(a) {
      idx <- est$df[[arm_var]] == a
      w_a <- est$df$w[idx]
      y_a <- est$df$Y_rmst[idx]
      mu  <- if (sum(w_a) > 0) stats::weighted.mean(y_a, w_a) else mean(y_a)
      data.frame(arm = a, rmst_estimate = mu, std_error = NA_real_,
                 ci_lower = NA_real_, ci_upper = NA_real_,
                 scale = "original", stringsAsFactors = FALSE)
   }))
   capped_frac <- if (is.finite(est$weight_cap))
      mean(est$df$w[est$df$w > 0] >= est$weight_cap, na.rm = TRUE) else NA_real_
   list(
      coefficient_table   = coef_tbl,
      treatment_effect    = trt_eff,
      arm_specific_rmst   = arm_rmst,
      variance_components = list(A_hat = est$A_hat, B_hat = est$B_hat,
                                 V_hat_n = est$V_hat_n, se_effect_n1 = est$se_beta_n1),
      censoring_weights   = list(
         raw_summary     = stats::quantile(est$df$w, c(0, .25, .5, .75, .99, 1), na.rm = TRUE),
         cap_value       = est$weight_cap,
         capped_fraction = capped_frac),
      diagnostics         = list(n_used = nrow(est$fit_data), n_events = sum(est$df$is_event),
                                 n_complete = sum(est$df$is_complete),
                                 convergence_ok = TRUE, singular_flag = FALSE),
      simulation_draws    = NULL
   )
}

# --- Jackknife pseudo-observations (shared by the bootstrap engines) ----------

#' Leave-one-out jackknife pseudo-observations for the RMST at L
#' @noRd
.rmst_jackknife_pseudo_obs <- function(time, status, L) {
   n <- length(time)
   if (n == 0) return(numeric(0))
   km_fit_full <- survival::survfit(survival::Surv(time, status) ~ 1)
   km_step_full <- stats::stepfun(km_fit_full$time, c(1, km_fit_full$surv))
   rmst_full <- tryCatch(stats::integrate(km_step_full, 0, L, subdivisions = 2000, stop.on.error = FALSE)$value, error = function(e) 0)
   rmst_loo <- vapply(seq_len(n), function(i) {
      if (n > 1) {
         km_fit_loo <- survival::survfit(survival::Surv(time[-i], status[-i]) ~ 1)
         km_step_loo <- stats::stepfun(km_fit_loo$time, c(1, km_fit_loo$surv))
         tryCatch(stats::integrate(km_step_loo, 0, L, subdivisions = 2000, stop.on.error = FALSE)$value, error = function(e) 0)
      } else { 0 }
   }, FUN.VALUE = numeric(1))
   n * rmst_full - (n - 1) * rmst_loo
}

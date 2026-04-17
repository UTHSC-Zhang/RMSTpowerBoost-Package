# rmst-interface.R â€” High-level wrappers rmst.power() and rmst.ss()
# with formula interface, routing logic, and S3 methods.

# â”€â”€ Internal helpers â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

#' Parse a Surv-based formula into components
#' @noRd
.parse_rmst_formula <- function(formula) {
  lhs <- formula[[2L]]
  if (!is.call(lhs) || !identical(tail(as.character(lhs[[1L]]), 1L), "Surv"))
    stop("LHS of formula must be Surv(time, status).", call. = FALSE)
  tt     <- stats::terms(formula)
  labels <- attr(tt, "term.labels")
  is_sm  <- grepl("^s\\(", labels)
  list(
    time_var     = as.character(lhs[[2L]]),
    status_var   = as.character(lhs[[3L]]),
    linear_terms = if (any(!is_sm)) labels[!is_sm] else NULL,
    smooth_terms = if (any(is_sm))
                     sub("^s\\(([^,)]+).*\\)$", "\\1", labels[is_sm]) else NULL
  )
}

#' Resolve strata argument to a character column name (or NULL)
#' @noRd
.resolve_strata <- function(strata) {
  if (is.null(strata))       return(NULL)
  if (inherits(strata, "formula")) {
    lbs <- attr(stats::terms(strata), "term.labels")
    if (length(lbs) != 1L)
      stop("strata formula must reference exactly one column, e.g. ~stratum.", call. = FALSE)
    return(lbs)
  }
  if (is.character(strata) && length(strata) == 1L) return(strata)
  stop("strata must be NULL, a one-sided formula (~col), or a single character string.", call. = FALSE)
}

#' Routing table: returns list(fn_power, fn_ss, model, method_used)
#' @noRd
.rmst_route <- function(dep_cens, has_smooth, strata_var, strata_type, type) {
  method_used <- type
  if (dep_cens) {
    if (type == "boot") {
      message("DC model only supports analytical method; switching type to 'analytical'.")
      method_used <- "analytical"
    }
    return(list(fn_power    = DC.power.analytical,
                fn_ss       = DC.ss.analytical,
                model       = "DC",
                method_used = method_used))
  }
  if (has_smooth) {
    if (type == "analytical") {
      message("GAM model requires bootstrap; switching type to 'boot'.")
      method_used <- "boot"
    }
    return(list(fn_power    = GAM.power.boot,
                fn_ss       = GAM.ss.boot,
                model       = "GAM",
                method_used = method_used))
  }
  if (is.null(strata_var)) {
    if (type == "analytical")
      return(list(fn_power    = linear.power.analytical,
                  fn_ss       = linear.ss.analytical,
                  model       = "linear",
                  method_used = method_used))
    else
      return(list(fn_power    = linear.power.boot,
                  fn_ss       = linear.ss.boot,
                  model       = "linear",
                  method_used = method_used))
  }
  st <- match.arg(strata_type, c("additive", "multiplicative"))
  if (st == "additive") {
    if (type == "analytical")
      return(list(fn_power    = additive.power.analytical,
                  fn_ss       = additive.ss.analytical,
                  model       = "additive",
                  method_used = method_used))
    else
      return(list(fn_power    = GAM.power.boot,
                  fn_ss       = GAM.ss.boot,
                  model       = "GAM",
                  method_used = method_used))
  }
  # multiplicative
  if (type == "analytical")
    return(list(fn_power    = MS.power.analytical,
                fn_ss       = MS.ss.analytical,
                model       = "multiplicative",
                method_used = method_used))
  else
    return(list(fn_power    = MS.power.boot,
                fn_ss       = MS.ss.boot,
                model       = "multiplicative",
                method_used = method_used))
}

#' Build argument list for the underlying power function call
#' @noRd
.build_args_power <- function(parsed, arm, sample_sizes, L, strata_var, alpha,
                               n_sim, parallel.cores, route, data,
                               strata_type, dep_cens) {
  base <- list(
    pilot_data   = data,
    time_var     = parsed$time_var,
    status_var   = parsed$status_var,
    arm_var      = arm,
    sample_sizes = sample_sizes,
    L            = L,
    alpha        = alpha
  )
  lt <- parsed$linear_terms
  if (!is.null(lt)) base$linear_terms <- lt

  m  <- route$model
  mt <- route$method_used

  if (m %in% c("additive", "multiplicative", "GAM")) {
    base$strata_var <- strata_var
  }
  if (m == "GAM") {
    if (!is.null(parsed$smooth_terms)) base$smooth_terms <- parsed$smooth_terms
    base$n_sim          <- n_sim
    base$parallel.cores <- parallel.cores
  } else if (m %in% c("additive", "multiplicative") && mt == "boot") {
    base$n_sim          <- n_sim
    base$parallel.cores <- parallel.cores
  } else if (m == "linear" && mt == "boot") {
    base$n_sim <- n_sim
    # linear.power.boot has no parallel.cores parameter
  }
  # DC analytical has no n_sim parameter
  base
}

#' Build argument list for the underlying SS function call
#' @noRd
.build_args_ss <- function(parsed, arm, target_power, L, strata_var, alpha,
                            n_sim, parallel.cores, route, data,
                            n_start, n_step, max_n, patience) {
  base <- list(
    pilot_data   = data,
    time_var     = parsed$time_var,
    status_var   = parsed$status_var,
    arm_var      = arm,
    target_power = target_power,
    L            = L,
    alpha        = alpha
  )
  lt <- parsed$linear_terms
  if (!is.null(lt)) base$linear_terms <- lt

  m  <- route$model
  mt <- route$method_used

  if (m %in% c("additive", "multiplicative", "GAM")) {
    base$strata_var <- strata_var
  }
  if (m == "GAM") {
    if (!is.null(parsed$smooth_terms)) base$smooth_terms <- parsed$smooth_terms
    base$n_sim          <- n_sim
    base$parallel.cores <- parallel.cores
  } else if (m %in% c("additive", "multiplicative") && mt == "boot") {
    base$n_sim          <- n_sim
    base$parallel.cores <- parallel.cores
  } else if (m == "linear" && mt == "boot") {
    base$n_sim <- n_sim
    # linear.ss.boot has no parallel.cores parameter
  }
  # DC analytical has no n_sim; all analytical SS functions have no patience
  # All SS functions use max_n_per_arm (not max_n)
  base$n_start       <- n_start
  base$n_step        <- n_step
  base$max_n_per_arm <- max_n
  if (mt == "boot") base$patience <- patience
  base
}

# â”€â”€ rmst.power â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

#' Power analysis for RMST-based models via formula interface
#'
#' @description
#' Routes a formula-based call to the matching RMST power routine.
#'
#' @param formula A formula of the form \code{Surv(time, status) ~ cov1 + cov2}.
#'   Use \code{s()} wrapping (mgcv style) for smooth terms: \code{Surv(time, status) ~ s(age)}.
#' @param data A \code{data.frame} containing the reference (pilot) data.
#' @param arm Character string naming the treatment arm column (binary 0/1).
#' @param sample_sizes Integer vector of per-arm (or per-stratum) sample sizes to evaluate.
#' @param L Numeric truncation time for RMST.
#' @param strata Character column name, one-sided formula (\code{~col}), or \code{NULL}
#'   (default). Ignored when \code{dep_cens = TRUE}.
#' @param strata_type One of \code{"additive"} (default) or \code{"multiplicative"}.
#'   Only used when \code{strata} is non-\code{NULL} and \code{dep_cens = FALSE}.
#' @param dep_cens Logical; use dependent-censoring model? Default \code{FALSE}.
#' @param type One of \code{"analytical"} or \code{"boot"}. Auto-switched with a message
#'   when the requested type is unavailable for the chosen model.
#' @param alpha Significance level. Default \code{0.05}.
#' @param n_sim Number of bootstrap replicates (boot methods only). Default \code{1000}.
#' @param parallel.cores Number of cores for parallel processing. Default \code{1}.
#'
#' @return An object of class \code{c("rmst_power", "list")} with elements
#'   \code{results_data}, \code{results_plot}, \code{results_summary},
#'   \code{model_output}, and \code{.meta}.
#'
#' @seealso \code{\link{rmst.ss}}, \code{\link{print.rmst_power}},
#'   \code{\link{summary.rmst_power}}, \code{\link{plot.rmst_power}}
#'
#' @examples
#' r <- rmst.power(Surv(time, status) ~ age,
#'                 data = aft_lognormal_L12_n150,
#'                 arm  = "arm",
#'                 sample_sizes = c(50, 100, 150),
#'                 L    = 12)
#' print(r)
#' s <- summary(r)
#' plot(r)
#'
#' @export
rmst.power <- function(formula,
                       data,
                       arm,
                       sample_sizes,
                       L,
                       strata      = NULL,
                       strata_type = c("additive", "multiplicative"),
                       dep_cens    = FALSE,
                       type        = c("analytical", "boot"),
                       alpha       = 0.05,
                       n_sim       = 1000L,
                       parallel.cores = 1L) {
  mc         <- match.call()
  type       <- match.arg(type)
  strata_type <- match.arg(strata_type)
  parsed     <- .parse_rmst_formula(formula)
  strata_var <- .resolve_strata(strata)
  has_smooth <- !is.null(parsed$smooth_terms)
  route      <- .rmst_route(dep_cens, has_smooth, strata_var, strata_type, type)

  args <- .build_args_power(parsed, arm, sample_sizes, L, strata_var, alpha,
                             n_sim, parallel.cores, route, data, strata_type, dep_cens)
  raw  <- do.call(route$fn_power, args)

  n_col <- if (route$model %in% c("multiplicative", "GAM") && !is.null(strata_var))
             "N_per_Stratum" else "N_per_Arm"

  structure(
    list(
      results_data    = raw$results_data,
      results_plot    = raw$results_plot,
      results_summary = raw$results_summary,
      model_output    = raw$model_output,
      .meta = list(
        model          = route$model,
        method         = route$method_used,
        type_requested = type,
        formula        = formula,
        arm            = arm,
        strata         = strata_var,
        strata_type    = strata_type,
        dep_cens       = dep_cens,
        truncation     = L,
        alpha          = alpha,
        n_col          = n_col,
        call           = mc
      )
    ),
    class = c("rmst_power", "list")
  )
}

# â”€â”€ rmst.ss â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

#' Sample size estimation for RMST-based models via formula interface
#'
#' @description
#' Routes a formula-based call to the matching RMST sample-size routine.
#'
#' @param formula A formula of the form \code{Surv(time, status) ~ cov1 + cov2}.
#' @param data A \code{data.frame} containing the reference (pilot) data.
#' @param arm Character string naming the treatment arm column (binary 0/1).
#' @param target_power Numeric target power (e.g., \code{0.80}).
#' @param L Numeric truncation time for RMST.
#' @param strata Character column name, one-sided formula (\code{~col}), or \code{NULL}.
#' @param strata_type One of \code{"additive"} (default) or \code{"multiplicative"}.
#' @param dep_cens Logical; use dependent-censoring model? Default \code{FALSE}.
#' @param type One of \code{"analytical"} or \code{"boot"}.
#' @param alpha Significance level. Default \code{0.05}.
#' @param n_sim Number of bootstrap replicates. Default \code{1000}.
#' @param parallel.cores Number of cores for parallel processing. Default \code{1}.
#' @param n_start Starting sample size for the search. Default \code{50}.
#' @param n_step Search increment. Default \code{25}.
#' @param max_n Maximum sample size to try. Default \code{2000}.
#' @param patience Number of consecutive non-improving steps before stopping. Default \code{5}.
#'
#' @return An object of class \code{c("rmst_ss", "list")} with elements
#'   \code{results_data}, \code{results_plot}, \code{results_summary},
#'   \code{model_output}, and \code{.meta}.
#'
#' @seealso \code{\link{rmst.power}}, \code{\link{print.rmst_ss}},
#'   \code{\link{summary.rmst_ss}}, \code{\link{plot.rmst_ss}}
#'
#' @examples
#' r <- rmst.ss(Surv(time, status) ~ age,
#'              data         = aft_lognormal_L12_n150,
#'              arm          = "arm",
#'              target_power = 0.80,
#'              L            = 12)
#' print(r)
#'
#' @export
rmst.ss <- function(formula,
                    data,
                    arm,
                    target_power,
                    L,
                    strata      = NULL,
                    strata_type = c("additive", "multiplicative"),
                    dep_cens    = FALSE,
                    type        = c("analytical", "boot"),
                    alpha       = 0.05,
                    n_sim       = 1000L,
                    parallel.cores = 1L,
                    n_start     = 50L,
                    n_step      = 25L,
                    max_n       = 2000L,
                    patience    = 5L) {
  mc          <- match.call()
  type        <- match.arg(type)
  strata_type <- match.arg(strata_type)
  parsed      <- .parse_rmst_formula(formula)
  strata_var  <- .resolve_strata(strata)
  has_smooth  <- !is.null(parsed$smooth_terms)
  route       <- .rmst_route(dep_cens, has_smooth, strata_var, strata_type, type)

  args <- .build_args_ss(parsed, arm, target_power, L, strata_var, alpha,
                          n_sim, parallel.cores, route, data,
                          n_start, n_step, max_n, patience)
  raw  <- do.call(route$fn_ss, args)

  n_col <- if (route$model %in% c("multiplicative", "GAM") && !is.null(strata_var))
             "N_per_Stratum" else "N_per_Arm"

  structure(
    list(
      results_data    = raw$results_data,
      results_plot    = raw$results_plot,
      results_summary = raw$results_summary,
      model_output    = raw$model_output,
      .meta = list(
        model          = route$model,
        method         = route$method_used,
        type_requested = type,
        formula        = formula,
        arm            = arm,
        strata         = strata_var,
        strata_type    = strata_type,
        dep_cens       = dep_cens,
        truncation     = L,
        alpha          = alpha,
        n_col          = n_col,
        call           = mc
      )
    ),
    class = c("rmst_ss", "list")
  )
}

# â”€â”€ S3 methods â€” print â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

.rmst_model_label <- function(meta) {
  m  <- meta$model
  mt <- if (m == "linear")         "Linear IPCW"
        else if (m == "additive")  "Additive Stratified"
        else if (m == "multiplicative") "Multiplicative Stratified"
        else if (m == "GAM")       "GAM (Additive)"
        else if (m == "DC")        "Dependent Censoring"
        else m
  md <- if (meta$method == "analytical") "Analytical" else "Bootstrap"
  paste0(mt, " (", md, ")")
}

#' Print an rmst_power result
#'
#' Prints the model metadata and power table returned by [rmst.power()].
#'
#' @param x An object returned by [rmst.power()].
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return The input object `x`, invisibly.
#'
#' @rdname rmst_power_methods
#' @export
print.rmst_power <- function(x, ...) {
  meta  <- x$.meta
  label <- .rmst_model_label(meta)
  cat("\u2500\u2500 RMST Power Analysis ", strrep("\u2500", max(0L, 50L - nchar(label))), "\n", sep = "")
  cat("  Model           :", label, "\n")
  cat("  Formula         :", deparse(meta$formula), "\n")
  cat("  Arm             :", meta$arm, "\n")
  if (!is.null(meta$strata))
    cat("  Strata          :", meta$strata, "\n")
  cat("  Truncation time :", meta$truncation, "\n")
  cat("  Alpha           :", meta$alpha, "\n")
  if (!identical(meta$type_requested, meta$method))
    cat("  Note            : type switched from '", meta$type_requested,
        "' to '", meta$method, "'\n", sep = "")
  cat("\nPower Results:\n")
  rd <- x$results_data
  print(rd, row.names = FALSE)
  te <- x$model_output$treatment_effect
  if (!is.null(te) && nrow(te) >= 1L) {
    cat("\nTreatment Effect (from reference data):\n")
    cat("  Estimand :", te$estimand[1L], "\n")
    cat("  Estimate :", round(te$estimate[1L], 4L), "\n")
  }
  invisible(x)
}

#' Print an rmst_ss result
#'
#' Prints the model metadata and sample-size table returned by [rmst.ss()].
#'
#' @param x An object returned by [rmst.ss()].
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return The input object `x`, invisibly.
#'
#' @rdname rmst_ss_methods
#' @export
print.rmst_ss <- function(x, ...) {
  meta  <- x$.meta
  label <- .rmst_model_label(meta)
  cat("\u2500\u2500 RMST Sample Size Estimation ", strrep("\u2500", max(0L, 40L - nchar(label))), "\n", sep = "")
  cat("  Model           :", label, "\n")
  cat("  Formula         :", deparse(meta$formula), "\n")
  cat("  Arm             :", meta$arm, "\n")
  if (!is.null(meta$strata))
    cat("  Strata          :", meta$strata, "\n")
  cat("  Truncation time :", meta$truncation, "\n")
  cat("  Alpha           :", meta$alpha, "\n")
  if (!identical(meta$type_requested, meta$method))
    cat("  Note            : type switched from '", meta$type_requested,
        "' to '", meta$method, "'\n", sep = "")
  cat("\nSample Size Result:\n")
  print(x$results_data, row.names = FALSE)
  te <- x$model_output$treatment_effect
  if (!is.null(te) && nrow(te) >= 1L) {
    cat("\nTreatment Effect (from reference data):\n")
    cat("  Estimand :", te$estimand[1L], "\n")
    cat("  Estimate :", round(te$estimate[1L], 4L), "\n")
  }
  invisible(x)
}

# â”€â”€ S3 methods â€” summary â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

#' Summarize an rmst_power result
#'
#' Returns a summary object for printing and further inspection.
#'
#' @param object An object returned by [rmst.power()].
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return An object of class `"summary.rmst_power"`.
#'
#' @rdname rmst_power_methods
#' @export
summary.rmst_power <- function(object, ...) {
  meta <- object$.meta
  structure(
    list(
      model_info = data.frame(
        Model           = meta$model,
        Method          = meta$method,
        Truncation_time = meta$truncation,
        Alpha           = meta$alpha,
        Arm_variable    = meta$arm,
        Strata          = if (!is.null(meta$strata)) meta$strata else "None",
        Dep_censoring   = meta$dep_cens,
        stringsAsFactors = FALSE
      ),
      power_results       = object$results_data,
      effect_summary      = object$results_summary,
      coefficient_table   = object$model_output$coefficient_table,
      treatment_effect    = object$model_output$treatment_effect,
      arm_specific_rmst   = object$model_output$arm_specific_rmst,
      variance_components = object$model_output$variance_components,
      censoring_weights   = object$model_output$censoring_weights,
      diagnostics         = object$model_output$diagnostics,
      simulation_draws    = object$model_output$simulation_draws
    ),
    class = "summary.rmst_power"
  )
}

#' Summarize an rmst_ss result
#'
#' Returns a summary object for printing and further inspection.
#'
#' @param object An object returned by [rmst.ss()].
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return An object of class `"summary.rmst_ss"`.
#'
#' @rdname rmst_ss_methods
#' @export
summary.rmst_ss <- function(object, ...) {
  meta <- object$.meta
  structure(
    list(
      model_info = data.frame(
        Model           = meta$model,
        Method          = meta$method,
        Truncation_time = meta$truncation,
        Alpha           = meta$alpha,
        Arm_variable    = meta$arm,
        Strata          = if (!is.null(meta$strata)) meta$strata else "None",
        Dep_censoring   = meta$dep_cens,
        stringsAsFactors = FALSE
      ),
      ss_results          = object$results_data,
      effect_summary      = object$results_summary,
      coefficient_table   = object$model_output$coefficient_table,
      treatment_effect    = object$model_output$treatment_effect,
      arm_specific_rmst   = object$model_output$arm_specific_rmst,
      variance_components = object$model_output$variance_components,
      censoring_weights   = object$model_output$censoring_weights,
      diagnostics         = object$model_output$diagnostics,
      simulation_draws    = object$model_output$simulation_draws
    ),
    class = "summary.rmst_ss"
  )
}

# â”€â”€ S3 methods â€” print.summary â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

.print_summary_rmst <- function(x, result_field, result_label) {
  cat("\u2500\u2500 Model Information ", strrep("\u2500", 42L), "\n", sep = "")
  print(x$model_info, row.names = FALSE)

  cat("\n\u2500\u2500", result_label, strrep("\u2500", max(0L, 55L - nchar(result_label))), "\n")
  print(x[[result_field]], row.names = FALSE)

  if (!is.null(x$effect_summary)) {
    cat("\n\u2500\u2500 Effect Summary ", strrep("\u2500", 44L), "\n", sep = "")
    print(x$effect_summary, row.names = FALSE)
  }

  if (!is.null(x$treatment_effect)) {
    cat("\n\u2500\u2500 Treatment Effect ", strrep("\u2500", 42L), "\n", sep = "")
    print(x$treatment_effect, row.names = FALSE)
  }

  if (!is.null(x$coefficient_table)) {
    cat("\n\u2500\u2500 Coefficient Table ", strrep("\u2500", 41L), "\n", sep = "")
    print(x$coefficient_table, row.names = FALSE)
  } else {
    cat("\n  coefficient_table: not available (bootstrap method)\n")
  }

  if (!is.null(x$arm_specific_rmst)) {
    cat("\n\u2500\u2500 Arm-Specific RMST ", strrep("\u2500", 40L), "\n", sep = "")
    print(x$arm_specific_rmst, row.names = FALSE)
    if (any(is.na(x$arm_specific_rmst$std_error)))
      cat("  Note: SE/CI not available for arm-specific RMST in analytical methods.\n")
  }

  if (!is.null(x$variance_components)) {
    cat("\n\u2500\u2500 Variance Components ", strrep("\u2500", 39L), "\n", sep = "")
    vc <- x$variance_components
    if (!is.null(vc$se_effect_n1))
      cat("  se_effect (n=1):", round(vc$se_effect_n1, 6L), "\n")
    if (!is.null(vc$A_hat))
      { cat("  A_hat :\n"); print(round(vc$A_hat, 6L)) }
    if (!is.null(vc$B_hat))
      { cat("  B_hat :\n"); print(round(vc$B_hat, 6L)) }
    if (!is.null(vc$V_hat_n))
      { cat("  V_hat_n :\n"); print(round(vc$V_hat_n, 6L)) }
  }

  if (!is.null(x$censoring_weights)) {
    cat("\n\u2500\u2500 Censoring Weights ", strrep("\u2500", 41L), "\n", sep = "")
    cw <- x$censoring_weights
    if (!is.null(cw$raw_summary)) {
      cat("  Distribution:\n")
      print(round(cw$raw_summary, 4L))
    }
    if (!is.null(cw$cap_value))
      cat("  Cap value       :", round(cw$cap_value, 4L), "\n")
    if (!is.null(cw$capped_fraction) && !is.na(cw$capped_fraction))
      cat("  Capped fraction :", round(cw$capped_fraction, 4L), "\n")
  }

  if (!is.null(x$diagnostics)) {
    cat("\n\u2500\u2500 Diagnostics ", strrep("\u2500", 47L), "\n", sep = "")
    dg <- x$diagnostics
    if (!is.null(dg$n_used))    cat("  n used          :", dg$n_used, "\n")
    if (!is.null(dg$n_events))  cat("  n events        :", dg$n_events, "\n")
    if (!is.null(dg$n_sim))     cat("  n simulations   :", dg$n_sim, "\n")
    if (!is.null(dg$convergence_ok)) cat("  convergence_ok  :", dg$convergence_ok, "\n")
    if (!is.null(dg$singular_flag))  cat("  singular_flag   :", dg$singular_flag, "\n")
  }

  if (!is.null(x$simulation_draws)) {
    cat("\n\u2500\u2500 Simulation Draws ", strrep("\u2500", 42L), "\n", sep = "")
    cat("  (first 6 rows; access via $simulation_draws for full data)\n")
    print(utils::head(x$simulation_draws, 6L), row.names = FALSE)
  }

  invisible(x)
}

#' @export
print.summary.rmst_power <- function(x, ...) {
  .print_summary_rmst(x, "power_results", "Power Results")
}

#' @export
print.summary.rmst_ss <- function(x, ...) {
  .print_summary_rmst(x, "ss_results", "Sample Size Result")
}

# â”€â”€ S3 methods â€” plot â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

#' Plot an rmst_power result
#'
#' Prints and returns the stored `ggplot2` object.
#'
#' @param x An object returned by [rmst.power()].
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return A `ggplot2` object, invisibly.
#'
#' @rdname rmst_power_methods
#' @export
plot.rmst_power <- function(x, ...) {
  print(x$results_plot)
  invisible(x$results_plot)
}

#' Plot an rmst_ss result
#'
#' Prints and returns the stored `ggplot2` object.
#'
#' @param x An object returned by [rmst.ss()].
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return A `ggplot2` object, invisibly.
#'
#' @rdname rmst_ss_methods
#' @export
plot.rmst_ss <- function(x, ...) {
  print(x$results_plot)
  invisible(x$results_plot)
}

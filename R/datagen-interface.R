# datagen-interface.R â€” Convenient data generation wrappers
# Provides covariate helpers, recipe_quick_ph(), rmst.sim(), and S3 methods.

# â”€â”€ Internal helper â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

.rmst_sim_model_label <- function(model) {
  switch(model,
    aft_lognormal  = "AFT Log-Normal",
    aft_weibull    = "AFT Weibull",
    aft_loglogistic = "AFT Log-Logistic",
    ph_exponential = "PH Exponential",
    ph_weibull     = "PH Weibull",
    ph_gompertz    = "PH Gompertz",
    cox_pwexp      = "PH Piecewise-Exponential (Cox)",
    model
  )
}

# â”€â”€ Covariate helpers â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

#' Continuous covariate definition
#'
#' @description
#' Builds a covariate definition list for use in \code{\link{rmst.sim}},
#' \code{\link{recipe_quick_aft}}, or \code{\link{recipe_quick_ph}}.
#'
#' @param name Column name in the generated data.
#' @param dist Distribution name. One of \code{"normal"} (default),
#'   \code{"lognormal"}, \code{"gamma"}, \code{"weibull"}, \code{"uniform"},
#'   \code{"beta"}, \code{"t"}.
#' @param ... Named distribution parameters passed as the \code{params} list
#'   (e.g., \code{mean = 0, sd = 1} for \code{dist = "normal"}).
#'
#' @return A named list suitable as one element of the \code{covariates} argument.
#'
#' @examples
#' covar_continuous("age", dist = "normal", mean = 50, sd = 10)
#' covar_continuous("bmi", dist = "lognormal", meanlog = 3.2, sdlog = 0.2)
#'
#' @seealso \code{\link{covar_binary}}, \code{\link{covar_categorical}},
#'   \code{\link{rmst.sim}}
#' @export
covar_continuous <- function(name, dist = "normal", ...) {
  dist <- match.arg(dist, c("normal","lognormal","gamma","weibull","uniform","beta","t"))
  list(name = name, type = "continuous", dist = dist, params = list(...))
}

#' Binary covariate definition
#'
#' @description
#' Builds a binary (Bernoulli) covariate definition list.
#'
#' @param name Column name in the generated data.
#' @param p Probability of value 1. Default \code{0.5}.
#'
#' @return A named list suitable as one element of the \code{covariates} argument.
#'
#' @examples
#' covar_binary("female", p = 0.52)
#'
#' @seealso \code{\link{covar_continuous}}, \code{\link{covar_categorical}},
#'   \code{\link{rmst.sim}}
#' @export
covar_binary <- function(name, p = 0.5) {
  if (!is.numeric(p) || length(p) != 1L || p < 0 || p > 1)
    stop("'p' must be a single numeric value in [0, 1].", call. = FALSE)
  list(name = name, type = "categorical", dist = "bernoulli", params = list(p = p))
}

#' Categorical covariate definition
#'
#' @description
#' Builds a multi-level categorical covariate definition list.
#'
#' @param name Column name in the generated data.
#' @param probs Numeric vector of category probabilities (must sum to 1).
#' @param labels Optional character vector of level labels. Length must equal
#'   \code{length(probs)}. Defaults to \code{"1"}, \code{"2"}, ...
#'
#' @return A named list suitable as one element of the \code{covariates} argument.
#'
#' @examples
#' covar_categorical("stage", probs = c(0.4, 0.35, 0.25),
#'                   labels = c("I", "II", "III"))
#'
#' @seealso \code{\link{covar_continuous}}, \code{\link{covar_binary}},
#'   \code{\link{rmst.sim}}
#' @export
covar_categorical <- function(name, probs, labels = NULL) {
  if (!is.numeric(probs) || any(probs < 0))
    stop("'probs' must be a non-negative numeric vector.", call. = FALSE)
  if (abs(sum(probs) - 1) > 1e-6)
    stop("'probs' must sum to 1.", call. = FALSE)
  if (!is.null(labels) && length(labels) != length(probs))
    stop("'labels' must have the same length as 'probs'.", call. = FALSE)
  labs <- labels %||% as.character(seq_along(probs))
  list(name = name, type = "categorical", dist = "categorical",
       params = list(prob = probs, labels = labs))
}

# â”€â”€ recipe_quick_ph â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

#' Quick PH recipe builder
#'
#' @description
#' Convenience wrapper that builds a recipe list for proportional-hazard (PH)
#' models, parallel to \code{\link{recipe_quick_aft}} for AFT models.
#' The returned recipe is validated and ready for \code{\link{simulate_from_recipe}}.
#'
#' @param n Total sample size (integer).
#' @param model PH model. One of \code{"ph_exponential"}, \code{"ph_weibull"},
#'   \code{"ph_gompertz"}, \code{"cox_pwexp"}.
#' @param baseline Named list of baseline hazard parameters:
#' \describe{
#'   \item{ph_exponential}{\code{list(rate = ...)}}
#'   \item{ph_weibull}{\code{list(shape = ..., scale = ...)}}
#'   \item{ph_gompertz}{\code{list(shape = ..., rate = ...)}}
#'   \item{cox_pwexp}{\code{list(rates = c(...), cuts = c(...))}}
#' }
#' @param treat_effect Numeric log-hazard ratio for the treatment arm.
#' @param covariates List of covariate definitions. Use \code{\link{covar_continuous}},
#'   \code{\link{covar_binary}}, or \code{\link{covar_categorical}}. Default \code{list()}.
#' @param target_censoring Target overall censoring fraction (0â€“1). Default \code{0.25}.
#' @param allocation Treatment allocation ratio string, e.g. \code{"1:1"} (default)
#'   or \code{"1:2"}.
#' @param seed Optional integer seed.
#'
#' @return A recipe list suitable for \code{\link{simulate_from_recipe}}.
#'
#' @examples
#' r <- recipe_quick_ph(100, "ph_weibull",
#'        baseline = list(shape = 1.5, scale = 10),
#'        treat_effect = -0.5,
#'        covariates = list(covar_continuous("age")),
#'        target_censoring = 0.30)
#' dat <- simulate_from_recipe(r, seed = 1)
#'
#' @seealso \code{\link{recipe_quick_aft}}, \code{\link{rmst.sim}},
#'   \code{\link{simulate_from_recipe}}
#' @export
recipe_quick_ph <- function(n,
                             model = c("ph_exponential", "ph_weibull",
                                       "ph_gompertz", "cox_pwexp"),
                             baseline,
                             treat_effect,
                             covariates       = list(),
                             target_censoring = 0.25,
                             allocation       = "1:1",
                             seed             = NULL) {
  model <- match.arg(model)
  validate_recipe(list(
    n          = as.integer(n),
    covariates = list(defs = covariates),
    treatment  = list(assignment = "randomization", allocation = allocation),
    event_time = list(
      model    = model,
      baseline = baseline,
      effects  = list(intercept = 0, treatment = treat_effect, covariates = NULL)
    ),
    censoring  = list(mode = "target_overall", target = target_censoring,
                      admin_time = Inf),
    seed       = seed
  ))
}

# â”€â”€ rmst.sim â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

#' Simulate survival data for RMST analysis
#'
#' @description
#' A unified, single-call wrapper for generating survival data suitable for use
#' as reference/pilot data in \code{\link{rmst.power}} and \code{\link{rmst.ss}}.
#' Supports all seven built-in event-time models (AFT and PH families).
#'
#' @details
#' Internally routes to \code{\link{recipe_quick_aft}} for AFT models and
#' \code{\link{recipe_quick_ph}} for PH models, then calls
#' \code{\link{simulate_from_recipe}}.
#'
#' @param n Total sample size (split by \code{allocation} ratio).
#' @param model Event-time model. One of:
#'   \code{"aft_lognormal"} (default), \code{"aft_weibull"},
#'   \code{"aft_loglogistic"}, \code{"ph_exponential"}, \code{"ph_weibull"},
#'   \code{"ph_gompertz"}, \code{"cox_pwexp"}.
#' @param baseline Named list of baseline parameters (model-specific; see
#'   \code{\link{recipe_quick_aft}} and \code{\link{recipe_quick_ph}}).
#' @param treat_effect Numeric treatment coefficient. Log-time scale for AFT
#'   models; log-hazard ratio for PH models. Default \code{0}.
#' @param covariates List of covariate definitions. Elements can be created with
#'   \code{\link{covar_continuous}}, \code{\link{covar_binary}}, or
#'   \code{\link{covar_categorical}}, or supplied as raw named lists in the
#'   existing recipe format. Default \code{list()} (no covariates).
#' @param target_censoring Target overall censoring fraction (0â€“1). Default \code{0.25}.
#' @param allocation Treatment allocation ratio string, e.g. \code{"1:1"} (default).
#' @param L Optional numeric truncation time. Stored as an attribute on the
#'   returned object for downstream use by \code{\link{rmst.power}} /
#'   \code{\link{rmst.ss}}. Does not affect data generation.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A \code{data.frame} of class \code{c("rmst_simdata", "data.frame")}
#'   with columns \code{time}, \code{status}, \code{arm} (when treatment is
#'   present), and one column per covariate. Attributes:
#' \describe{
#'   \item{recipe}{The validated recipe list used for generation.}
#'   \item{L}{The truncation time if supplied, else \code{NULL}.}
#'   \item{achieved_censoring}{Actual censoring fraction achieved.}
#' }
#'
#' @examples
#' df <- rmst.sim(
#'   n            = 150,
#'   model        = "aft_lognormal",
#'   baseline     = list(mu = 2.2, sigma = 0.5),
#'   treat_effect = -0.3,
#'   covariates   = list(covar_continuous("age"), covar_binary("female")),
#'   L            = 12,
#'   seed         = 42
#' )
#' print(df)
#' s <- summary(df)
#'
#' @seealso \code{\link{rmst.power}}, \code{\link{rmst.ss}},
#'   \code{\link{recipe_quick_ph}}, \code{\link{recipe_quick_aft}},
#'   \code{\link{covar_continuous}}
#' @export
rmst.sim <- function(n,
                     model            = "aft_lognormal",
                     baseline,
                     treat_effect     = 0,
                     covariates       = list(),
                     target_censoring = 0.25,
                     allocation       = "1:1",
                     L                = NULL,
                     seed             = NULL) {
  aft_quick  <- c("aft_lognormal", "aft_weibull")
  aft_models <- c(aft_quick, "aft_loglogistic")
  ph_models  <- c("ph_exponential", "ph_weibull", "ph_gompertz", "cox_pwexp")
  model      <- match.arg(model, c(aft_models, ph_models))

  if (model %in% aft_quick) {
    recipe <- recipe_quick_aft(n, model = model, baseline = baseline,
                               treat_effect = treat_effect,
                               covariates = covariates,
                               target_censoring = target_censoring,
                               allocation = allocation, seed = seed)
  } else if (model == "aft_loglogistic") {
    # recipe_quick_aft only accepts lognormal/weibull; build directly for loglogistic
    recipe <- validate_recipe(list(
      n          = as.integer(n),
      covariates = list(defs = covariates),
      treatment  = list(assignment = "randomization", allocation = allocation),
      event_time = list(
        model    = "aft_loglogistic",
        baseline = baseline,
        effects  = list(intercept = 0, treatment = treat_effect, covariates = NULL)
      ),
      censoring  = list(mode = "target_overall", target = target_censoring,
                        admin_time = Inf),
      seed       = seed
    ))
  } else {
    recipe <- recipe_quick_ph(n, model = model, baseline = baseline,
                              treat_effect = treat_effect,
                              covariates = covariates,
                              target_censoring = target_censoring,
                              allocation = allocation, seed = seed)
  }

  dat <- simulate_from_recipe(recipe, seed = seed)

  out <- structure(
    dat,
    class             = c("rmst_simdata", "data.frame"),
    recipe            = recipe,
    L                 = L,
    achieved_censoring = attr(dat, "achieved_censoring", exact = TRUE)
  )
  out
}

# â”€â”€ S3 methods â€” print â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

#' @export
print.rmst_simdata <- function(x, ...) {
  recipe <- attr(x, "recipe", exact = TRUE)
  model  <- recipe$event_time$model
  label  <- .rmst_sim_model_label(model)
  n_tot  <- nrow(x)
  n_ev   <- sum(x$status == 1L)
  n_cens <- n_tot - n_ev
  L_val  <- attr(x, "L", exact = TRUE)

  cat("\u2500\u2500 Simulated RMST Dataset ", strrep("\u2500", 36L), "\n", sep = "")
  cat("  Model          :", label, "\n")

  if ("arm" %in% names(x)) {
    n1 <- sum(x$arm == 1L)
    n0 <- n_tot - n1
    cat(sprintf("  N (total)      : %d  |  Arm 0: %d  |  Arm 1: %d\n", n_tot, n0, n1))
  } else {
    cat("  N (total)      :", n_tot, "\n")
  }

  cat(sprintf("  Events         : %d (%.1f%%)\n", n_ev,   100 * n_ev   / n_tot))
  cat(sprintf("  Censored       : %d (%.1f%%)\n", n_cens, 100 * n_cens / n_tot))

  if (!is.null(L_val))
    cat("  Truncation time:", L_val, " (stored as attribute)\n")

  cat("\nFirst 6 rows:\n")
  print(utils::head(as.data.frame(x), 6L))
  invisible(x)
}

# â”€â”€ S3 methods â€” summary â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

#' @export
summary.rmst_simdata <- function(object, ...) {
  recipe <- attr(object, "recipe", exact = TRUE)

  meta <- list(
    dataset_id          = "sim_001",
    scenario_id         = 1L,
    rep                 = 1L,
    seed_used           = recipe$seed %||% NA_integer_,
    n                   = nrow(object),
    n_treat             = if ("arm" %in% names(object)) sum(object$arm == 1L) else NA_integer_,
    n_control           = if ("arm" %in% names(object)) sum(object$arm == 0L) else NA_integer_,
    event_rate          = mean(object$status == 1L),
    achieved_censoring  = attr(object, "achieved_censoring", exact = TRUE) %||%
                            mean(object$status == 0L),
    model               = recipe$event_time$model,
    baseline            = recipe$event_time$baseline,
    effects             = recipe$event_time$effects,
    treatment           = recipe$treatment %||% NULL,
    censoring           = recipe$censoring %||% NULL,
    covariates          = if (!is.null(recipe$covariates$defs))
                            lapply(recipe$covariates$defs,
                                   function(d) d[c("name","type","dist","params")]) else NULL,
    allocation          = (recipe$treatment %||% list())$allocation %||% NA_character_,
    params              = list(),
    files               = list(),
    created_at          = as.character(Sys.time())
  )

  set <- list(data = as.data.frame(object), meta = meta)
  desc <- describe_generation(set)

  structure(
    list(
      header     = desc$header,
      effects    = desc$effects,
      treatment  = desc$treatment,
      censoring  = desc$censoring,
      covariates = desc$covariates,
      baseline   = desc$baseline
    ),
    class = "summary.rmst_simdata"
  )
}

#' @export
print.summary.rmst_simdata <- function(x, ...) {
  .print_section <- function(label, df) {
    if (is.null(df)) return(invisible(NULL))
    cat("\n\u2500\u2500", label, strrep("\u2500", max(0L, 55L - nchar(label))), "\n")
    print(df, row.names = FALSE)
  }
  cat("\u2500\u2500 Simulated Dataset Summary ", strrep("\u2500", 33L), "\n", sep = "")
  if (!is.null(x$header)) print(x$header, row.names = FALSE)
  .print_section("Baseline Parameters", x$baseline)
  .print_section("Treatment Effects",   x$effects)
  .print_section("Treatment Assignment", x$treatment)
  .print_section("Censoring",           x$censoring)
  .print_section("Covariates",          x$covariates)
  invisible(x)
}

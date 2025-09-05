# R/describe_generation.R

#' Summarize a generated dataset and its simulation mechanism
#'
#' Given one element returned by [load_recipe_sets()] (a list with \code{data} and
#' a rich \code{meta} list), this builds tidy tables describing the data-generating
#' process: model family, baseline parameters, linear predictor coefficients
#' (intercept / treatment / covariates), treatment assignment, and censoring.
#'
#' @param set A single element from [load_recipe_sets()], i.e. \code{list(data=<df>, meta=<list>)}.
#' @return A list of data frames:
#' \itemize{
#'   \item \code{header}: \code{n}, \code{L} (analysis horizon; stored as attr \code{"tau"}), \code{model},
#'         \code{event_rate}, \code{achieved_censoring}.
#'   \item \code{baseline}: flattened baseline parameters.
#'   \item \code{effects}: coefficients for intercept/treatment and covariates (and formula terms if used).
#'   \item \code{treatment}: assignment mechanism and knobs.
#'   \item \code{censoring}: censoring mode/target/admin time.
#'   \item \code{covariates}: each generated covariate with its distribution and parameters.
#'   \item \code{files}: paths to on-disk files (csv/rds/rdata) when available.
#' }
#' @export
#' @examples
#' \dontrun{
#' sets <- load_recipe_sets("checks/manifest.rds")
#' spec <- describe_generation(sets[[1]])
#' spec$header
#' spec$baseline
#' spec$effects
#' spec$treatment
#' spec$censoring
#' spec$covariates
#' spec$files
#' }
describe_generation <- function(set) {
   `%||%` <- function(x, y) if (is.null(x)) y else x

   .pretty_model <- function(m) {
      switch(tolower(m),
             "cox_pwexp"       = "PH piecewise-exponential",
             "ph_exponential"  = "PH exponential",
             "ph_weibull"      = "PH Weibull",
             "ph_gompertz"     = "PH Gompertz",
             "aft_lognormal"   = "AFT log-normal",
             "aft_weibull"     = "AFT Weibull",
             "aft_loglogistic" = "AFT log-logistic",
             m
      )
   }

   .flatten_params <- function(p) {
      if (is.null(p) || !length(p)) return(NA_character_)
      paste(sprintf("%s=%s", names(p), vapply(p, function(v){
         if (is.list(v)) return("[list]")
         if (is.numeric(v) && length(v) > 6) return(sprintf("[%.0f vals]", length(v)))
         paste(v, collapse = ",")
      }, character(1))), collapse = "; ")
   }

   describe_covariates <- function(meta) {
      defs <- meta$covariates %||% list()
      if (!length(defs)) return(data.frame())
      do.call(rbind, lapply(defs, function(d) {
         data.frame(
            name       = d$name %||% NA_character_,
            type       = d$type %||% NA_character_,
            dist       = d$dist %||% NA_character_,
            parameters = .flatten_params(d$params),
            transform  = if (length(d$transform)) paste(d$transform, collapse = " + ") else NA_character_,
            stringsAsFactors = FALSE, check.names = FALSE
         )
      }))
   }

   describe_effects <- function(effects, data = NULL) {
      if (is.null(effects)) return(data.frame())
      rows <- list()
      rows[["(Intercept)"]] <- effects$intercept %||% 0
      rows[["treatment"]]   <- effects$treatment %||% 0

      # (A) named covariate coefficients
      if (!is.null(effects$covariates)) {
         for (nm in names(effects$covariates)) rows[[nm]] <- effects$covariates[[nm]]
      }

      # (B) formula + beta → model.matrix columns
      if (!is.null(effects$formula) && !is.null(effects$beta)) {
         X0 <- if (!is.null(data)) data else data.frame()
         mm <- stats::model.matrix(stats::as.formula(effects$formula), data = X0)
         b  <- effects$beta
         if (length(b) != ncol(mm)) stop("effects$beta length mismatch with model.matrix columns")
         for (j in seq_along(b)) rows[[colnames(mm)[j]]] <- b[[j]]
      }

      data.frame(term = names(rows), coefficient = unlist(rows), row.names = NULL, check.names = FALSE)
   }

   describe_treatment <- function(meta) {
      tr <- meta$treatment %||% list()
      mode <- tr$assignment %||% NA_character_
      out <- list(
         assignment = mode,
         allocation = tr$allocation %||% NA_character_
      )
      if (identical(mode, "stratified")) out$stratify_by <- paste(tr$stratify_by, collapse = ", ")
      if (identical(mode, "logistic_ps")) {
         out$ps_formula <- if (!is.null(tr$ps_model$formula)) paste(deparse(tr$ps_model$formula), collapse = "") else NA_character_
         out$ps_beta    <- if (!is.null(tr$ps_model$beta)) paste(tr$ps_model$beta, collapse = ",") else NA_character_
      }
      as.data.frame(out, check.names = FALSE)
   }

   describe_censoring <- function(meta) {
      cz <- meta$censoring %||% list()
      data.frame(
         mode        = cz$mode %||% NA_character_,
         target      = cz$target %||% NA_real_,
         admin_time  = cz$admin_time %||% NA_real_,
         stringsAsFactors = FALSE, check.names = FALSE
      )
   }

   # ---- main body ----
   d <- set$data
   m <- set$meta

   L <- m$tau %||% attr(d, "tau") %||% NA_real_   # shown as "L" in docs
   header <- data.frame(
      n                  = NROW(d),
      L                  = L,
      model              = .pretty_model(m$model %||% NA_character_),
      event_rate         = mean(d$status == 1L, na.rm = TRUE),
      achieved_censoring = (m$achieved_censoring %||% mean(d$status == 0L, na.rm = TRUE)),
      stringsAsFactors = FALSE, check.names = FALSE
   )

   baseline <- data.frame(
      baseline = .flatten_params(m$baseline),
      stringsAsFactors = FALSE, check.names = FALSE
   )

   covs_tbl <- describe_covariates(m)
   eff_tbl  <- describe_effects(m$effects %||% list(), data = d)
   tr_tbl   <- describe_treatment(m)
   cens_tbl <- describe_censoring(m)

   files <- data.frame(
      csv   = m$files$csv   %||% NA_character_,
      rds   = m$files$rds   %||% NA_character_,
      rdata = m$files$rdata %||% NA_character_,
      stringsAsFactors = FALSE, check.names = FALSE
   )

   list(
      header    = header,
      baseline  = baseline,
      effects   = eff_tbl,
      treatment = tr_tbl,
      censoring = cens_tbl,
      covariates = covs_tbl,
      files     = files
   )
}

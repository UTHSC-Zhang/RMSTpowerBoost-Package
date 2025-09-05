
# R/recipe_convenience.R

#' Quick AFT recipe builder (list-only)
#'
#' Convenience constructor for common AFT scenarios (lognormal/Weibull).
#' @param n Sample size.
#' @param tau Restriction time (analysis horizon).
#' @param model One of \code{"aft_lognormal"} or \code{"aft_weibull"}.
#' @param baseline Baseline parameter list (see model).
#' @param treat_effect Numeric treatment coefficient (on log-time scale).
#' @param covariates Covariate definitions (list of defs).
#' @param target_censoring Target overall censoring fraction (0-1).
#' @param allocation Allocation ratio string (e.g., "1:1").
#' @param seed Optional seed.
#' @return A recipe list suitable for \code{\link{simulate_from_recipe}}.
#' @examples
#' covs <- list(list(name="x", type="continuous", dist="normal", params=list(mean=0, sd=1)))
#' r <- recipe_quick_aft(120, 12, "aft_lognormal",
#'        baseline=list(mu=2.3, sigma=0.5), treat_effect=-0.2,
#'        covariates=covs, target_censoring=0.25, allocation="1:1")
#' dat <- simulate_from_recipe(r, seed = 1)
#' @export
recipe_quick_aft <- function(n, tau, model = c("aft_lognormal","aft_weibull"),
                             baseline, treat_effect, covariates,
                             target_censoring = 0.25, allocation = "1:1",
                             seed = NULL) {
  model <- match.arg(model)
  list(
    n = as.integer(n),
    covariates = list(defs = covariates),
    treatment = list(assignment = "randomization", allocation = allocation),
    event_time = list(
      model = model,
      baseline = baseline,
      effects = list(intercept = 0, treatment = treat_effect, covariates = NULL),
      tau = tau
    ),
    censoring = list(mode = "target_overall", target = target_censoring, admin_time = Inf),
    seed = seed
  )
}

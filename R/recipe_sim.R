
# R/recipe_sim.R
# List-only simulation engine (no YAML). Supports AFT/PH parametric models,
# piecewise-exponential, optional frailty, richer covariate distributions,
# explicit intercept handling, and robust censoring calibration.

#' Validate a simulation recipe (list-only schema)
#'
#' Ensures the list schema is consistent and fills reasonable defaults.
#' @param recipe A named list (see the vignette "Data Generation Scheme").
#' @return A validated recipe list.
#' @examples
#' r <- recipe_quick_aft(100, 12, "aft_lognormal",
#'   baseline = list(mu=2.2, sigma=0.5),
#'   treat_effect = -0.2,
#'   covariates = list(list(name="x", type="continuous", dist="normal", params=list(mean=0, sd=1))),
#'   target_censoring = 0.25, allocation = "1:1")
#' r2 <- validate_recipe(r)
#' @export
validate_recipe <- function(recipe) {
  if (is.null(recipe) || !is.list(recipe)) stop("`recipe` must be a list (YAML is no longer supported).")
  # n
  if (is.null(recipe$n) || !is.numeric(recipe$n) || length(recipe$n) != 1L || recipe$n <= 0)
    stop("`recipe$n` must be a positive number.")
  recipe$n <- as.integer(recipe$n)

  # covariates
  if (is.null(recipe$covariates) || is.null(recipe$covariates$defs) || !is.list(recipe$covariates$defs))
    recipe$covariates <- list(defs = list())

  # treatment (optional)
  if (!is.null(recipe$treatment)) {
    tr <- recipe$treatment
    if (is.null(tr$assignment)) stop("treatment$assignment must be provided when treatment block exists.")
    if (!tr$assignment %in% c("randomization","stratified","logistic_ps"))
      stop("Unknown treatment$assignment: ", tr$assignment)
    if (tr$assignment %in% c("randomization","stratified")) {
      tr$allocation <- tr$allocation %||% "1:1"
      recipe$treatment$allocation <- tr$allocation
      if (tr$assignment == "stratified") {
        if (is.null(tr$stratify_by)) stop("treatment$stratify_by required for stratified assignment.")
      }
    } else { # logistic_ps
      if (is.null(tr$ps_model$formula)) stop("treatment$ps_model$formula is required for logistic_ps.")
    }
  }

  # event_time
  et <- recipe$event_time
  if (is.null(et) || is.null(et$model)) stop("event_time$model must be provided.")
  allowed <- c("aft_lognormal","aft_weibull","aft_loglogistic",
               "ph_exponential","ph_weibull","ph_gompertz","cox_pwexp")
  if (!et$model %in% allowed) stop("Unsupported event_time$model: ", et$model)
  if (is.null(et$baseline) || !is.list(et$baseline)) stop("event_time$baseline must be a list.")
  et$effects <- et$effects %||% list()
  et$effects$intercept <- et$effects$intercept %||% 0
  et$effects$treatment <- et$effects$treatment %||% 0
  if (!is.null(et$effects$formula) && is.null(et$effects$beta))
    stop("effects$formula provided but effects$beta missing.")
  if (!is.null(et$effects$beta) && is.null(et$effects$formula))
    stop("effects$beta provided but effects$formula missing.")
  et$tau <- et$tau %||% Inf
  recipe$event_time <- et

  # frailty (optional)
  if (!is.null(et$frailty)) {
    fr <- et$frailty
    if (is.null(fr$type) || !fr$type %in% c("gamma","lognormal"))
      stop("frailty$type must be 'gamma' or 'lognormal'.")
    if (is.null(fr$var) || !is.numeric(fr$var) || fr$var < 0)
      stop("frailty$var must be a nonnegative number (variance parameter).")
    if (is.null(fr$group) || !is.character(fr$group) || length(fr$group)!=1L)
      stop("frailty$group must be a single existing covariate name.")
  }

  # censoring
  cz <- recipe$censoring %||% list(mode = "target_overall", target = 0.25, admin_time = Inf)
  if (is.null(cz$mode) || !cz$mode %in% c("target_overall","explicit"))
    stop("censoring$mode must be 'target_overall' or 'explicit'.")
  if (cz$mode == "target_overall") {
    cz$target <- cz$target %||% 0.25
    cz$admin_time <- cz$admin_time %||% Inf
  } else {
    # explicit: any of administrative, random, dependent may be present
    # e.g., random=list(dist="exponential", params=list(rate=0.03))
  }
  recipe$censoring <- cz

  recipe$seed <- recipe$seed %||% NULL
  recipe
}

# null-coalescing helper
`%||%` <- function(x, y) if (is.null(x)) y else x

# ---------- Covariates ----------
# Apply simple "center(a)" and "scale(b)" transforms
.apply_transforms <- function(x, tf) {
  if (is.null(tf)) return(x)
  for (t in tf) {
    if (grepl("^center\\(", t)) {
      a <- as.numeric(sub(".*center\\(([^)]+)\\).*", "\\1", t))
      x <- (x - a)
    } else if (grepl("^scale\\(", t)) {
      b <- as.numeric(sub(".*scale\\(([^)]+)\\).*", "\\1", t))
      x <- x / b
    }
  }
  x
}

.rdraw <- function(def, n) {
  p <- def$params %||% list()
  switch(def$dist,
    normal    = stats::rnorm(n, p$mean %||% 0, p$sd %||% 1),
    lognormal = stats::rlnorm(n, p$meanlog %||% 0, p$sdlog %||% 1),
    gamma     = stats::rgamma(n, shape = p$shape, scale = p$scale %||% 1),
    weibull   = stats::rweibull(n, shape = p$shape, scale = p$scale %||% 1),
    uniform   = stats::runif(n, p$min %||% 0, p$max %||% 1),
    beta      = stats::rbeta(n, p$shape1, p$shape2),
    t         = stats::rt(n, df = p$df %||% 5),
    bernoulli = stats::rbinom(n, 1, p$p %||% 0.5),
    categorical = {
      k <- length(p$prob); idx <- base::sample.int(k, n, replace = TRUE, prob = p$prob)
      if (!is.null(p$labels)) factor(p$labels[idx], levels = p$labels) else idx
    },
    ordinal = {
      k <- length(p$prob); idx <- base::sample.int(k, n, replace = TRUE, prob = p$prob)
      labs <- p$labels %||% as.character(seq_len(k))
      factor(labs[idx], levels = labs, ordered = TRUE)
    },
    stop("Unsupported covariate dist: ", def$dist)
  )
}

#' Generate covariate matrix/data frame from a recipe
#' @keywords internal
#' @param n sample size
#' @param covariates list(defs = list(...))
#' @return data.frame of covariates
#' @examples
#' \dontrun{
#' # Example covariate defs:
#' defs <- list(
#'   list(name = "x", type = "continuous", dist = "normal",
#'        params = list(mean = 0, sd = 1)),
#'   list(
#'     name = "z", type = "categorical", dist = "categorical",
#'     params = list(prob = c(0.3, 0.7), labels = c("A", "B"))
#'   )
#' )
#' X <- gen_covariates(10, list(defs = defs))

#' }
gen_covariates <- function(n, covariates) {
  defs <- covariates$defs %||% list()
  if (!length(defs)) return(data.frame())
  out <- vector("list", length(defs))
  names(out) <- vapply(defs, function(d) d$name, character(1))
  for (i in seq_along(defs)) {
    d <- defs[[i]]
    xi <- .rdraw(d, n)
    xi <- .apply_transforms(xi, d$transform)
    out[[i]] <- xi
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

# ---------- Treatment assignment ----------
.parse_allocation <- function(s) {
  parts <- strsplit(s, ":", fixed = TRUE)[[1]]
  if (length(parts) != 2) stop("allocation must be of the form 'a:b'")
  a <- as.numeric(parts[1]); b <- as.numeric(parts[2])
  p1 <- a/(a+b); c(p0 = 1-p1, p1 = p1)
}

.assign_treatment <- function(n, X, tr) {
  if (is.null(tr)) return(NULL)
  if (tr$assignment == "randomization") {
    p <- .parse_allocation(tr$allocation %||% "1:1")
    return(stats::rbinom(n, 1, prob = p["p1"]))
  }
  if (tr$assignment == "stratified") {
    p <- .parse_allocation(tr$allocation %||% "1:1")
    z <- X[, tr$stratify_by, drop = FALSE]
    # ensure factors for strata
    for (nm in names(z)) if (!is.factor(z[[nm]])) z[[nm]] <- factor(z[[nm]])
    strata <- interaction(z, drop = TRUE)
    arm <- integer(n)
    for (lv in levels(strata)) {
      idx <- which(strata == lv)
      arm[idx] <- stats::rbinom(length(idx), 1, prob = p["p1"])
    }
    return(arm)
  }
  if (tr$assignment == "logistic_ps") {
    form <- tr$ps_model$formula
    mm <- stats::model.matrix(stats::as.formula(form), data = X)
    eta <- as.vector(mm %*% rep(1, ncol(mm))) # default if no coefficients provided
    if (!is.null(tr$ps_model$beta)) {
      if (length(tr$ps_model$beta) != ncol(mm)) stop("ps_model$beta length mismatch with model.matrix")
      eta <- as.vector(mm %*% tr$ps_model$beta)
    }
    p1 <- stats::plogis(eta)
    return(stats::rbinom(n, 1, p1))
  }
  stop("Unknown treatment assignment.")
}

# ---------- Linear predictor and frailty ----------
.build_lp <- function(effects, X, arm) {
  b0 <- effects$intercept %||% 0
  bt <- effects$treatment %||% 0
  lp_cov <- 0
  if (!is.null(effects$formula) && !is.null(effects$beta)) {
    mm <- stats::model.matrix(stats::as.formula(effects$formula), data = X)
    if (length(effects$beta) != ncol(mm))
      stop("length(effects$beta) must equal ncol(model.matrix)")
    lp_cov <- as.vector(mm %*% effects$beta)
  } else if (!is.null(effects$covariates)) {
    for (nm in names(effects$covariates)) {
      if (is.null(X[[nm]])) stop("effects refers to unknown covariate: ", nm)
      lp_cov <- lp_cov + effects$covariates[[nm]] * X[[nm]]
    }
  }
  b0 + bt * (arm %||% 0) + lp_cov
}

.apply_frailty <- function(eta, model, frailty, groups) {
  if (is.null(frailty)) return(eta)
  type <- frailty$type; v <- frailty$var; gname <- frailty$group
  if (is.null(groups[[gname]])) stop("frailty group variable '", gname, "' not found in data.")
  g <- groups[[gname]]; ug <- unique(g)
  if (type == "gamma") {
    if (!grepl("^ph_|^cox_", model)) stop("gamma frailty supported for PH/Cox models only.")
    shape <- 1/v; scale <- v
    vj <- stats::rgamma(length(ug), shape = shape, scale = scale); names(vj) <- as.character(ug)
    eta + log(vj[as.character(g)])
  } else if (type == "lognormal") {
    zj <- stats::rnorm(length(ug), mean = 0, sd = sqrt(v)); names(zj) <- as.character(ug)
    eta + zj[as.character(g)]
  } else {
    stop("Unsupported frailty type: ", type)
  }
}

# ---------- Event-time simulators ----------
.sim_pwexp <- function(rates, cuts, lp) {
  # rates per segment; multiply by exp(lp)
  n <- length(lp)
  t <- numeric(n); alive <- rep(TRUE, n)
  last <- 0
  seg_rates <- c(rates, tail(rates,1))
  seg_cuts  <- c(cuts, Inf)
  for (s in seq_along(seg_rates)) {
    idx <- which(alive)
    if (!length(idx)) break
    r <- seg_rates[s] * exp(lp[idx])
    dt <- stats::rexp(length(idx), rate = r)
    t[idx] <- t[idx] + dt
    alive_idx <- idx[ t[idx] < seg_cuts[s] ]
    # those exceeding the cut move to next segment; time carries over
    alive <- rep(FALSE, n); alive[alive_idx] <- TRUE
    # clamp times at the cut for those who cross, then continue loop
    cross_idx <- setdiff(idx, alive_idx)
    if (length(cross_idx)) t[cross_idx] <- seg_cuts[s]
  }
  t
}

.sim_time <- function(model, baseline, eta, n) {
  if (model == "aft_lognormal") {
    mu <- baseline$mu; sigma <- baseline$sigma
    return( exp(mu + eta + stats::rnorm(n, 0, sigma)) )
  }
  if (model == "aft_weibull") {
    k <- baseline$shape; lam <- baseline$scale
    u <- stats::runif(n)
    return( lam * exp(eta) * (-log(1 - u))^(1/k) )
  }
  if (model == "aft_loglogistic") {
    k <- baseline$shape; lam <- baseline$scale
    u <- stats::runif(n)
    return( lam * exp(eta) * (u/(1-u))^(1/k) )
  }
  if (model == "ph_exponential") {
    lam0 <- baseline$rate
    return( stats::rexp(n, rate = lam0 * exp(eta)) )
  }
  if (model == "ph_weibull") {
    k <- baseline$shape; lam <- baseline$scale
    u <- stats::runif(n)
    return( lam * (-log(u) / exp(eta))^(1/k) )
  }
  if (model == "ph_gompertz") {
    a <- baseline$rate; b <- baseline$gamma
    u <- stats::runif(n)
    return( (1/b) * log(1 - (b/a) * log(u) / exp(eta)) )
  }
  if (model == "cox_pwexp") {
    rates <- baseline$rates; cuts <- baseline$cuts %||% numeric(0)
    return( .sim_pwexp(rates, cuts, lp = eta) )
  }
  stop("Unsupported model: ", model)
}

# ---------- Censoring ----------
.achieved_cens_exp <- function(T_event, rate, admin_time = Inf) {
  n <- length(T_event)
  admin <- if (is.null(admin_time) || !is.finite(admin_time)) Inf else admin_time
  if (is.na(rate) || rate <= 0) {
    C_rand <- rep(Inf, n)
  } else {
    C_rand <- stats::rexp(n, rate = rate)
  }
  C <- pmin(C_rand, admin)
  status <- as.integer(T_event <= C)
  mean(status == 0, na.rm = TRUE)
}

.solve_rate_for_target <- function(T_event, target, admin_time = Inf, tol = 0.002) {
  admin <- if (is.null(admin_time) || !is.finite(admin_time)) Inf else admin_time
  target <- max(min(target, 0.999), 0.001)
  cmin <- mean(T_event > admin, na.rm = TRUE)
  if (!is.finite(cmin) || is.na(cmin)) cmin <- 0
  if (cmin >= target - tol) return(0)
  f <- function(rate) .achieved_cens_exp(T_event, rate, admin) - target
  lo <- 1e-12; hi <- 0.1
  for (k in 1:60) {
    val <- f(hi)
    if (is.finite(val) && val >= 0) break
    hi <- hi * 2
    if (hi > 1e6) break
  }
  val_hi <- f(hi)
  if (!is.finite(val_hi) || val_hi < 0) return(hi)
  uniroot(function(r) f(r), c(lo, hi), tol = tol)$root
}

#' Simulate a dataset from a validated recipe (list-only)
#' @param recipe A validated recipe list (use \code{validate_recipe()}).
#' @param seed Optional integer seed to override recipe$seed.
#' @return A data.frame with columns time, status, arm (if treatment present), plus covariates.
#'         Attributes: "tau", "achieved_censoring".
#' @examples
#' covs <- list(list(name="x", type="continuous", dist="normal", params=list(mean=0, sd=1)))
#' rec <- recipe_quick_aft(120, 12, "aft_lognormal",
#'         baseline=list(mu=2.2, sigma=0.5), treat_effect=-0.2,
#'         covariates=covs, target_censoring=0.25)
#' dat <- simulate_from_recipe(rec, seed = 11)
#' @export
simulate_from_recipe <- function(recipe, seed = NULL) {
  if (!is.list(recipe)) stop("`recipe` must be a list (no YAML paths).")
  recipe <- validate_recipe(recipe)
  if (!is.null(seed)) {
    old <- .Random.seed; on.exit({ if (exists("old")) .Random.seed <<- old }, add = TRUE)
    set.seed(seed)
  } else if (!is.null(recipe$seed)) {
    old <- .Random.seed; on.exit({ if (exists("old")) .Random.seed <<- old }, add = TRUE)
    set.seed(recipe$seed)
  }
  n <- recipe$n

  # Covariates
  X <- gen_covariates(n, recipe$covariates)

  # Treatment
  arm <- .assign_treatment(n, X, recipe$treatment %||% NULL)
  if (!is.null(arm)) X$arm <- arm

  # Linear predictor + frailty
  eta <- .build_lp(recipe$event_time$effects, X, arm)
  # If frailty references a group var, ensure it's in X (e.g., cluster)
  eta <- .apply_frailty(eta, model = recipe$event_time$model,
                        frailty = recipe$event_time$frailty %||% NULL,
                        groups = X)

  # Event times
  T_event <- .sim_time(recipe$event_time$model, recipe$event_time$baseline, eta, n)

  # Censoring
  cz <- recipe$censoring
  admin <- cz$admin_time %||% Inf
  if (cz$mode == "target_overall") {
    rate <- .solve_rate_for_target(T_event, cz$target, admin_time = admin, tol = 0.002)
    C_rand <- if (rate <= 0) rep(Inf, n) else stats::rexp(n, rate = rate)
    C <- pmin(C_rand, admin)
  } else {
    C <- rep(Inf, n)
    if (!is.null(cz$administrative$time)) {
      C <- pmin(C, cz$administrative$time)
    }
    if (!is.null(cz$random)) {
      rd <- cz$random
      if (rd$dist == "exponential") {
        C <- pmin(C, stats::rexp(n, rate = rd$params$rate))
      } else {
        stop("explicit random censoring supports dist='exponential' only for now.")
      }
    }
    if (!is.null(cz$dependent)) {
      # simple linear predictor on covariates: hazard = base * exp(formula(X))
      f <- cz$dependent$formula %||% NULL
      if (!is.null(f)) {
        mm <- stats::model.matrix(stats::as.formula(f), data = X)
        gamma <- cz$dependent$beta %||% rep(0, ncol(mm))
        if (length(gamma) != ncol(mm)) stop("censoring dependent beta length mismatch.")
        lp_c <- as.vector(mm %*% gamma)
        base <- cz$dependent$base %||% 0.02
        C <- pmin(C, stats::rexp(n, rate = base * exp(lp_c)))
      }
    }
  }

  time <- pmin(T_event, C)
  status <- as.integer(T_event <= C)

  out <- data.frame(time = as.numeric(time), status = as.integer(status))
  if (!is.null(arm)) out$arm <- as.integer(arm)
  # bind covariates last for readability
  if (ncol(X)) {
    # drop helper 'arm' from X if duplicated
    if ("arm" %in% names(X)) X$arm <- NULL
    out <- cbind(out, X)
  }

  attr(out, "tau") <- recipe$event_time$tau %||% Inf
  attr(out, "achieved_censoring") <- mean(out$status == 0)
  out
}

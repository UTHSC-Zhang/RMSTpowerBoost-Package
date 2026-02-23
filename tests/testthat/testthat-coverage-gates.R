make_small_pilot <- function(n = 40, stratified = FALSE) {
  set.seed(777)
  dat <- data.frame(
    time = stats::rexp(n, 0.2) + 0.1,
    status = stats::rbinom(n, 1, 0.7),
    arm = rep(0:1, length.out = n),
    x1 = stats::rnorm(n)
  )
  if (stratified) {
    dat$region <- factor(rep(c("A", "B"), length.out = n))
  }
  dat
}

test_that("run_app covers env override and missing app path branches", {
  testthat::skip_if_not_installed("shiny")

  shiny_ns <- asNamespace("shiny")
  orig_run_app <- get("runApp", envir = shiny_ns)
  called <- NULL
  assignInNamespace("runApp", function(appDir, ...) {
    called <<- appDir
    appDir
  }, ns = "shiny")
  on.exit(assignInNamespace("runApp", orig_run_app, ns = "shiny"), add = TRUE)

  td <- tempfile("appdir_")
  dir.create(td, recursive = TRUE, showWarnings = FALSE)
  old <- Sys.getenv("RMSTPOWERBOOST_APP_DIR", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("RMSTPOWERBOOST_APP_DIR") else Sys.setenv(RMSTPOWERBOOST_APP_DIR = old)
  }, add = TRUE)

  Sys.setenv(RMSTPOWERBOOST_APP_DIR = td)
  expect_identical(
    testthat::with_mocked_bindings(
      RMSTpowerBoost::run_app(),
      .ensure_app_dependencies = function(...) invisible(TRUE),
      .package = "RMSTpowerBoost"
    ),
    td
  )
  expect_identical(called, td)

  Sys.setenv(RMSTPOWERBOOST_APP_DIR = file.path(tempdir(), "does-not-exist-rmst-app"))
  expect_error(
    testthat::with_mocked_bindings(
      RMSTpowerBoost::run_app(),
      .ensure_app_dependencies = function(...) invisible(TRUE),
      .package = "RMSTpowerBoost"
    ),
    "Could not find the app directory"
  )
})

test_that("describe_generation covers model labels, treatment modes, and effects branches", {
  dat <- make_small_pilot(20, stratified = TRUE)
  attr(dat, "tau") <- 2.5

  model_keys <- c("cox_pwexp", "ph_exponential", "ph_weibull", "ph_gompertz", "aft_lognormal", "aft_weibull", "aft_loglogistic")
  for (mk in model_keys) {
    set <- list(
      data = dat,
      meta = list(
        model = mk,
        baseline = list(rate = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), nested = list(a = 1)),
        effects = list(intercept = 0.1, treatment = -0.2, covariates = list(x1 = 0.5)),
        treatment = list(assignment = "stratified", allocation = "1:1", stratify_by = c("region")),
        censoring = list(mode = "target_overall", target = 0.2, admin_time = 5),
        covariates = list(list(name = "x1", type = "continuous", dist = "normal", params = list(mean = 0, sd = 1), transform = c("center(0)"))
        )
      )
    )
    out <- RMSTpowerBoost::describe_generation(set)
    expect_true(is.data.frame(out$header))
    expect_true(is.data.frame(out$treatment))
    expect_true(any(grepl("stratified", out$treatment$assignment)))
    expect_true(any(grepl("\\[7 vals\\]|\\[list\\]", out$baseline$baseline)))
  }

  set_logit <- list(
    data = dat,
    meta = list(
      model = "aft_lognormal",
      baseline = list(mu = 2, sigma = 0.5),
      effects = list(intercept = 0, treatment = 0),
      treatment = list(assignment = "logistic_ps", allocation = "1:1", ps_model = list(formula = ~x1, beta = c(0.1, 0.2))),
      censoring = list(mode = "explicit")
    )
  )
  out_logit <- RMSTpowerBoost::describe_generation(set_logit)
  expect_true("ps_formula" %in% names(out_logit$treatment))
  expect_true("ps_beta" %in% names(out_logit$treatment))

  # effects$formula + beta mismatch branch
  set_bad <- set_logit
  set_bad$meta$effects <- list(intercept = 0, treatment = 0, formula = "~ x1", beta = c(0.1, 0.2, 0.3))
  expect_error(RMSTpowerBoost::describe_generation(set_bad), "effects\\$beta length mismatch")
})

test_that("load/rebuild manifest covers csv branch, no-path branch, and vary token branch", {
  td <- tempfile("manifest_cov_")
  dir.create(td, recursive = TRUE, showWarnings = FALSE)

  dat <- make_small_pilot(20, stratified = FALSE)
  csv_path <- file.path(td, "d.csv")
  utils::write.csv(dat, csv_path, row.names = FALSE)

  man <- data.frame(
    scenario_id = 1, rep = 1, seed = 1L,
    file_rds = NA_character_, file_rdata = NA_character_, file_csv = csv_path, file_txt = NA_character_,
    stringsAsFactors = FALSE
  )
  saveRDS(man, file.path(td, "manifest.rds"))
  sets <- RMSTpowerBoost::load_recipe_sets(file.path(td, "manifest.rds"))
  expect_equal(length(sets), 1)
  expect_true(is.data.frame(sets[[1]]$data))

  bad <- man
  bad$file_csv <- NA_character_
  saveRDS(bad, file.path(td, "manifest.rds"))
  expect_error(RMSTpowerBoost::load_recipe_sets(file.path(td, "manifest.rds")), "No data file path")

  empty <- man[0, ]
  saveRDS(empty, file.path(td, "manifest.rds"))
  expect_error(RMSTpowerBoost::load_recipe_sets(file.path(td, "manifest.rds")), "manifest has 0 rows")

  rec <- RMSTpowerBoost::recipe_quick_aft(
    n = 10,
    model = "aft_lognormal",
    baseline = list(mu = 2, sigma = 0.4),
    treat_effect = -0.1,
    covariates = list(list(name = "x1", type = "continuous", dist = "normal", params = list(mean = 0, sd = 1))),
    target_censoring = 0.1,
    allocation = "1:1"
  )
  d1 <- RMSTpowerBoost::simulate_from_recipe(rec, seed = 1)
  saveRDS(d1, file.path(td, "sc1_r1_n10.rds"))
  rb <- RMSTpowerBoost::rebuild_manifest(
    rec,
    vary = list(n = c(10, 20)),
    out_dir = td,
    filename_template = "sc{scenario_id}_r{rep}_n{n}"
  )
  expect_true(nrow(rb) >= 1)
  expect_true(any(grepl("^p__", names(rb))))
})

test_that("rebuild_manifest covers missing-dir, no-match, and file-format branches", {
  rec <- RMSTpowerBoost::recipe_quick_aft(
    n = 12,
    model = "aft_lognormal",
    baseline = list(mu = 2, sigma = 0.4),
    treat_effect = -0.1,
    covariates = list(list(name = "x1", type = "continuous", dist = "normal", params = list(mean = 0, sd = 1))),
    target_censoring = 0.1,
    allocation = "1:1"
  )

  missing_dir <- file.path(tempdir(), "rmst_missing_dir_cov")
  expect_error(
    RMSTpowerBoost::rebuild_manifest(rec, vary = list(), out_dir = missing_dir),
    "out_dir does not exist"
  )

  empty_dir <- tempfile("rmst_empty_rebuild_")
  dir.create(empty_dir, recursive = TRUE, showWarnings = FALSE)
  expect_error(
    RMSTpowerBoost::rebuild_manifest(rec, vary = list(), out_dir = empty_dir),
    "no files matched"
  )

  dat <- RMSTpowerBoost::simulate_from_recipe(rec, seed = 5)
  fmt_dir <- tempfile("rmst_fmt_rebuild_")
  dir.create(fmt_dir, recursive = TRUE, showWarnings = FALSE)

  utils::write.csv(dat, file.path(fmt_dir, "sc1_r1.csv"), row.names = FALSE)
  rb_csv <- RMSTpowerBoost::rebuild_manifest(rec, vary = list(), out_dir = fmt_dir)
  expect_true(nrow(rb_csv) >= 1)

  unlink(file.path(fmt_dir, "manifest.rds"), force = TRUE)
  dat_obj <- dat
  save(dat_obj, file = file.path(fmt_dir, "sc1_r1.RData"))
  rb_rdata <- RMSTpowerBoost::rebuild_manifest(rec, vary = list(), out_dir = fmt_dir)
  expect_true(nrow(rb_rdata) >= 1)

  unlink(file.path(fmt_dir, "manifest.rds"), force = TRUE)
  utils::write.table(dat, file.path(fmt_dir, "sc1_r1.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  rb_txt <- RMSTpowerBoost::rebuild_manifest(rec, vary = list(), out_dir = fmt_dir)
  expect_true(nrow(rb_txt) >= 1)

  # Vary with zero-length value can produce empty grid and triggers fallback to base recipe.
  rb_empty_grid <- RMSTpowerBoost::rebuild_manifest(
    rec,
    vary = list("event_time.effects.treatment" = numeric(0)),
    out_dir = fmt_dir
  )
  expect_true(nrow(rb_empty_grid) >= 1)
})

test_that("linear boot functions cover validation and no-estimate summary branches", {
  pilot <- make_small_pilot(20, stratified = FALSE)
  expect_error(RMSTpowerBoost::linear.power.boot(pilot, "time", "status", "arm", sample_sizes = NULL, L = 2), "sample_sizes")
  expect_error(RMSTpowerBoost::linear.ss.boot(pilot, "time", "status", "arm", target_power = c(0.8, 0.9), L = 2), "single numeric value")

  # Force null-estimate branches by using nonnumeric event times that make survfit fail.
  bad <- pilot
  bad$time <- as.character(round(bad$time, 3))
  p <- RMSTpowerBoost::linear.power.boot(bad, "time", "status", "arm", sample_sizes = c(6), L = 2, n_sim = 2, alpha = 0.1)
  expect_true(is.null(p$results_summary))

  s <- suppressWarnings(RMSTpowerBoost::linear.ss.boot(
    bad, "time", "status", "arm",
    target_power = 0.95, L = 2, n_sim = 2, alpha = 0.1,
    patience = 1, n_start = 4, n_step = 2, max_n_per_arm = 6
  ))
  expect_true(is.null(s$results_summary))
})

test_that("multiplicative analytical and bootstrap cover additional warning/error branches", {
  dat <- make_small_pilot(60, stratified = TRUE)
  expect_warning(
    RMSTpowerBoost::MS.ss.analytical(
      dat, "time", "status", "arm", "region",
      target_power = 0.999, L = 2, n_start = 4, n_step = 2, max_n_per_arm = 6
    ),
    "not achieved by max N"
  )

  # low-data guard in analytical power
  tiny <- dat[1:6, ]
  tiny$status <- c(1, 0, 0, 0, 0, 0)
  expect_error(
    RMSTpowerBoost::MS.power.analytical(tiny, "time", "status", "arm", "region", sample_sizes = c(4), L = 2),
    "Not enough data points"
  )

  expect_error(RMSTpowerBoost::MS.power.boot(dat, "time", "status", "arm", "region", sample_sizes = NULL, L = 2), "sample_sizes")
  expect_error(RMSTpowerBoost::MS.ss.boot(dat, "time", "status", "arm", "region", target_power = c(0.8, 0.9), L = 2), "single numeric value")

  # Trigger max_n warning branch rather than stagnation.
  expect_warning(
    RMSTpowerBoost::MS.ss.boot(
      dat, "time", "status", "arm", "region",
      target_power = 0.999, L = 2, n_sim = 2, alpha = 0.1,
      parallel.cores = 1, patience = 99, n_start = 4, n_step = 2, max_n_per_arm = 6
    ),
    "not achieved by max N"
  )
})

test_that("GAM bootstrap functions cover validation, parallel, and search terminal branches", {
  dat <- make_small_pilot(30, stratified = TRUE)

  expect_error(
    RMSTpowerBoost::GAM.power.boot(dat, "time", "status", "arm", sample_sizes = NULL, L = 2),
    "sample_sizes"
  )
  expect_error(
    RMSTpowerBoost::GAM.ss.boot(dat, "time", "status", "arm", target_power = NULL, L = 2),
    "target_power"
  )

  # Parallel dependency guard branch.
  testthat::with_mocked_bindings(
    expect_error(
      RMSTpowerBoost::GAM.power.boot(
        dat, "time", "status", "arm",
        sample_sizes = c(4), L = 2, parallel.cores = 2
      ),
      "required for parallel processing"
    ),
    requireNamespace = function(...) FALSE,
    .package = "base"
  )
  testthat::with_mocked_bindings(
    expect_error(
      RMSTpowerBoost::GAM.ss.boot(
        dat, "time", "status", "arm",
        target_power = 0.8, L = 2, parallel.cores = 2
      ),
      "required for parallel processing"
    ),
    requireNamespace = function(...) FALSE,
    .package = "base"
  )

  # Exercise pseudo-observation n==1 and n==0 branches (expected unstable fit errors).
  expect_error(
    RMSTpowerBoost::GAM.power.boot(
      dat, "time", "status", "arm", strata_var = "region",
      sample_sizes = c(1), linear_terms = "x1",
      L = 2, n_sim = 1, alpha = 0.1, parallel.cores = 1
    ),
    "infinite or missing|Canceling all iterations"
  )
  expect_error(
    RMSTpowerBoost::GAM.power.boot(
      dat, "time", "status", "arm", strata_var = "region",
      sample_sizes = c(0), linear_terms = "x1",
      L = 2, n_sim = 1, alpha = 0.1, parallel.cores = 1
    ),
    "second argument must be a list|Canceling all iterations"
  )

  # Exercise explicit parallel plan path with a stable small sample size.
  gp <- RMSTpowerBoost::GAM.power.boot(
    dat, "time", "status", "arm", strata_var = "region",
    sample_sizes = c(4), linear_terms = "x1",
    L = 2, n_sim = 1, alpha = 0.1, parallel.cores = 2
  )
  expect_true(is.list(gp))
  expect_named(gp, c("results_data", "results_plot", "results_summary"))

  # Stagnation branch.
  expect_warning(
    RMSTpowerBoost::GAM.ss.boot(
      dat, "time", "status", "arm", strata_var = "region",
      target_power = 0.99, linear_terms = "x1",
      L = 2, n_sim = 1, alpha = 0.1, parallel.cores = 1,
      patience = 1, n_start = 4, n_step = 1, max_n_per_arm = 8
    ),
    "stagnation"
  )

  # max_n_per_arm terminal warning branch.
  expect_warning(
    RMSTpowerBoost::GAM.ss.boot(
      dat, "time", "status", "arm", strata_var = "region",
      target_power = 0.999, linear_terms = "x1",
      L = 2, n_sim = 1, alpha = 0.1, parallel.cores = 1,
      patience = 99, n_start = 4, n_step = 1, max_n_per_arm = 5
    ),
    "not achieved by max N"
  )

  expect_warning(
    RMSTpowerBoost::GAM.ss.boot(
      dat, "time", "status", "arm", strata_var = "region",
      target_power = 0.9, linear_terms = "x1",
      L = 2, n_sim = 1, alpha = 0.1, parallel.cores = 1,
      patience = 1, n_start = 1, n_step = 1, max_n_per_arm = 2
    ),
    "stagnation|Canceling all iterations|NaNs produced"
  )
  expect_error(
    RMSTpowerBoost::GAM.ss.boot(
      dat, "time", "status", "arm", strata_var = "region",
      target_power = 0.9, linear_terms = "x1",
      L = 2, n_sim = 1, alpha = 0.1, parallel.cores = 1,
      patience = 1, n_start = 0, n_step = 1, max_n_per_arm = 1
    ),
    "second argument must be a list|Canceling all iterations"
  )
})

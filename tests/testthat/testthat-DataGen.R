test_that("example datasets are present and sane", {
   # rmst_weibull_200
   expect_silent(data(rmst_weibull_200, package = "RMSTSS"))
   expect_true(exists("rmst_weibull_200"), info = "rmst_weibull_200 missing")
   d <- rmst_weibull_200
   expect_true(is.data.frame(d))
   expect_true(all(c("time","status") %in% names(d)))
   expect_true(all(d$time > 0))
   expect_true(all(d$status %in% c(0L,1L)))

   # rmst_pwexp_strat_300 (if shipped)
   suppressWarnings(try(data(rmst_pwexp_strat_300, package = "RMSTSS"), silent = TRUE))
   if (exists("rmst_pwexp_strat_300")) {
      d2 <- rmst_pwexp_strat_300
      expect_true(is.data.frame(d2))
      expect_true(all(c("time","status") %in% names(d2)))
      expect_true(all(d2$time > 0))
      expect_true(all(d2$status %in% c(0L,1L)))
   }

   # rmst_small_50 (if shipped)
   suppressWarnings(try(data(rmst_small_50, package = "RMSTSS"), silent = TRUE))
   if (exists("rmst_small_50")) {
      d3 <- rmst_small_50
      expect_true(is.data.frame(d3))
      expect_true(all(c("time","status") %in% names(d3)))
      expect_true(all(d3$time > 0))
      expect_true(all(d3$status %in% c(0L,1L)))
   }
})

test_that("gen_covariates returns expected types for labeled vs Bernoulli", {
   defs <- list(
      list(name="x", type="continuous", dist="normal", params=list(mean=0, sd=1)),
      list(name="g", type="categorical", dist="bernoulli", params=list(p=0.5)),
      list(name="cat", type="categorical", dist="categorical",
           params=list(prob=c(0.2,0.3,0.5), labels=c("A","B","C")))
   )
   X <- gen_covariates(100, list(defs = defs))
   expect_true(is.numeric(X$x))
   expect_true(all(X$g %in% c(0L,1L)))     # Bernoulli: keep 0/1 numeric
   expect_true(is.factor(X$cat))           # labeled categorical: factor
   expect_identical(levels(X$cat), c("A","B","C"))
})

test_that("simulate_from_recipe (AFT) returns valid data and attributes", {
   covs <- list(
      list(name="age", type="continuous", dist="normal", params=list(mean=62, sd=10),
           transform=c("center(60)","scale(10)")),
      list(name="sex", type="categorical", dist="bernoulli", params=list(p=0.45))
   )
   rec <- recipe_quick_aft(
      n = 300, tau = 24, model = "aft_weibull",
      baseline = list(shape = 1.3, scale = 12),
      treat_effect = -0.25, covariates = covs,
      target_censoring = 0.25, allocation = "1:1", seed = 42
   )
   # Ensure target is feasible: long admin time if supported by your schema
   try(rec$censoring$admin_time <- 1e6, silent = TRUE)

   dat <- simulate_from_recipe(rec, seed = 42)
   expect_equal(nrow(dat), 300)
   expect_true(all(dat$time > 0))
   expect_true(all(dat$status %in% c(0L,1L)))
   if ("arm" %in% names(dat)) expect_true(all(dat$arm %in% c(0L,1L)))
   # Attributes (if your simulator sets them)
   ac <- attr(dat, "achieved_censoring", exact = TRUE)
   tau <- attr(dat, "tau", exact = TRUE)
   expect_true(is.numeric(tau) || is.null(tau))
   expect_true(is.numeric(ac)  || is.null(ac))
})

test_that("generate_recipe_sets + load_recipe_sets round-trip", {
   skip_on_cran()
   skip_on_os("solaris")

   covs <- list(list(name="x", type="continuous", dist="normal", params=list(mean=0, sd=1)))
   base <- recipe_quick_aft(
      n = 60, tau = 12, model = "aft_lognormal",
      baseline = list(mu = 2.7, sigma = 0.6),
      treat_effect = -0.2, covariates = covs,
      target_censoring = 0.20, allocation = "1:1", seed = 123
   )
   out_dir <- file.path(tempdir(), "rmstss-sets-ci")
   unlink(out_dir, recursive = TRUE, force = TRUE)

   man <- generate_recipe_sets(
      base_recipe = base,
      vary = list(n = c(50, 60), "event_time.effects.treatment" = c(-0.1, -0.2)),
      out_dir = out_dir,
      formats = c("rds","csv","txt"),
      n_reps = 2,
      seed_base = 2025,
      filename_template = "n{n}_te{event_time.effects.treatment}_sc{scenario_id}_r{rep}"
   )
   expect_true(file.exists(file.path(out_dir, "manifest.rds")))
   m <- readRDS(file.path(out_dir, "manifest.rds"))
   expect_gt(nrow(m), 0)

   sets <- load_recipe_sets(file.path(out_dir, "manifest.rds"))
   expect_true(length(sets) == nrow(m))
   s1 <- sets[[1]]
   expect_true(is.data.frame(s1$data))
   expect_true(all(s1$data$time > 0))
   expect_true(all(s1$data$status %in% c(0L,1L)))
})

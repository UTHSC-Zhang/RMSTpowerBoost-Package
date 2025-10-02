# tests/testthat/test-all-data.R

context("Data generation & metadata (no L, PH naming)")

test_that("AFT lognormal simulation returns valid data and attributes", {
   set.seed(1001)
   covs <- list(
      list(name="age", type="continuous",  dist="normal",    params=list(mean=60, sd=10)),
      list(name="sex", type="categorical", dist="bernoulli", params=list(p=0.5)),
      list(name="x1",  type="continuous",  dist="normal",    params=list(mean=0, sd=1)),
      list(name="x2",  type="continuous",  dist="gamma",     params=list(shape=2, scale=1))
   )

   rec <- recipe_quick_aft(
      n = 200,
      model = "aft_lognormal",
      baseline = list(mu = 2.7, sigma = 0.6),
      treat_effect = -0.2,
      covariates = covs,
      target_censoring = 0.25,
      allocation = "1:1",
      seed = 42
   )
   dat <- simulate_from_recipe(rec, seed = 42)

   expect_true(is.data.frame(dat))
   expect_true(all(c("time","status") %in% names(dat)))
   expect_true(all(dat$time > 0))
   expect_true(all(dat$status %in% c(0L,1L)))
   expect_false(anyNA(dat$time))
   ac <- attr(dat, "achieved_censoring")
   expect_true(is.numeric(ac) && is.finite(ac) && ac >= 0 && ac <= 1)
})

test_that("PH Weibull with gamma frailty simulates positive times", {
   set.seed(2025)
   covs <- list(
      list(name="age",     type="continuous",  dist="normal",    params=list(mean=62, sd=10)),
      list(name="sex",     type="categorical", dist="bernoulli", params=list(p=0.45)),
      list(name="cluster", type="categorical", dist="categorical",
           params=list(prob=c(0.3,0.4,0.3), labels=c("c1","c2","c3")))
   )

   rec <- validate_recipe(list(
      n = 300, covariates = list(defs = covs),
      treatment = list(assignment="randomization", allocation="1:1"),
      event_time = list(
         model="ph_weibull",
         baseline=list(shape=1.2, scale=10),
         effects=list(intercept=0, treatment=-0.3, covariates=list(age=0.01)),
         frailty=list(type="gamma", var=0.5, group="cluster")
      ),
      censoring = list(mode="target_overall", target=0.25, admin_time=36)
   ))

   dat <- simulate_from_recipe(rec, seed = 31415)
   expect_true(nrow(dat) == 300)
   expect_true(all(is.finite(dat$time)) && all(dat$time > 0))
   expect_true(mean(dat$status == 1L) > 0)
})

test_that("generate_recipe_sets + load_recipe_sets round trip (no L)", {
   set.seed(777)
   covs <- list(
      list(name="x",   type="continuous",  dist="normal",    params=list(mean=0, sd=1)),
      list(name="sex", type="categorical", dist="bernoulli", params=list(p=0.5)),
      list(name="x2",  type="continuous",  dist="gamma",     params=list(shape=2, scale=1))
   )
   rec <- recipe_quick_aft(
      n = 80,
      model = "aft_lognormal",
      baseline = list(mu = 2.5, sigma = 0.5),
      treat_effect = -0.1,
      covariates = covs,
      target_censoring = 0.15,
      allocation = "1:1"
   )
   od <- file.path(tempdir(), "rmstss-rt-one")
   dir.create(od, showWarnings = FALSE, recursive = TRUE)
   man  <- generate_recipe_sets(rec, out_dir = od, formats = c("csv","rds"), n_reps = 2L, seed_base = 2025)
   sets <- load_recipe_sets(file.path(od, "manifest.rds"))

   expect_true(nrow(man) == length(sets))
   expect_true(all(vapply(sets, function(s) is.list(s$meta), logical(1))))
   expect_true(all(vapply(sets, function(s) is.data.frame(s$data), logical(1))))
   # achieved_censoring present
   expect_true(all(is.finite(vapply(sets, function(s) attr(s$data, "achieved_censoring"), numeric(1)))))
})

test_that("describe_generation returns readable components when exported", {
   # Only run if describe_generation is exported; otherwise skip gracefully
   if (!("describe_generation" %in% getNamespaceExports("RMSTpowerBoost"))) {
      testthat::skip("describe_generation() not exported in this build")
   }
   covs <- list(list(name="x", type="continuous", dist="normal", params=list(mean=0, sd=1)))
   rec  <- recipe_quick_aft(
      n = 50, model = "aft_lognormal",
      baseline = list(mu=2.3, sigma=0.5),
      treat_effect = -0.15, covariates = covs,
      target_censoring = 0.2, allocation = "1:1"
   )
   od  <- file.path(tempdir(), "rmstss-desc-one")
   dir.create(od, showWarnings = FALSE, recursive = TRUE)
   generate_recipe_sets(rec, out_dir = od, formats = c("rds"))
   sets <- load_recipe_sets(file.path(od, "manifest.rds"))
   spec <- describe_generation(sets[[1]])
   expect_true(is.data.frame(spec$effects))
   expect_true(is.data.frame(spec$covariates))
})

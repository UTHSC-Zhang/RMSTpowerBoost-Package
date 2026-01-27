test_that("load_recipe_sets handles rds, rdata, and txt inputs", {
  td <- file.path(tempdir(), "rmst-load-branches")
  dir.create(td, showWarnings = FALSE, recursive = TRUE)

  df <- data.frame(
    time = c(1, 2, 3, 4),
    status = c(1, 0, 1, 0),
    arm = c(0, 1, 0, 1)
  )

  rds_path <- file.path(td, "set1.rds")
  saveRDS(df, rds_path)

  dat_obj <- df
  rdata_path <- file.path(td, "set2.RData")
  save(dat_obj, file = rdata_path)

  txt_path <- file.path(td, "set3.txt")
  utils::write.table(df, txt_path, sep = "\t", row.names = FALSE)

  man <- data.frame(
    scenario_id = c(1, 1, 1),
    rep = c(1, 2, 3),
    seed = c(NA_integer_, NA_integer_, NA_integer_),
    file_rds = c(rds_path, NA_character_, NA_character_),
    file_rdata = c(NA_character_, rdata_path, NA_character_),
    file_csv = c(NA_character_, NA_character_, NA_character_),
    file_txt = c(NA_character_, NA_character_, txt_path),
    stringsAsFactors = FALSE
  )

  saveRDS(man, file.path(td, "manifest.rds"))

  sets <- load_recipe_sets(file.path(td, "manifest.rds"))
  expect_true(length(sets) == 3)
  expect_true(all(vapply(sets, function(s) is.data.frame(s$data), logical(1))))
  expect_true(all(vapply(sets, function(s) is.list(s$meta), logical(1))))
  expect_true(all(is.finite(vapply(sets, function(s) attr(s$data, "achieved_censoring"), numeric(1)))))
})

test_that("rebuild_manifest reconstructs metadata from existing files", {
  set.seed(404)
  covs <- list(
    list(name = "x", type = "continuous", dist = "normal", params = list(mean = 0, sd = 1))
  )

  rec <- recipe_quick_aft(
    n = 30,
    model = "aft_lognormal",
    baseline = list(mu = 2.4, sigma = 0.6),
    treat_effect = -0.1,
    covariates = covs,
    target_censoring = 0.1,
    allocation = "1:1",
    seed = 101
  )

  dat <- simulate_from_recipe(rec, seed = 101)

  out_dir <- file.path(tempdir(), "rmst-rebuild-manifest")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  saveRDS(dat, file.path(out_dir, "sc1_r1.rds"))

  man <- rebuild_manifest(rec, vary = list(), out_dir = out_dir)
  expect_true(nrow(man) == 1)
  expect_true(file.exists(file.path(out_dir, "manifest.rds")))
  expect_true(is.list(man$meta[[1]]))
})

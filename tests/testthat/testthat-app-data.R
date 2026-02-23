test_that("app bundled datasets exist and are valid data.frames", {
  data_dir <- system.file("shiny_app", "data", package = "RMSTpowerBoost")
  expect_true(dir.exists(data_dir))

  rda_files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE, ignore.case = TRUE)
  rds_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
  expect_gt(length(rda_files) + length(rds_files), 0)

  loaded <- list()
  for (fp in rda_files) {
    e <- new.env(parent = emptyenv())
    load(fp, envir = e)
    for (nm in ls(e, all.names = TRUE)) {
      obj <- get(nm, envir = e, inherits = FALSE)
      if (is.data.frame(obj)) loaded[[nm]] <- obj
      if (is.list(obj) && !is.data.frame(obj)) {
        for (k in names(obj)) if (is.data.frame(obj[[k]])) loaded[[k]] <- obj[[k]]
      }
    }
  }
  for (fp in rds_files) {
    obj <- readRDS(fp)
    nm <- tools::file_path_sans_ext(basename(fp))
    if (is.data.frame(obj)) loaded[[nm]] <- obj
    if (is.list(obj) && !is.data.frame(obj)) {
      for (k in names(obj)) if (is.data.frame(obj[[k]])) loaded[[k]] <- obj[[k]]
    }
  }

  expect_gt(length(loaded), 0)
  for (nm in names(loaded)) {
    expect_true(is.data.frame(loaded[[nm]]))
    expect_true(all(c("time", "status", "arm") %in% names(loaded[[nm]])))
  }
})

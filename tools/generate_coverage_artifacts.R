#!/usr/bin/env Rscript

if (!requireNamespace("covr", quietly = TRUE)) {
  stop("Package 'covr' is required. Install it with install.packages('covr').", call. = FALSE)
}

cov <- covr::package_coverage(quiet = FALSE)
covr::to_cobertura(cov, filename = file.path("coverage", "cobertura.xml"))
# saveRDS(cov, file = file.path("coverage", "coverage.rds"))

message("Coverage artifacts written to ./coverage:")
message(" - coverage/cobertura.xml")
# message(" - coverage/coverage.rds")

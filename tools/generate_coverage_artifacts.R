#!/usr/bin/env Rscript

if (!requireNamespace("covr", quietly = TRUE)) {
  stop("Package 'covr' is required. Install it with install.packages('covr').", call. = FALSE)
}

# Run only the testthat suite, as CI does: tests/ also holds a manual
# validation script that needs truncnorm and devtools.
cov <- covr::package_coverage(
  type = "none",
  code = "testthat::test_package('RMSTpowerBoost')",
  quiet = FALSE
)
covr::to_cobertura(cov, filename = file.path("coverage", "cobertura.xml"))
# saveRDS(cov, file = file.path("coverage", "coverage.rds"))

message("Coverage artifacts written to ./coverage:")
message(" - coverage/cobertura.xml")
# message(" - coverage/coverage.rds")

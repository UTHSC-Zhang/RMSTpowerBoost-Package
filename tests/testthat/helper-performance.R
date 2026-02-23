is_covr_run <- function() {
  nzchar(Sys.getenv("R_COVR")) || tolower(Sys.getenv("COVR")) == "true"
}

covr_n <- function(covr_value, default_value) {
  if (is_covr_run()) covr_value else default_value
}

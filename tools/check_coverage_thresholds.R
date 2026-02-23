#!/usr/bin/env Rscript

if (!requireNamespace("covr", quietly = TRUE)) {
  stop("Package 'covr' is required.", call. = FALSE)
}

threshold <- 90
cov_file <- file.path("coverage", "coverage.rds")
if (!file.exists(cov_file)) {
  stop("Coverage artifact not found: ", cov_file, call. = FALSE)
}

cov <- readRDS(cov_file)
overall <- covr::percent_coverage(cov)
df <- as.data.frame(cov)

by_file <- aggregate(
  value ~ filename,
  data = df,
  FUN = function(x) 100 * mean(x > 0)
)
names(by_file)[2] <- "pct"
by_file <- by_file[grepl("^R/.+\\.R$", by_file$filename), , drop = FALSE]
by_file <- by_file[order(by_file$pct), , drop = FALSE]

cat(sprintf("Overall coverage: %.2f%%\n", overall))
cat("Per-file coverage (R/*.R):\n")
print(by_file, row.names = FALSE)

low <- by_file[by_file$pct <= threshold, , drop = FALSE]
if (overall <= threshold || nrow(low) > 0) {
  cat("\nCoverage gate failed.\n")
  if (overall <= threshold) {
    cat(sprintf("- Overall coverage %.2f%% is not > %d%%.\n", overall, threshold))
  }
  if (nrow(low) > 0) {
    cat("- Files not > 90%:\n")
    for (i in seq_len(nrow(low))) {
      cat(sprintf("  * %s: %.2f%%\n", low$filename[i], low$pct[i]))
    }
  }
  quit(status = 1)
}

cat(sprintf("\nCoverage gate passed: overall > %d%% and all R/*.R > %d%%.\n", threshold, threshold))

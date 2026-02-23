#!/usr/bin/env Rscript

app_path <- file.path("inst", "shiny_app", "app.R")
if (!file.exists(app_path)) {
  stop("Missing synced app file: ", app_path, call. = FALSE)
}

txt <- readLines(app_path, warn = FALSE, encoding = "UTF-8")
bad_patterns <- c(
  "\\bsource\\s*\\(",
  "\\.app\\s*\\(",
  "list\\.files\\(app_r_dir",
  "Run the app from the RMSTpowerBoost-App repo root"
)

hits <- lapply(bad_patterns, function(p) grep(p, txt, perl = TRUE))
names(hits) <- bad_patterns
has_bad <- any(vapply(hits, length, integer(1)) > 0)

if (has_bad) {
  message("Synced app validation failed. Forbidden patterns found:")
  for (p in names(hits)) {
    if (length(hits[[p]]) > 0) {
      message(" - ", p, " at lines: ", paste(hits[[p]], collapse = ", "))
    }
  }
  quit(status = 1)
}

message("Synced app validation passed: no forbidden app-repo-only patterns detected.")

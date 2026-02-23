#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
app_repo <- if (length(args) >= 1) args[[1]] else "../RMSTpowerBoost-App"

full_args <- commandArgs(FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg) == 1) {
  sub("^--file=", "", file_arg)
} else if ("-f" %in% full_args) {
  full_args[match("-f", full_args) + 1]
} else {
  file.path("tools", "sync_app_from_repo.R")
}
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

app_repo <- normalizePath(app_repo, winslash = "/", mustWork = FALSE)
if (!dir.exists(app_repo)) {
  stop("App repo not found: ", app_repo, call. = FALSE)
}

copy_file <- function(from, to) {
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(from, to, overwrite = TRUE)
  if (!ok) stop("Failed to copy: ", from, " -> ", to, call. = FALSE)
}

transform_app_text <- function(lines) {
  out <- lines

  # Ensure RMSTpowerBoost is included in packages vector
  has_pkg <- any(grepl("\"RMSTpowerBoost\"", out, fixed = TRUE))
  if (!has_pkg) {
    pkg_line <- grep("^\\s*\"kableExtra\".*\"tibble\"\\s*$", out)
    if (length(pkg_line) == 1) {
      out[pkg_line] <- sub(
        "\"tibble\"",
        "\"tibble\", \"RMSTpowerBoost\"",
        out[pkg_line],
        fixed = TRUE
      )
    }
  }

  # Replace app-repo sourcing block with package-require block
  start_idx <- grep("^# ------------------ Source your R/ scripts", out)
  if (length(start_idx) == 1) {
    end_idx <- grep("^`%\\|\\|%` <- function", out)
    if (length(end_idx) == 1 && end_idx > start_idx) {
      replacement <- c(
        "# ------------------ Use package functions directly ------------------",
        "if (!requireNamespace(\"RMSTpowerBoost\", quietly = TRUE)) {",
        "  stop(\"Package 'RMSTpowerBoost' is required for analysis methods.\", call. = FALSE)",
        "}"
      )
      out <- c(out[seq_len(start_idx - 1)], replacement, out[end_idx:length(out)])
    }
  }

  replacements <- c(
    "DC.power.analytical.app(" = "RMSTpowerBoost::DC.power.analytical(",
    "DC.ss.analytical.app(" = "RMSTpowerBoost::DC.ss.analytical(",
    "linear.power.analytical.app(" = "RMSTpowerBoost::linear.power.analytical(",
    "linear.ss.analytical.app(" = "RMSTpowerBoost::linear.ss.analytical(",
    "additive.power.analytical.app(" = "RMSTpowerBoost::additive.power.analytical(",
    "additive.ss.analytical.app(" = "RMSTpowerBoost::additive.ss.analytical(",
    "MS.power.analytical.app(" = "RMSTpowerBoost::MS.power.analytical(",
    "MS.ss.analytical.app(" = "RMSTpowerBoost::MS.ss.analytical(",
    "additive.power.boot.app(" = "RMSTpowerBoost::GAM.power.boot(",
    "additive.ss.boot.app(" = "RMSTpowerBoost::GAM.ss.boot(",
    "linear.power.boot.app(" = "RMSTpowerBoost::linear.power.boot(",
    "linear.ss.boot.app(" = "RMSTpowerBoost::linear.ss.boot(",
    "MS.power.boot.app(" = "RMSTpowerBoost::MS.power.boot(",
    "MS.ss.boot.app(" = "RMSTpowerBoost::MS.ss.boot("
  )

  for (k in names(replacements)) {
    out <- gsub(k, replacements[[k]], out, fixed = TRUE)
  }

  # Remove arguments that only existed in app-wrapper APIs
  out <- out[!grepl("dep_cens_status_var\\s*=", out)]
  out <- out[!grepl("^\\s*parallel\\.cores\\s*=", out)]

  # Remove app-repo-only runtime guard that enforces app repo root.
  guard_start <- grep("^is_app_repo <- ", out)
  if (length(guard_start) >= 1) {
    i <- guard_start[1]
    j <- i
    depth <- 0L
    started_if <- FALSE
    while (j <= length(out)) {
      line <- out[j]
      if (grepl("^if \\(interactive\\(\\) && !is_app_repo\\) \\{", line)) {
        started_if <- TRUE
      }
      if (started_if) {
        depth <- depth + lengths(regmatches(line, gregexpr("\\{", line, perl = TRUE)))
        depth <- depth - lengths(regmatches(line, gregexpr("\\}", line, perl = TRUE)))
        if (depth <= 0L && grepl("^\\}", line)) break
      }
      j <- j + 1L
    }
    if (j <= length(out)) {
      prefix <- if (i > 1L) out[1:(i - 1L)] else character(0)
      suffix <- if (j < length(out)) out[(j + 1L):length(out)] else character(0)
      out <- c(prefix, suffix)
    }
  }

  out
}

src_app <- file.path(app_repo, "app.R")
src_css <- file.path(app_repo, "www", "custom.css")
src_html <- file.path(app_repo, "www", "pipeline.html")

if (!file.exists(src_app)) stop("Missing app source: ", src_app, call. = FALSE)
if (!file.exists(src_css)) stop("Missing css source: ", src_css, call. = FALSE)
if (!file.exists(src_html)) stop("Missing html source: ", src_html, call. = FALSE)

dst_app <- file.path("inst", "shiny_app", "app.R")
dst_css <- file.path("inst", "shiny_app", "www", "custom.css")
dst_html <- file.path("inst", "shiny_app", "www", "pipeline.html")

dir.create(file.path("inst", "shiny_app", "www"), recursive = TRUE, showWarnings = FALSE)

app_lines <- readLines(src_app, warn = FALSE, encoding = "UTF-8")
app_lines <- transform_app_text(app_lines)
writeLines(app_lines, dst_app, useBytes = TRUE)
copy_file(src_css, dst_css)
copy_file(src_html, dst_html)

# Copy one app-alignment test (recipe simulator) and retarget namespace.
src_test <- file.path(app_repo, "tests", "testthat", "test-recipe-sim.R")
dst_test <- file.path("tests", "testthat", "testthat-sync-recipe-sim.R")
if (file.exists(src_test)) {
  tlines <- readLines(src_test, warn = FALSE, encoding = "UTF-8")
  tlines <- gsub("RMSTpowerBoostApp:::", "RMSTpowerBoost:::", tlines, fixed = TRUE)
  writeLines(tlines, dst_test, useBytes = TRUE)
}

message("Synced app artifacts from: ", app_repo)
message(" - ", dst_app)
message(" - ", dst_css)
message(" - ", dst_html)
if (file.exists(dst_test)) message(" - ", dst_test)

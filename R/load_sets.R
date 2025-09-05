
# R/load_sets.R

#' Load datasets from a recipe-sets manifest
#'
#' Reads a \code{manifest.rds} created by \code{generate_recipe_sets()}, loads one
#' dataset per row (preferring \code{rds} → \code{rdata} → \code{csv} → \code{txt}),
#' restores attributes \code{"tau"} and \code{"achieved_censoring"}, and returns a
#' named list of \code{list(data = <data.frame>, meta = <list>)}.
#' @param manifest Path to \code{manifest.rds} or a manifest \code{data.frame}.
#' @param prefer_formats Format preference order. Default \code{c("rds","rdata","csv","txt")}.
#' @param name_col Optional manifest column to use as list names; otherwise \code{key_template} is used.
#' @param key_template Template for names when \code{name_col} is \code{NULL}. Tokens: \code{"{scenario_id}"}, \code{"{rep}"}.
#' @param stringsAsFactors Passed to CSV/TXT readers.
#' @param skip_missing If TRUE, skip rows with missing files; else error.
#' @return A named list where each element is \code{list(data=..., meta=...)}.
#' @export
load_recipe_sets <- function(manifest,
                             prefer_formats = c("rds","rdata","csv","txt"),
                             name_col = NULL,
                             key_template = "sc{scenario_id}_r{rep}",
                             stringsAsFactors = FALSE,
                             skip_missing = TRUE) {
  if (is.character(manifest) && length(manifest) == 1L) {
    mf <- readRDS(manifest)
  } else if (is.data.frame(manifest)) {
    mf <- manifest
  } else stop("`manifest` must be a path to manifest.rds or a data.frame")

  req <- c("scenario_id","rep","tau","achieved_censoring","n","file_txt","file_csv","file_rds","file_rdata")
  missing_cols <- setdiff(req, names(mf))
  if (length(missing_cols)) stop("Manifest missing columns: ", paste(missing_cols, collapse=", "))

  prefer_formats <- tolower(prefer_formats)
  allowed <- c("rds","rdata","csv","txt")
  if (!all(prefer_formats %in% allowed))
    stop("prefer_formats must be subset of: ", paste(allowed, collapse=", "))

  choose_path <- function(row) {
    for (fmt in prefer_formats) {
      col <- switch(fmt, rds="file_rds", rdata="file_rdata", csv="file_csv", txt="file_txt")
      p <- row[[col]]
      if (!is.na(p) && nzchar(p) && file.exists(p)) return(list(path=p, fmt=fmt))
    }
    NULL
  }
  fill_key <- function(row) {
    out <- key_template
    out <- gsub("{scenario_id}", as.character(row$scenario_id), out, fixed = TRUE)
    out <- gsub("{rep}",         as.character(row$rep),         out, fixed = TRUE)
    out
  }
  load_one <- function(path, fmt, row) {
    meta <- list(seed = if ("seed" %in% names(row)) row$seed else NA_integer_,
                 tau = row$tau, achieved_censoring = row$achieved_censoring,
                 scenario_id = row$scenario_id, rep = row$rep, n = row$n,
                 file = path, format = fmt)
    pcols <- grep("^p__", names(row), value = TRUE)
    if (length(pcols)) {
      params <- as.list(row[pcols]); names(params) <- sub("^p__", "", names(params))
      meta$params <- params
    } else meta$params <- list()

    if (fmt == "rds") {
      dat <- readRDS(path)
      attr(dat, "tau") <- meta$tau; attr(dat, "achieved_censoring") <- meta$achieved_censoring
      return(list(data = dat, meta = meta))
    }
    if (fmt == "rdata") {
      env <- new.env(parent = emptyenv()); loaded <- load(path, envir = env)
      dat <- if (exists("dat_obj", envir = env, inherits = FALSE)) get("dat_obj", envir = env, inherits = FALSE) else {
        cand <- Filter(is.data.frame, mget(loaded, envir = env, inherits = FALSE))
        if (!length(cand)) stop("RData did not contain a data.frame: ", path)
        cand[[1]]
      }
      if (exists("meta", envir = env, inherits = FALSE)) {
        meta_in <- get("meta", envir = env, inherits = FALSE)
        for (nm in setdiff(names(meta_in), c("tau","achieved_censoring"))) meta[[nm]] <- meta_in[[nm]]
      }
      attr(dat, "tau") <- meta$tau; attr(dat, "achieved_censoring") <- meta$achieved_censoring
      return(list(data = dat, meta = meta))
    }
    if (fmt == "csv") {
      dat <- utils::read.csv(path, stringsAsFactors = stringsAsFactors)
      attr(dat, "tau") <- meta$tau; attr(dat, "achieved_censoring") <- meta$achieved_censoring
      return(list(data = dat, meta = meta))
    }
    if (fmt == "txt") {
      dat <- utils::read.table(path, header = TRUE, sep = "\t",
                               stringsAsFactors = stringsAsFactors, check.names = FALSE)
      attr(dat, "tau") <- meta$tau; attr(dat, "achieved_censoring") <- meta$achieved_censoring
      return(list(data = dat, meta = meta))
    }
    stop("Unhandled format: ", fmt)
  }

  out <- list()
  for (i in seq_len(nrow(mf))) {
    row <- mf[i, , drop = FALSE]
    sel <- choose_path(row)
    if (is.null(sel)) {
      msg <- sprintf("No readable file for scenario %s rep %s", row$scenario_id, row$rep)
      if (isTRUE(skip_missing)) { warning(msg); next } else stop(msg)
    }
    key <- if (!is.null(name_col) && name_col %in% names(row)) as.character(row[[name_col]]) else fill_key(row)
    out[[key]] <- load_one(sel$path, sel$fmt, row)
  }
  out
}

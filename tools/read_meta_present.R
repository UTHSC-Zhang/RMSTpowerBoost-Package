# tools/read_meta_present.R
# Read checks/manifest.rds and present a compact summary.

#' Present dataset metadata from checks/manifest.rds
#'
#' Reads the manifest produced by generate_data_and_meta.R and returns a tidy
#' summary data.frame. Optionally prints a compact view.
#'
#' @param manifest_path Path to 'checks/manifest.rds' (default).
#' @param print Logical; print a compact table to console (default TRUE).
#' @return A data.frame with key, n, L, model, event_rate, achieved_censoring, allocation, file paths.
#' @examples
#' \dontrun{
#' present_manifest()                        # default path
#' df <- present_manifest(print = FALSE)     # programmatic use
#' }
present_manifest <- function(manifest_path = file.path("checks", "manifest.rds"),
                             print = TRUE) {
   `%||%` <- function(x, y) if (is.null(x)) y else x

   if (!file.exists(manifest_path)) stop("Manifest not found: ", manifest_path)
   man <- readRDS(manifest_path)
   if (!nrow(man)) {
      if (print) cat("Manifest is empty:", normalizePath(manifest_path), "\n")
      return(invisible(man))
   }

   # Pull a few nested bits from meta
   unpack <- function(m) {
      c(
         model = m$model %||% NA_character_,
         L     = as.numeric(m$L %||% NA_real_),
         alloc = (m$treatment$allocation %||% NA_character_)
      )
   }
   add_cols <- do.call(rbind, lapply(man$meta, unpack))
   add_cols <- as.data.frame(add_cols, stringsAsFactors = FALSE)
   suppressWarnings(add_cols$L <- as.numeric(add_cols$L))

   out <- data.frame(
      key     = man$key,
      n       = man$n,
      L       = ifelse(is.na(man$L), add_cols$L, man$L),
      model   = add_cols$model,
      event_rate = round(man$event_rate, 3),
      achieved_censoring = round(man$achieved_censoring, 3),
      allocation = add_cols$alloc,
      file_csv = man$file_csv,
      file_rda = man$file_rda,
      stringsAsFactors = FALSE
   )

   if (print) {
      cat("\nDatasets summary (", normalizePath(manifest_path), "):\n\n", sep = "")
      show <- out
      show$file_csv <- basename(show$file_csv)
      show$file_rda <- basename(show$file_rda)
      print(show, row.names = FALSE, right = FALSE)
      cat("\nTip: each row’s full metadata is in the 'meta' list-column of the manifest.\n")
   }
   invisible(out)
}
present_manifest()

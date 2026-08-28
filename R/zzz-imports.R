# R/zzz-imports.R
# Centralized roxygen imports to satisfy R CMD check for non-exported helpers.

#' @importFrom stats uniroot setNames
#' @importFrom utils head tail
NULL

# @noRd
.rmst_verbose_message <- function(verbose, ...) {
   if (isTRUE(verbose)) {
      message(..., appendLF = TRUE)
   }
}

# Null-coalescing operator shared across the package (single definition).
# @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

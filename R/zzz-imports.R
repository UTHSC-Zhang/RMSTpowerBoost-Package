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

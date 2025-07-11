#' Launch the RMSTdesign Shiny Application
#'
#' @description
#' A function to easily launch the interactive Shiny application bundled with
#' the package.
#'
#' @export
#' @importFrom shiny runApp
#'
#' @examples
#' \dontrun{
#'   RMSTSS::run_app()
#' }
run_app <- function() {
   r_files <- list.files(path = "R/", pattern = "\\.R$", full.names = TRUE)
   if (length(r_files) == 0) {
      stop("No R files found in the 'R/' directory of the RMSTSS package.", call. = FALSE)
   }
   sapply(r_files, source)
   # Check if the shiny package is installed
   if (!requireNamespace("shiny", quietly = TRUE)) {
      stop("The 'shiny' package is required to run the RMSTSS application. Please install it using `install.packages('shiny')`.", call. = FALSE)
   }
   # Get the path to the Shiny app directorr
   app_dir <- system.file("shiny_app", package = "RMSTSS")

   # Check if the directory exists
   if (app_dir == "") {
      stop(
         "Could not find the app directory. ",
         "Try re-installing the `RMSTdesign-Package`.",
         call. = FALSE
      )
   }

   # Launch the app
   shiny::runApp(app_dir, display.mode = "normal")
}


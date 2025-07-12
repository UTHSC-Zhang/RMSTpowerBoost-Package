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
   # First, load the package itself to make all its functions available
   # to the Shiny app's environment. This is a crucial step.
   library(RMSTSS)

   # Get the path to the Shiny app directory within the installed package
   app_dir <- system.file("shiny_app", package = "RMSTSS")

   # Check if the directory exists
   if (app_dir == "") {
      stop(
         "Could not find the app directory. ",
         "Try re-installing the `RMSTSS` package.",
         call. = FALSE
      )
   }

   # Launch the app
   shiny::runApp(app_dir, display.mode = "normal")
}

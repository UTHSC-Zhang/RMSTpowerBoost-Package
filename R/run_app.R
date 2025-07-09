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
#'   run_app()
#' }
run_app <- function() {
   # Find the app's directory within the installed package
   app_dir <- system.file("shiny_app", package = "RMSTSS")

   # Check if the directory exists
   if (app_dir == "") {
      stop(
         "Could not find the app directory. ",
         "Try re-installing the `RMSTSS-Package`.",
         call. = FALSE
      )
   }

   # Launch the app
   shiny::runApp(app_dir, display.mode = "normal")
}

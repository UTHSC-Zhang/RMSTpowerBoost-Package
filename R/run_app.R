#' Launch the RMSTdesign Shiny Application
#'
#' @description
#' A function to easily launch the interactive Shiny application bundled with
#' the package.
#'
#' @export
#' @importFrom shiny runApp
#' @importFrom bslib bs_themer
#' @importFrom DT datatable renderDataTable
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_light ylim
#' @importFrom kableExtra kbl kable_styling
#' @importFrom magrittr %>%
#' @importFrom plotly ggplotly renderPlotly
#' @importFrom rmarkdown render
#' @importFrom shiny reactive req showNotification renderUI tagList fluidRow
#' @importFrom shiny column selectInput selectizeInput radioButtons textInput
#' @importFrom shiny sliderInput wellPanel h4 numericInput observe observeEvent
#' @importFrom shiny updateSelectInput reactiveVal validate need withProgress
#' @importFrom shiny setProgress renderText downloadButton downloadHandler
#' @importFrom shiny removeNotification hr p bindCache bindEvent
#' @importFrom shinyjs toggle toggleState reset
#' @importFrom stats as.formula complete.cases pchisq sd na.omit
#' @importFrom survival survfit Surv survdiff
#' @importFrom survminer ggsurvplot
#' @examples
#' \dontrun{
#'   RMSTSS::run_app()
#' }
run_app <- function() {
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

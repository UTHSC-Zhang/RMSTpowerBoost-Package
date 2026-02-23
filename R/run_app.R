# Internal app dependency package list.
# @noRd
.app_dependency_packages <- function() {
   c(
      "shiny", "shinyjs", "bslib", "DT", "plotly",
      "kableExtra", "survminer", "rmarkdown", "mice"
   )
}

# @noRd
.get_missing_app_dependencies <- function() {
   deps <- .app_dependency_packages()
   deps[!vapply(deps, requireNamespace, logical(1), quietly = TRUE)]
}

# @noRd
.install_command_for <- function(pkgs) {
   paste0("install.packages(c(", paste(shQuote(pkgs), collapse = ", "), "))")
}

# @noRd
.prompt_install_app_dependencies <- function(missing) {
   if (!interactive()) return(FALSE)
   choice <- utils::menu(
      choices = c("Install missing app dependencies", "Cancel launch"),
      title = paste0(
         "run_app() requires additional packages:\n",
         paste(missing, collapse = ", "),
         "\nInstall now?"
      )
   )
   identical(choice, 1L)
}

# @noRd
.install_app_dependencies <- function(missing, repos) {
   utils::install.packages(missing, repos = repos)
}

# @noRd
.ensure_app_dependencies <- function(install_missing = TRUE, repos = getOption("repos")) {
   missing <- .get_missing_app_dependencies()
   if (!length(missing)) return(invisible(TRUE))

   install_cmd <- .install_command_for(missing)
   if (!isTRUE(install_missing)) {
      stop(
         "Missing app dependencies: ", paste(missing, collapse = ", "), ". ",
         "Install them with ", install_cmd, ".",
         call. = FALSE
      )
   }

   if (!.prompt_install_app_dependencies(missing)) {
      stop(
         "App launch canceled because dependencies are missing: ",
         paste(missing, collapse = ", "), ". ",
         "Install them with ", install_cmd, ".",
         call. = FALSE
      )
   }

   .install_app_dependencies(missing, repos = repos)
   still_missing <- .get_missing_app_dependencies()
   if (length(still_missing)) {
      stop(
         "App dependencies are still missing after install: ",
         paste(still_missing, collapse = ", "), ". ",
         "Try running ", .install_command_for(still_missing), ".",
         call. = FALSE
      )
   }

   invisible(TRUE)
}

#' Launch the RMSTdesign Shiny Application
#'
#' @description
#' A helper to launch the interactive Shiny application bundled with the package.
#' App-specific dependencies are installed lazily when needed.
#'
#' @param install_missing Logical; if `TRUE`, prompt to install missing app
#'   dependencies.
#' @param repos CRAN mirror(s) passed to `utils::install.packages()` when
#'   installing missing app dependencies.
#'
#' @return Invisible return value from `shiny::runApp()`.
#' @export
#'
#' @examples
#' \dontrun{
#'   RMSTpowerBoost::run_app()
#' }
run_app <- function(install_missing = TRUE, repos = getOption("repos")) {
   .ensure_app_dependencies(install_missing = install_missing, repos = repos)

   app_dir <- Sys.getenv("RMSTPOWERBOOST_APP_DIR", unset = NA_character_)
   if (is.na(app_dir)) {
      app_dir <- system.file("shiny_app", package = "RMSTpowerBoost")
   }

   if (app_dir == "" || !dir.exists(app_dir)) {
      stop(
         "Could not find the app directory. ",
         "Try re-installing the `RMSTpowerBoost` package.",
         call. = FALSE
      )
   }

   shiny::runApp(app_dir, display.mode = "normal")
}

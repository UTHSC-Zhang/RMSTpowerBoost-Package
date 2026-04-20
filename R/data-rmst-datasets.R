# R/data-aft-and-ph-examples.R

#' Simulated dataset: AFT log-normal, L = 12, n = 150
#'
#' Columns:
#' \itemize{
#'   \item time: observed time
#'   \item status: event indicator (1 = event, 0 = censored)
#'   \item arm: treatment arm (0 = control, 1 = treatment)
#'   \item age: baseline age
#'   \item gender: factor with two levels
#' }
#' @format A data frame with 150 rows and 5 variables.
#' @docType data
#' @name aft_lognormal_L12_n150
#' @keywords datasets
#' @examples
#' data(aft_lognormal_L12_n150)
#' table(aft_lognormal_L12_n150$arm)
NULL

#' Simulated dataset: AFT Weibull, L = 24, n = 200
#'
#' Columns: time, status, arm, age, gender.
#' @format A data frame with 200 rows and 5 variables.
#' @docType data
#' @name aft_weibull_L24_n200
#' @keywords datasets
#' @examples
#' data(aft_weibull_L24_n200)
#' mean(aft_weibull_L24_n200$status == 0)
NULL

#' Simulated dataset: PH piecewise-exponential, L = 18, n = 250
#'
#' Columns: time, status, arm, age, gender.
#' @format A data frame with 250 rows and 5 variables.
#' @docType data
#' @name ph_pwexp_L18_n250
#' @keywords datasets
#' @examples
#' data(ph_pwexp_L18_n250)
#' summary(ph_pwexp_L18_n250$time)
NULL

#' Simulated dataset: PH Weibull, L = 24, n = 300
#'
#' Columns: time, status, arm, age, gender.
#' @format A data frame with 300 rows and 5 variables.
#' @docType data
#' @name ph_weibull_L24_n300
#' @keywords datasets
#' @examples
#' data(ph_weibull_L24_n300)
#' prop.table(table(ph_weibull_L24_n300$status))
NULL


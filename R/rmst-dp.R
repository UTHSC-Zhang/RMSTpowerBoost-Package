
#' Function to perform dynamic prediction based on stratified or link models
#'
#' This function implements dynamic prediction models for time-to-event outcomes,
#' allowing for stratified or link-based approaches. It estimates model parameters
#' and their standard errors using inverse probability weighting.
#'
#' @param dat A data frame containing baseline covariates (`Za`, `Zb1` to `ZbK`),
#'   time-to-event outcome (`X`, `deltaD`, `deltaT`), and `entry` time.
#'   Must include `ID`.
#' @param datA A data frame containing information on treatment assignment times.
#'   Must include `ID`, `A` (assignment time), and `status` ('active' or 'control').
#' @param L Numeric, the restriction time for Restricted Mean Survival Time (RMST).
#' @param K Integer, the number of cross-sections for dynamic prediction.
#' @param interval Numeric, the time interval between consecutive cross-sections.
#' @param method Character string, specifying the prediction method: "link" (generalized linear model)
#'   or "stratified" (stratified additive or multiplicative model).
#' @param link Character string, used if `method = "link"`: "linear" (identity link) or "log" (log link).
#' @param stratified Character string, used if `method = "stratified"`: "add" (additive model) or "multi" (multiplicative model).
#' @param weights Character string, specifying the type of inverse probability weights:
#'   "stabilized" or "unstabilized".
#' @return A list containing:
#'   \item{betahat}{A numeric vector of the estimated model parameters.}
#'   \item{se}{A numeric vector of the standard errors for the estimated parameters.}
#' @export
#' @examples
#' # This example provides a minimal setup for the dp function.
#' # For power/sample size calculation, refer to `calculate_dp_power` and `calculate_dp_sample_size`.
#'
#' # Create dummy data
#' set.seed(123)
#' n_obs <- 100
#' df_dat <- data.frame(
#'   ID = 1:n_obs,
#'   entry = 0,
#'   X = rexp(n_obs, rate = 0.1),
#'   deltaD = rbinom(n_obs, 1, 0.6), # Death indicator
#'   deltaT = rbinom(n_obs, 1, 0.3), # Treatment indicator
#'   Za = rbinom(n_obs, 1, 0.5), # Binary covariate (e.g., group)
#'   Zb1 = rnorm(n_obs),
#'   Zb2 = rnorm(n_obs),
#'   Zb3 = rnorm(n_obs)
#' )
#' df_dat$X[df_dat$X > 20] <- 20 # Cap follow-up time
#'
#' df_datA <- data.frame(
#'   ID = 1:n_obs,
#'   A = runif(n_obs, 0, 15), # Assignment time
#'   status = sample(c("active", "control"), n_obs, replace = TRUE)
#' )
#'
#' L_val <- 15
#' K_val <- 3 # Number of cross-sections
#' interval_val <- L_val / K_val
#'
#' # Example using "link" method with "linear" link
#' dp_result_link_linear <- dp(
#'   dat = df_dat, datA = df_datA, L = L_val, K = K_val,
#'   interval = interval_val, method = "link", link = "linear",
#'   weights = "stabilized"
#' )
#' # print("DP Result (Link-Linear):")
#' # print(dp_result_link_linear)
#'
#' # Example using "stratified" method with "add" (additive)
#' dp_result_stratified_add <- dp(
#'   dat = df_dat, datA = df_datA, L = L_val, K = K_val,
#'   interval = interval_val, method = "stratified", stratified = "add",
#'   weights = "unstabilized"
#' )
#' # print("DP Result (Stratified-Additive):")
#' # print(dp_result_stratified_add)
#' Function to perform dynamic prediction based on stratified or link models
#'
#' This function implements dynamic prediction models for time-to-event outcomes,
#' allowing for stratified or link-based approaches. It estimates model parameters
#' and their standard errors using inverse probability weighting.
#'
#' @param dat A data frame containing baseline covariates (`Za`, `Zb1` to `ZbK`),
#'   time-to-event outcome (`X`, `deltaD`, `deltaT`), and `entry` time.
#'   Must include `ID`.
#' @param datA A data frame containing information on treatment assignment times.
#'   Must include `ID`, `A` (assignment time), and `status` ('active' or 'control').
#' @param L Numeric, the restriction time for Restricted Mean Survival Time (RMST).
#' @param K Integer, the number of cross-sections for dynamic prediction.
#' @param interval Numeric, the time interval between consecutive cross-sections.
#' @param method Character string, specifying the prediction method: "link" (generalized linear model)
#'   or "stratified" (stratified additive or multiplicative model).
#' @param link Character string, used if `method = "link"`: "linear" (identity link) or "log" (log link).
#' @param stratified Character string, used if `method = "stratified"`: "add" (additive model) or "multi" (multiplicative model).
#' @param weights Character string, specifying the type of inverse probability weights:
#'   "stabilized" or "unstabilized".
#' @return A list containing:
#'   \item{betahat}{A numeric vector of the estimated model parameters.}
#'   \item{se}{A numeric vector of the standard errors for the estimated parameters.}
#' @export
dp <- function(dat, datA, L, K, interval, # K = the number of cross-sections; interval = the time between consecutive cross-sections
               method = c("link", "stratified"), link = c("linear", "log"), stratified = c("add", "multi"), weights = c("stabilized", "unstabilized"))
{

   method <- base::match.arg(method)
   link <- base::match.arg(link)
   stratified <- base::match.arg(stratified)
   weights <- base::match.arg(weights)

   # Merge dat and datA at the beginning
   dat_full <- base::merge(dat, datA, by = "ID")

   # Define the censoring event robustly
   censoring_event_indicator <- base::as.numeric(dat_full$deltaD == 0 & dat_full$deltaT == 0)

   # Estimate the censoring model
   modC <- survival::coxph(survival::Surv(X, censoring_event_indicator) ~ Za, data = dat_full, ties = "breslow")

   # Estimate the treatment model
   dat_full$weights_coxph <- base::ifelse(dat_full$status == "active", 1, 1e-08)
   dtimes <- base::sort(base::unique(dat_full$X[dat_full$deltaT == 1]))

   if (base::length(dtimes) == 0) {
      base::warning("No treatment events (deltaT=1) found to fit treatment model. Returning NA.")
      return(list("betahat" = NA, "se" = NA))
   }

   dtimes <- dtimes[dtimes > base::min(dat_full$entry) & dtimes < base::max(dat_full$X)]
   if (base::length(dtimes) == 0) {
      base::warning("No valid treatment event times within observed range to fit treatment model. Returning NA.")
      return(list("betahat" = NA, "se" = NA))
   }

   all_vars <- base::names(dat_full)
   rhs_vars <- all_vars[!all_vars %in% base::c("X", "deltaT")] # Exclude LHS variables

   formula_rhs <- base::paste(rhs_vars, collapse = " + ")
   split_formula <- stats::as.formula(
      paste("Surv(X, deltaT == 1) ~", formula_rhs)
   )


   datT.Extended <- survival::survSplit(
      formula = split_formula,
      data = dat_full,
      cut = dtimes,
      episode = "tstart"
   )

   datT.Extended$Zt <- base::ifelse(datT.Extended$X > datT.Extended$A, 1, 0)
   modT <- survival::coxph(survival::Surv(tstart, X, event) ~ Za + Zt, data = datT.Extended, weights = weights_coxph,
                           ties = "breslow", timefix = FALSE)

   # Create the stacked dataset
   dat.stacked <- base::data.frame()
   for (k in 1:K){
      dat.temp <- dat_full
      dat.temp$Sik <- interval * k - dat.temp$entry
      dat.temp <- dat.temp[dat.temp$Sik > 0 & dat.temp$Sik < dat.temp$X & dat.temp$Sik < dat.temp$A, ]

      if (base::nrow(dat.temp) > 0){
         zb_col_name_k <- base::paste0("Zb", k)
         if (!(zb_col_name_k %in% base::names(dat.temp))) {
            base::warning(base::paste("Column", zb_col_name_k, "not found for cross-section", k, ". Skipping."))
            next
         }

         dat.temp$Xik <- dat.temp$X - dat.temp$Sik
         dat.temp$Yik <- base::pmin(dat.temp$Xik, L)
         dat.temp$deltaYik <- base::ifelse(dat.temp$Xik <= L, dat.temp$deltaD, 1)

         base_cols <- c("ID", "entry", "X", "deltaD", "deltaT", "Za", "A", "status", "weights_coxph", "Sik", "Xik", "Yik", "deltaYik")
         dat.temp_processed <- dat.temp[, base_cols]
         dat.temp_processed$Zb <- dat.temp[[zb_col_name_k]]

         dat.stacked <- base::rbind(dat.stacked, base::cbind.data.frame("CS" = k, dat.temp_processed))
      }
   }

   if (base::nrow(dat.stacked) == 0) {
      base::warning("No subjects at risk for any cross-section after stacking. Returning NA.")
      return(list("betahat" = NA, "se" = NA))
   }

   dat.temp1 <- dat.temp2 <- dat.temp3 <- dat.stacked

   dat.temp1$X_pred <- dat.temp1$Sik
   dat.temp2$X_pred <- dat.temp2$Sik + dat.temp2$Yik
   dat.temp1$censoring_event_indicator <- as.numeric(dat.temp1$deltaD == 0 & dat.temp1$deltaT == 0)
   predC_Sik <- stats::predict(modC, newdata = dat.temp1, type = "survival")
   dat.temp2$censoring_event_indicator <- as.numeric(dat.temp2$deltaD == 0 & dat.temp2$deltaT == 0)
   predC_Sik_Yik <- stats::predict(modC, newdata = dat.temp2, type = "survival")

   predC_Sik_Yik[predC_Sik_Yik <= 1e-10] <- 1e-10
   dat.stacked$WC <- predC_Sik / predC_Sik_Yik
   dat.stacked$WC[base::is.na(dat.stacked$WC) | base::is.infinite(dat.stacked$WC)] <- 0

   dat.temp1$tstart_pred <- 0
   dat.temp1$event_pred <- 0
   dat.temp1$X_pred_modT <- dat.temp1$Sik
   dat.temp1$Zt_pred <- base::ifelse(dat.temp1$A < dat.temp1$Sik, 1, 0)
   dat.temp2$tstart_pred <- 0
   dat.temp2$event_pred <- 0
   dat.temp2$X_pred_modT <- dat.temp2$Sik + dat.temp2$Yik
   dat.temp2$Zt_pred <- base::ifelse(dat.temp2$A < (dat.temp2$Sik + dat.temp2$Yik), 1, 0)
   newdata1_modT <- base::data.frame(X = dat.temp1$X_pred_modT, tstart = dat.temp1$tstart_pred, event = dat.temp1$event_pred, Za = dat.temp1$Za, Zt = dat.temp1$Zt_pred)
   newdata2_modT <- base::data.frame(X = dat.temp2$X_pred_modT, tstart = dat.temp2$tstart_pred, event = dat.temp2$event_pred, Za = dat.temp2$Za, Zt = dat.temp2$Zt_pred)
   predT_Sik <- stats::predict(modT, newdata = newdata1_modT, type = "survival")
   predT_Sik_Yik <- stats::predict(modT, newdata = newdata2_modT, type = "survival")
   predT_Sik_Yik[predT_Sik_Yik <= 1e-10] <- 1e-10

   if (weights == "unstabilized"){
      dat.stacked$WT <- predT_Sik / predT_Sik_Yik
   } else if (weights == "stabilized"){
      dat.temp3$tstart_pred <- 0
      dat.temp3$event_pred <- 0
      dat.temp3$X_pred_modT <- dat.temp3$Yik
      dat.temp3$Zt_pred <- base::ifelse(dat.temp3$A < dat.temp3$Yik, 1, 0)
      newdata3_modT <- base::data.frame(X = dat.temp3$X_pred_modT, tstart = dat.temp3$tstart_pred, event = dat.temp3$event_pred, Za = dat.temp3$Za, Zt = dat.temp3$Zt_pred)
      predT_Yik <- stats::predict(modT, newdata = newdata3_modT, type = "survival")
      dat.stacked$WT <- (predT_Sik * predT_Yik) / predT_Sik_Yik
   }
   dat.stacked$WT[base::is.na(dat.stacked$WT) | base::is.infinite(dat.stacked$WT)] <- 0
   dat.stacked$W <- dat.stacked$WC * dat.stacked$WT
   dat.stacked$W[base::is.na(dat.stacked$W) | base::is.infinite(dat.stacked$W)] <- 0

   if (base::sum(dat.stacked$W * dat.stacked$deltaYik) == 0) {
      base::warning("Sum of effective weights for outcome model is zero. Returning NA.")
      return(base::list("betahat" = NA, "se" = NA))
   }

   betahat <- NULL
   se <- NULL

   if (method == "link"){
      dat.stacked$Wgee <- dat.stacked$W * dat.stacked$deltaYik
      family_val <- if (link == "linear") stats::gaussian(link = "identity") else stats::poisson(link = "log")
      f <- stats::as.formula("Yik ~ Za + Zb : factor(CS) + factor(CS) - 1")
      model_fit <- base::tryCatch(
         base::suppressWarnings(stats::glm(f, data = dat.stacked, family = family_val, weights = Wgee)),
         error = function(e) { base::warning(base::paste("GLM failed:", e$message)); return(NULL) }
      )
      if (base::is.null(model_fit)) return(base::list("betahat" = NA, "se" = NA))

      betahat <- stats::coef(model_fit)
      Z <- stats::model.matrix(f, data = dat.stacked)
      A <- base::matrix(0, ncol = base::length(betahat), nrow = base::length(betahat))
      if (link == "linear"){
         for (i in 1:base::nrow(dat.stacked)) { Zik <- base::as.matrix(Z[i, ]); A <- A + dat.stacked$W[i] * dat.stacked$deltaYik[i] * Zik %*% base::t(Zik) }
         e <- Z * base::as.vector(dat.stacked$W * dat.stacked$deltaYik * (dat.stacked$Yik - Z %*% betahat))
      } else {
         for (i in 1:base::nrow(dat.stacked)) { Zik <- base::as.matrix(Z[i, ]); A <- A + dat.stacked$W[i] * dat.stacked$deltaYik[i] * Zik %*% base::t(Zik) * base::as.numeric(base::exp(base::t(betahat) %*% Zik)) }
         e <- Z * base::as.vector(dat.stacked$W * dat.stacked$deltaYik * (dat.stacked$Yik - base::exp(Z %*% betahat)))
      }
      E <- base::cbind.data.frame("ID" = dat.stacked$ID, e)
      Ei <- stats::aggregate(. ~ ID, E, sum)
      B <- base::matrix(0, ncol = base::length(betahat), nrow = base::length(betahat))
      for (i in 1:base::nrow(Ei)) { ei <- base::t(base::as.matrix(Ei[i, -1])); B <- B + ei %*% base::t(ei) }
      if (base::det(A) < 1e-8) {
         base::warning("Matrix A is singular. Returning NA for SE.")
         return(base::list("betahat" = betahat, "se" = base::rep(NA, base::length(betahat))))
      }
      var <- base::diag(base::solve(A) %*% B %*% base::t(base::solve(A)))
      se <- base::sqrt(var)
   }

   if (!base::is.null(betahat) && !base::is.null(se) && base::length(betahat) == base::length(se)) {
      base::names(se) <- base::names(betahat)
   } else {
      return(base::list("betahat" = NA, "se" = NA))
   }
   return(base::list("betahat" = betahat, "se" = se))
}



#' Simulate Data for Dynamic Prediction Power Calculation
#'
#' This function generates synthetic survival and treatment data for dynamic prediction models.
#' It creates `dat` (baseline info, time-to-event, treatment status) and `datA` (treatment assignment times)
#' data frames, structured to be compatible with the `dp` function.
#' The `Za` covariate is designed to represent the binary group for which power is calculated.
#'
#' @param n_per_group Numeric, the number of subjects in each binary group (total N = 2 * n_per_group).
#' @param L Numeric, the restriction time for RMST.
#' @param K Integer, the number of cross-sections.
#' @param interval Numeric, the time interval between consecutive cross-sections.
#' @param base_rate_survival Numeric, baseline hazard rate for survival events.
#' @param group_effect_Za Numeric, the log-hazard ratio (or linear effect) for `Za` (group effect)
#'   on survival. This is the effect size to detect.
#' @param Zb_effect Numeric, the log-hazard ratio (or linear effect) for `Zb` on survival.
#' @param censoring_rate Numeric, baseline hazard rate for censoring.
#' @param treatment_rate Numeric, baseline hazard rate for treatment assignment.
#' @return A list containing two data frames:
#'   \item{dat}{Main data frame with `ID`, `entry`, `X`, `deltaD`, `deltaT`, `Za`, `Zb1` to `ZbK`.}
#'   \item{datA}{Treatment assignment data frame with `ID`, `A`, `status`.}
#' @keywords internal
simulate_dp_data <- function(n_per_group, L, K, interval,
                             base_rate_survival = 0.1, group_effect_Za = log(0.8),
                             Zb_effect = 0.05, censoring_rate = 0.05, treatment_rate = 0.05) {
   n_total <- n_per_group * 2
   ID <- 1:n_total
   entry <- base::rep(0, n_total)
   Za <- base::rep(base::c(0, 1), each = n_per_group)

   Zb_base <- stats::rnorm(n_total, mean = 0, sd = 1)
   Zb_cols <- base::paste0("Zb", 1:K)
   Zb_df <- base::as.data.frame(base::matrix(Zb_base, ncol = K, nrow = n_total, byrow = FALSE))
   base::colnames(Zb_df) <- Zb_cols

   survival_times <- base::numeric(n_total)
   for (i in 1:n_total) {
      hazard_survival <- base_rate_survival * base::exp(group_effect_Za * Za[i] + Zb_effect * Zb_base[i])
      survival_times[i] <- stats::rexp(1, rate = hazard_survival)
   }

   censoring_times <- stats::rexp(n_total, rate = censoring_rate)
   assignment_times <- stats::rexp(n_total, rate = treatment_rate)
   cap_time <- base::max(L, base::max(survival_times, censoring_times)) + 5
   assignment_times[assignment_times > cap_time] <- stats::runif(base::sum(assignment_times > cap_time), 0, cap_time - 5)

   X <- base::pmin(survival_times, censoring_times)
   deltaD <- base::as.numeric(survival_times <= censoring_times)

   df_datA_status <- base::sample(base::c("active", "control"), n_total, replace = TRUE, prob = base::c(0.5, 0.5))
   df_datA <- base::data.frame(ID = ID, A = assignment_times, status = df_datA_status)

   dat <- base::data.frame(ID = ID, entry = entry, X = X, deltaD = deltaD, deltaT = 0, Za = Za)
   dat <- base::cbind(dat, Zb_df)

   for (i in 1:n_total) {
      if (df_datA$status[i] == "active" && df_datA$A[i] < dat$X[i]) {
         dat$deltaT[i] <- 1
      }
   }
   base::return(base::list(dat = dat, datA = df_datA))
}

#' Perform Dynamic Prediction Test and Return P-value
#'
#' This helper function simulates data, runs the `dp` (dynamic prediction) function,
#' and performs a Wald test on the `Za` coefficient (representing the group effect)
#' to return a p-value.
#'
#' @param n_per_group Numeric, number of subjects in each group for simulation.
#' @param L Numeric, restriction time for RMST.
#' @param K Integer, number of cross-sections.
#' @param interval Numeric, time interval between cross-sections.
#' @param base_rate_survival Numeric, baseline hazard for survival.
#' @param group_effect_Za Numeric, effect of `Za` (group) on survival.
#' @param Zb_effect Numeric, effect of `Zb` on survival.
#' @param censoring_rate Numeric, baseline hazard for censoring.
#' @param treatment_rate Numeric, baseline hazard for treatment.
#' @param dp_method Character, `method` argument for `dp` function ("link" or "stratified").
#' @param dp_link Character, `link` argument for `dp` function ("linear" or "log").
#' @param dp_stratified Character, `stratified` argument for `dp` function ("add" or "multi").
#' @param dp_weights Character, `weights` argument for `dp` function ("stabilized" or "unstabilized").
#' @return A numeric p-value from the Wald test, or `NA` if the estimation fails.
#' @keywords internal
perform_dp_test <- function(n_per_group, L, K, interval,
                            base_rate_survival, group_effect_Za, Zb_effect,
                            censoring_rate, treatment_rate,
                            dp_method, dp_link, dp_stratified, dp_weights) {

   sim_data <- simulate_dp_data(n_per_group, L, K, interval,
                                base_rate_survival, group_effect_Za, Zb_effect,
                                censoring_rate, treatment_rate)

   result <- dp(dat = sim_data$dat, datA = sim_data$datA, L = L, K = K, interval = interval,
                method = dp_method, link = dp_link, stratified = dp_stratified,
                weights = dp_weights)

   if (base::is.null(result) || base::any(base::is.na(result$betahat)) || base::any(base::is.na(result$se))) {
      return(NA)
   }

   za_beta_name <- "Za"
   if (!za_beta_name %in% base::names(result$betahat)) {
      base::warning(base::paste("Coefficient for", za_beta_name, "not found. Failed run."))
      return(NA)
   }

   betahat_za <- result$betahat[za_beta_name]
   se_za <- result$se[za_beta_name]

   if (base::is.na(se_za) || se_za <= 1e-9 || base::is.infinite(se_za)) {
      return(NA)
   }

   test_statistic <- betahat_za / se_za
   p_value <- 2 * stats::pnorm(base::abs(test_statistic), lower.tail = FALSE)

   base::return(p_value)
}

#' Calculate Power for Dynamic Prediction Model
#'
#' This function calculates the statistical power of a test for the group effect (`Za`)
#' within a dynamic prediction model using Monte Carlo simulation.
#'
#' @param n_per_group Numeric, the number of subjects in each binary group.
#' @param L Numeric, the restriction time for RMST.
#' @param K Integer, the number of cross-sections.
#' @param interval Numeric, the time interval between consecutive cross-sections.
#' @param base_rate_survival Numeric, baseline hazard rate for survival events in simulation.
#' @param group_effect_Za Numeric, the log-hazard ratio (or linear effect) for `Za` (group effect)
#'   on survival under the alternative hypothesis. This is the effect size to detect.
#'   e.g., `log(0.8)` for Group 1 having 0.8 times the hazard of Group 0.
#' @param Zb_effect Numeric, the log-hazard ratio (or linear effect) for `Zb` on survival in simulation.
#' @param censoring_rate Numeric, baseline hazard rate for censoring in simulation.
#' @param treatment_rate Numeric, baseline hazard rate for treatment assignment in simulation.
#' @param dp_method Character, `method` argument for `dp` function ("link" or "stratified").
#' @param dp_link Character, `link` argument for `dp` function ("linear" or "log").
#' @param dp_stratified Character, `stratified` argument for `dp` function ("add" or "multi").
#' @param dp_weights Character, `weights` argument for `dp` function ("stabilized" or "unstabilized").
#' @param alpha Numeric, significance level for the hypothesis test (default: 0.05).
#' @param n_sim Integer, number of Monte Carlo simulations to run (default: 1000).
#' @return A numeric value representing the estimated power of the test (proportion of rejected
#'   null hypotheses), or `NA` if all simulations failed.
#' @export
#' @examples
#' # Example: Calculate power for a specific dynamic prediction scenario
#' set.seed(456)
#' power_result_dp <- calculate_dp_power(
#'   n_per_group = 75,
#'   L = 10, K = 2, interval = 5,
#'   base_rate_survival = 0.1, group_effect_Za = log(0.7), # Detect 30% reduction in hazard for Group 1
#'   Zb_effect = 0.1, censoring_rate = 0.03, treatment_rate = 0.04,
#'   dp_method = "link", dp_link = "linear", dp_stratified = "add", # stratified is ignored if method="link"
#'   dp_weights = "stabilized",
#'   alpha = 0.05, n_sim = 50 # Reduced n_sim for quick example
#' )
#' # print(paste("Estimated DP Power:", round(power_result_dp, 4)))
calculate_dp_power <- function(n_per_group, L, K, interval,
                               base_rate_survival, group_effect_Za, Zb_effect,
                               censoring_rate, treatment_rate,
                               dp_method = c("link", "stratified"), dp_link = c("linear", "log"),
                               dp_stratified = c("add", "multi"), dp_weights = c("stabilized", "unstabilized"),
                               alpha = 0.05, n_sim = 1000) {

   dp_method <- base::match.arg(dp_method)
   dp_link <- base::match.arg(dp_link)
   dp_stratified <- base::match.arg(dp_stratified)
   dp_weights <- base::match.arg(dp_weights)

   rejections <- 0
   num_failed_sims <- 0

   for (i in 1:n_sim) {
      p_value <- perform_dp_test(n_per_group, L, K, interval,
                                 base_rate_survival, group_effect_Za, Zb_effect,
                                 censoring_rate, treatment_rate,
                                 dp_method, dp_link, dp_stratified, dp_weights)

      if (!base::is.na(p_value)) {
         if (p_value < alpha) {
            rejections <- rejections + 1
         }
      } else {
         num_failed_sims <- num_failed_sims + 1
      }
   }

   total_successful_sims <- n_sim - num_failed_sims
   if (total_successful_sims == 0) {
      base::warning("All simulations failed. Power cannot be estimated.")
      power <- NA
   } else {
      power <- rejections / total_successful_sims
   }

   base::message(base::paste("Completed", n_sim, "simulations."))
   base::message(base::paste(num_failed_sims, "simulations failed."))
   base::message(base::paste("Estimated Power:", base::round(power, 4)))

   base::return(power)
}

#' Calculate Sample Size for Dynamic Prediction Model
#'
#' This function estimates the required sample size (per group) for a desired power
#' to detect a specified group effect within a dynamic prediction model.
#' It uses an iterative search approach.
#'
#' @param target_power Numeric, the desired statistical power (e.g., 0.8).
#' @param L Numeric, the restriction time for RMST.
#' @param K Integer, the number of cross-sections.
#' @param interval Numeric, the time interval between consecutive cross-sections.
#' @param base_rate_survival Numeric, baseline hazard rate for survival events in simulation.
#' @param group_effect_Za Numeric, the log-hazard ratio (or linear effect) for `Za` (group effect)
#'   on survival under the alternative hypothesis. This is the effect size to detect.
#' @param Zb_effect Numeric, the log-hazard ratio (or linear effect) for `Zb` on survival in simulation.
#' @param censoring_rate Numeric, baseline hazard rate for censoring in simulation.
#' @param treatment_rate Numeric, baseline hazard rate for treatment assignment in simulation.
#' @param dp_method Character, `method` argument for `dp` function ("link" or "stratified").
#' @param dp_link Character, `link` argument for `dp` function ("linear" or "log").
#' @param dp_stratified Character, `stratified` argument for `dp` function ("add" or "multi").
#' @param dp_weights Character, `weights` argument for `dp` function ("stabilized" or "unstabilized").
#' @param alpha Numeric, significance level for the hypothesis test (default: 0.05).
#' @param n_sim_per_iter Integer, number of Monte Carlo simulations per iteration of sample size search (default: 100).
#' @param tol Numeric, tolerance for power difference (default: 0.01).
#' @param max_iter Integer, maximum iterations for sample size search (default: 50).
#' @return An integer representing the estimated required sample size per group, or `NA` if a solution is not found.
#' @export
#' @examples
#' # Example: Estimate sample size for 80% power
#' set.seed(789)
#' sample_size_dp <- calculate_dp_sample_size(
#'   target_power = 0.8,
#'   L = 10, K = 2, interval = 5,
#'   base_rate_survival = 0.1, group_effect_Za = log(0.7),
#'   Zb_effect = 0.1, censoring_rate = 0.03, treatment_rate = 0.04,
#'   dp_method = "link", dp_link = "linear", dp_stratified = "add",
#'   dp_weights = "stabilized",
#'   alpha = 0.05, n_sim_per_iter = 30, # Reduced n_sim for quick example
#'   max_iter = 10
#' )
#' # print(paste("Estimated DP Sample Size per group:", sample_size_dp))
calculate_dp_sample_size <- function(target_power, L, K, interval,
                                     base_rate_survival, group_effect_Za, Zb_effect,
                                     censoring_rate, treatment_rate,
                                     dp_method = c("link", "stratified"), dp_link = c("linear", "log"),
                                     dp_stratified = c("add", "multi"), dp_weights = c("stabilized", "unstabilized"),
                                     alpha = 0.05, n_sim_per_iter = 100,
                                     tol = 0.01, max_iter = 50) {

   dp_method <- base::match.arg(dp_method)
   dp_link <- base::match.arg(dp_link)
   dp_stratified <- base::match.arg(dp_stratified)
   dp_weights <- base::match.arg(dp_weights)

   current_n_per_group <- 50
   power_achieved <- 0
   iter <- 0

   base::message("Starting sample size estimation...")

   while (base::abs(power_achieved - target_power) > tol && iter < max_iter) {
      iter <- iter + 1
      base::message(base::paste("Iteration", iter, ": Testing n_per_group =", current_n_per_group))

      current_power <- calculate_dp_power(
         n_per_group = current_n_per_group, L = L, K = K, interval = interval,
         base_rate_survival = base_rate_survival, group_effect_Za = group_effect_Za,
         Zb_effect = Zb_effect, censoring_rate = censoring_rate, treatment_rate = treatment_rate,
         dp_method = dp_method, dp_link = dp_link, dp_stratified = dp_stratified,
         dp_weights = dp_weights, alpha = alpha, n_sim = n_sim_per_iter
      )

      if (base::is.na(current_power)) {
         base::warning("Power calculation failed. Adjusting n and re-trying.")
         current_n_per_group <- current_n_per_group + 25
         if (current_n_per_group > 10000) {
            base::warning("Sample size search diverging. Aborting.")
            return(NA)
         }
         next
      }

      power_achieved <- current_power

      # Add safeguard for power_achieved being zero
      power_ratio <- if(power_achieved > 0) (target_power / power_achieved) else 4 # Increase by factor of 2 if power is 0

      if (power_achieved < target_power) {
         current_n_per_group <- base::ceiling(current_n_per_group * power_ratio^(1/2))
      } else {
         current_n_per_group <- base::floor(current_n_per_group * power_ratio^(1/2))
      }

      if (current_n_per_group < 10) current_n_per_group <- 10
      current_n_per_group <- base::max(current_n_per_group, 5)
   }

   if (iter == max_iter) {
      base::warning("Max iterations reached. Consider increasing max_iter.")
      return(NA)
   }

   base::message(base::paste("Estimated sample size per group for", target_power * 100, "% power:", current_n_per_group))
   base::return(current_n_per_group)
}


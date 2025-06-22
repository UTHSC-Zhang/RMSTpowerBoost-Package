
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
dp <- function(dat, datA, L, K, interval, # K = the number of cross-sections; interval = the time between consecutive cross-sections
               method = c("link", "stratified"), link = c("linear", "log"), stratified = c("add", "multi"), weights = c("stabilized", "unstabilized")){

   method <- match.arg(method)
   link <- match.arg(link)
   stratified <- match.arg(stratified)
   weights <- match.arg(weights)

   # CRITICAL FIX: Merge dat and datA at the beginning so 'A' is available throughout
   dat_full <- merge(dat, datA, by = "ID")

   # estimate the censoring model - use dat_full
   modC <- coxph(Surv(X, 1 - (deltaD + deltaT)) ~ Za, data = dat_full, ties = "breslow")

   # estimate the treatment model - use dat_full for the base of datT.Extended
   # The `weights_coxph` definition should use `dat_full` as well.
   dat_full$weights_coxph <- ifelse(dat_full$status == "active", 1, 1e-08)
   dtimes <- sort(unique(with(dat_full, X[deltaT == 1])))
   if (length(dtimes) == 0) { # Handle case where no treatment events occur
      warning("No treatment events (deltaT=1) found to fit treatment model. Returning NA.")
      return(list("betahat" = NA, "se" = NA))
   }
   # Ensuring dtimes is appropriate; if max X < min dtimes, survSplit might error.
   # Filter dtimes to be within observed X range to avoid issues with survSplit.
   dtimes <- dtimes[dtimes > min(dat_full$entry) & dtimes < max(dat_full$X)]
   if (length(dtimes) == 0) {
      warning("No valid treatment event times within observed range to fit treatment model. Returning NA.")
      return(list("betahat" = NA, "se" = NA))
   }

   # Now use dat_full instead of datT for survSplit base
   datT.Extended <- survSplit(Surv(X, deltaT == 1) ~ ., dat_full, cut = dtimes, episode = "tstart")
   datT.Extended$Zt <- ifelse(datT.Extended$X > datT.Extended$A, 1, 0) # 'A' is now available
   modT <- coxph(Surv(tstart, X, event) ~ Za + Zt, data = datT.Extended, weights = weights_coxph,
                 ties = "breslow", timefix = FALSE)

   # create the stacked dataset
   dat.stacked <- data.frame()
   for (k in 1:K){
      # Now dat.temp correctly inherits 'A' from dat_full
      dat.temp <- dat_full # This is the crucial change
      dat.temp$Sik <- interval * k - dat.temp$entry
      dat.temp <- dat.temp %>%
         filter(Sik > 0 & Sik < X & Sik < A) # 'A' is now available

      if (nrow(dat.temp) > 0){
         # Dynamically get the Zb_k column for this cross-section
         zb_col_name_k <- paste0("Zb", k)
         if (!(zb_col_name_k %in% names(dat.temp))) {
            warning(paste("Column", zb_col_name_k, "not found for cross-section", k, ". Skipping this cross-section or check data structure."))
            next # Skip this iteration if the Zb column is genuinely missing
         }

         # Select only the relevant columns and rename Zb_k to Zb
         # This makes a new dataframe with ONLY the columns we need + the correct Zb.
         dat.temp_processed <- dat.temp %>%
            mutate(Xik = X - Sik) %>%
            rowwise() %>%
            mutate(Yik = min(Xik, L), deltaYik = ifelse(Xik <= L, deltaD, 1)) %>%
            ungroup() %>%
            # Select base columns, and then the specific Zb_k as Zb
            # Use any_of() to be robust to columns not existing if future data has fewer Zb columns
            select(ID, entry, X, deltaD, deltaT, Za, A, status, weights_coxph,
                   Sik, Xik, Yik, deltaYik, Zb = !!sym(zb_col_name_k)) # Rename Zb_k to Zb

         dat.stacked <- bind_rows(dat.stacked, cbind.data.frame("CS" = k, dat.temp_processed))
      }
   }

   if (nrow(dat.stacked) == 0) {
      warning("No subjects at risk for any cross-section after stacking. Returning NA.")
      return(list("betahat" = NA, "se" = NA))
   }

   # calculate the inverse probability weights
   dat.temp1 = dat.temp2 = dat.temp3 = dat.stacked # These now contain A, Za, Zt implicitly

   # For censoring weights (WC):
   dat.temp1$X_pred <- dat.temp1$Sik # Time for initial survival prediction
   dat.temp2$X_pred <- dat.temp2$Sik + dat.temp2$Yik # Time for later survival prediction

   predC_Sik <- predict(modC, newdata = dat.temp1, type = "survival")
   predC_Sik_Yik <- predict(modC, newdata = dat.temp2, type = "survival")

   predC_Sik_Yik[predC_Sik_Yik <= 1e-10] <- 1e-10 # Prevent division by zero
   dat.stacked$WC <- predC_Sik / predC_Sik_Yik
   dat.stacked$WC[is.na(dat.stacked$WC) | is.infinite(dat.stacked$WC)] <- 0 # Set problematic weights to 0

   # For treatment weights (WT):
   dat.temp1$tstart_pred <- 0
   dat.temp1$event_pred <- 0 # No event (no treatment) by this time
   dat.temp1$X_pred_modT <- dat.temp1$Sik
   dat.temp1$Zt_pred <- ifelse(dat.temp1$A < dat.temp1$Sik, 1, 0) # Zt is 1 if A < current pred time (Sik)

   dat.temp2$tstart_pred <- 0
   dat.temp2$event_pred <- 0 # No event (no treatment) by this time
   dat.temp2$X_pred_modT <- dat.temp2$Sik + dat.temp2$Yik
   dat.temp2$Zt_pred <- ifelse(dat.temp2$A < (dat.temp2$Sik + dat.temp2$Yik), 1, 0)

   # Create newdata for predict.coxph with correct column names for modT
   # modT expects 'X', 'tstart', 'event', 'Za', 'Zt'
   newdata1_modT <- data.frame(
      X = dat.temp1$X_pred_modT,
      tstart = dat.temp1$tstart_pred,
      event = dat.temp1$event_pred,
      Za = dat.temp1$Za,
      Zt = dat.temp1$Zt_pred
   )
   newdata2_modT <- data.frame(
      X = dat.temp2$X_pred_modT,
      tstart = dat.temp2$tstart_pred,
      event = dat.temp2$event_pred,
      Za = dat.temp2$Za,
      Zt = dat.temp2$Zt_pred
   )

   predT_Sik <- predict(modT, newdata = newdata1_modT, type = "survival")
   predT_Sik_Yik <- predict(modT, newdata = newdata2_modT, type = "survival")

   predT_Sik_Yik[predT_Sik_Yik <= 1e-10] <- 1e-10 # Prevent division by zero

   if (weights == "unstabilized"){
      dat.stacked$WT <- predT_Sik / predT_Sik_Yik
   } else if (weights == "stabilized"){
      dat.temp3$tstart_pred <- 0
      dat.temp3$event_pred <- 0 # No event by this time
      dat.temp3$X_pred_modT <- dat.temp3$Yik # Time for numerator prediction
      dat.temp3$Zt_pred <- ifelse(dat.temp3$A < dat.temp3$Yik, 1, 0) # X here is Yik

      newdata3_modT <- data.frame(
         X = dat.temp3$X_pred_modT,
         tstart = dat.temp3$tstart_pred,
         event = dat.temp3$event_pred,
         Za = dat.temp3$Za,
         Zt = dat.temp3$Zt_pred
      )
      predT_Yik <- predict(modT, newdata = newdata3_modT, type = "survival")
      dat.stacked$WT <- (predT_Sik * predT_Yik) / predT_Sik_Yik
   }
   dat.stacked$WT[is.na(dat.stacked$WT) | is.infinite(dat.stacked$WT)] <- 0 # Set problematic weights to 0

   dat.stacked$W <- dat.stacked$WC * dat.stacked$WT
   dat.stacked$W[is.na(dat.stacked$W) | is.infinite(dat.stacked$W)] <- 0 # Final check on combined weights

   # Ensure weights are not all zero, or sum to zero which would cause issues later
   # Use W * deltaYik as the relevant weight for the outcome model
   if (sum(dat.stacked$W * dat.stacked$deltaYik) == 0) {
      warning("Sum of effective weights for outcome model is zero. Returning NA.")
      return(list("betahat" = NA, "se" = NA))
   }


   # point estimate and standard error of betahat
   if (method == "link"){
      dat.stacked$Wgee <- dat.stacked$W * dat.stacked$deltaYik # Use W * deltaYik as weights
      if (link == "linear"){
         family_val = gaussian(link = "identity")
      } else if (link == "log"){
         family_val = poisson(link = "log")
      }
      f <- as.formula("Yik ~ Za + Zb : factor(CS) + factor(CS) - 1")

      model_fit <- tryCatch({
         suppressWarnings(glm(f, data = dat.stacked, family = family_val, weights = Wgee))
      }, error = function(e) {
         warning(paste("GLM fitting failed for link method:", e$message))
         return(NULL)
      })

      if (is.null(model_fit)) {
         return(list("betahat" = NA, "se" = NA))
      }
      betahat <- coef(model_fit)

      A <- matrix(0, ncol = length(betahat), nrow = length(betahat))
      Z <- model.matrix(f, data = dat.stacked)

      if (link == "linear"){
         for (i in 1:nrow(dat.stacked)){
            Zik <- as.matrix(Z[i, ])
            A = A + dat.stacked$W[i] * dat.stacked$deltaYik[i] * Zik %*% t(Zik)
         }
         e <- Z * as.vector(dat.stacked$W * dat.stacked$deltaYik * (dat.stacked$Yik - as.matrix(Z) %*% betahat))
      } else if (link == "log"){
         for (i in 1:nrow(dat.stacked)){
            Zik <- as.matrix(Z[i, ])
            A = A + dat.stacked$W[i] * dat.stacked$deltaYik[i] * Zik %*% t(Zik) * as.numeric(exp(t(betahat) %*% Zik))
         }
         e <- Z * as.vector(dat.stacked$W * dat.stacked$deltaYik * (dat.stacked$Yik - exp(as.matrix(Z) %*% betahat)))
      }
      E <- cbind.data.frame("ID" = dat.stacked$ID, e)
      Ei <- aggregate(. ~ ID, E, sum) # Sum contributions by ID
      B <- matrix(0, ncol = length(betahat), nrow = length(betahat))
      for (i in 1:nrow(Ei)){
         ei <- t(as.matrix(Ei[i, -1])) # Ensure ei is column vector for outer product
         B = B + ei %*% t(ei)
      }
      # Handle singular A matrix for variance calculation
      if (det(A) == 0) {
         warning("Matrix A is singular for link method. Cannot estimate variance. Returning NA.")
         return(list("betahat" = betahat, "se" = rep(NA, length(betahat))))
      }
      var <- diag(solve(A) %*% B %*% t(solve(A)))
      se <- sqrt(var)

   } else if (method == "stratified"){
      dat.stacked$Wt <- dat.stacked$W * dat.stacked$deltaYik
      dat.temp_strat <- dat.stacked[!(is.na(dat.stacked$Wt)|dat.stacked$Wt == 0), ]

      if (nrow(dat.temp_strat) == 0) {
         warning("No valid weighted observations for stratified method. Returning NA.")
         return(list("betahat" = NA, "se" = NA))
      }

      if (stratified == "add"){
         Zbar_list <- lapply(split(dat.temp_strat, dat.temp_strat$CS), function(x){
            if(sum(x$Wt) == 0) return(c(Za=0, Zb=0))
            c(Za = weighted.mean(x$Za, x$Wt), Zb = weighted.mean(x$Zb, x$Wt))
         })
         Zbar <- do.call(rbind, Zbar_list)
         rownames(Zbar) <- unique(dat.temp_strat$CS)

         # Corrected: Use values from dat.temp_strat for centering relative to Zbar
         # and then ensure Zres is aligned with dat.stacked for subsequent calculations
         Za_centered_full <- dat.stacked$Za - Zbar[match(dat.stacked$CS, rownames(Zbar)), "Za"]
         Zb_centered_full <- dat.stacked$Zb - Zbar[match(dat.stacked$CS, rownames(Zbar)), "Zb"]
         Zres_full <- as.matrix(data.frame(Za=Za_centered_full, Zb=Zb_centered_full)) # This Zres is for the full dat.stacked

         model_fit_lm <- tryCatch({
            # Ensure only rows present in dat.temp_strat contribute to the lm
            lm_data <- dat.stacked %>% filter(ID %in% dat.temp_strat$ID) # Filter dat.stacked
            # Ensure Zres_full is indexed by row number, not ID, if it was created for full dat.stacked
            # The safer way is to create Zres only for the lm_data part.
            lm_Za_centered <- lm_data$Za - Zbar[match(lm_data$CS, rownames(Zbar)), "Za"]
            lm_Zb_centered <- lm_data$Zb - Zbar[match(lm_data$CS, rownames(Zbar)), "Zb"]
            lm_Zres <- as.matrix(data.frame(Za=lm_Za_centered, Zb=lm_Zb_centered))

            lm(lm_data$Yik ~ lm_Zres[, "Za"] + lm_Zres[, "Zb"] - 1, weights = lm_data$Wt)
         }, error = function(e) {
            warning(paste("LM fitting failed for stratified-additive method:", e$message))
            return(NULL)
         })

         if (is.null(model_fit_lm)) {
            return(list("betahat" = NA, "se" = NA))
         }
         betahat <- coef(model_fit_lm)
         names(betahat) <- c("Za", "Zb")

         A <- matrix(0, ncol = length(betahat), nrow = length(betahat))
         # Loop through the *original* dat.stacked rows, but apply weights from Wt for filtered subjects
         for (i in 1:nrow(dat.stacked)){
            if (dat.stacked$ID[i] %in% dat.temp_strat$ID) { # Only include subjects that passed initial filter
               # Re-calculate Zres for this specific row for A matrix contribution
               current_Za_centered <- dat.stacked$Za[i] - Zbar[match(dat.stacked$CS[i], rownames(Zbar)), "Za"]
               current_Zb_centered <- dat.stacked$Zb[i] - Zbar[match(dat.stacked$CS[i], rownames(Zbar)), "Zb"]
               current_Zres <- as.matrix(c(current_Za_centered, current_Zb_centered))

               A = A + current_Zres %*% t(current_Zres) * dat.stacked$Wt[i] / nrow(dat_full) # Use nrow(dat_full) for overall N
            }
         }
         if (det(A) == 0) {
            warning("Matrix A is singular for stratified additive method. Cannot estimate variance. Returning NA.")
            return(list("betahat" = betahat, "se" = rep(NA, length(betahat))))
         }

         # Estimate stratum-specific mean RMST after removing covariate effects
         # Ensure Za and Zb used here are from the original dat.stacked, not centered ones
         dat.stacked$Yres_temp <- dat.stacked$Yik - (dat.stacked$Za * betahat["Za"] + dat.stacked$Zb * betahat["Zb"])
         mu0k_list <- lapply(split(dat.stacked, dat.stacked$CS), function(x){
            if(sum(x$Wt) == 0) return(NA)
            weighted.mean(x$Yres_temp, x$Wt)
         })
         mu0k <- unlist(mu0k_list)

         # For e: need Zres for full dat.stacked and scalar part based on that
         e_matrix_part <- Zres_full # Use the Zres_full calculated before lm_data filter
         e_scalar_part <- as.vector(dat.stacked$Wt * (dat.stacked$Yik - mu0k[match(dat.stacked$CS, names(mu0k))] -
                                                         (dat.stacked$Za * betahat["Za"] + dat.stacked$Zb * betahat["Zb"])))
         e <- e_matrix_part * e_scalar_part
         e[is.na(e)] <- 0

         E <- cbind.data.frame("ID" = dat.stacked$ID, e)
         Ei <- aggregate(. ~ ID, E, sum)
         B <- matrix(0, ncol = length(betahat), nrow = length(betahat))
         for (i in 1:nrow(Ei)){
            ei <- t(as.matrix(Ei[i, -1]))
            B = B + ei %*% t(ei) / nrow(dat_full) # Use nrow(dat_full)
         }
         var <- diag(solve(A) %*% B %*% t(solve(A))) / nrow(dat_full) # Use nrow(dat_full)
         se <- sqrt(var)

      } else if (stratified == "multi"){
         dat.stacked$W.new <- dat.stacked$Wt * dat.stacked$Yik
         dat.temp <- dat.stacked
         dat.temp <- dat.temp[!(is.na(dat.temp$W.new)|dat.temp$W.new == 0), ]

         if (nrow(dat.temp) == 0) {
            warning("No valid weighted observations for stratified-multi method. Returning NA.")
            return(list("betahat" = NA, "se" = NA))
         }

         model_fit_cox <- tryCatch({
            coxph(Surv(rep(1, nrow(dat.temp)), rep(1, nrow(dat.temp))) ~ Za + Zb + offset(-log(Yik)) + strata(CS),
                  data = dat.temp, weights = W.new, ties = "breslow", id = ID)
         }, error = function(e) {
            warning(paste("Coxph fitting failed for stratified-multi method:", e$message))
            return(NULL)
         })

         if (is.null(model_fit_cox)) {
            return(list("betahat" = NA, "se" = NA))
         }
         betahat <- coefficients(model_fit_cox)

         Z <- dat.stacked[, c("Za", "Zb")]
         if (!is.numeric(Z$Za) || !is.numeric(Z$Zb)) {
            warning("Za or Zb not numeric for stratified-multi variance calculation. Returning NA.")
            return(list("betahat" = NA, "se" = NA))
         }

         Z_matrix <- as.matrix(Z)
         Zbar_Wt <- dat.stacked$Wt * exp(Z_matrix %*% betahat)
         S0 <- tapply(Zbar_Wt, as.factor(dat.stacked$CS), sum)
         S0[is.na(S0)] <- 1e-10

         S1_agg <- aggregate(Z_matrix * Zbar_Wt, by = list(dat.stacked$CS), sum)
         S1 <- as.matrix(S1_agg[,-1])
         rownames(S1) <- S1_agg[,1]
         S1[is.na(S1)] <- 0

         S2_list <- lapply(split(data.frame(Z=Z_matrix, W_new_val=Zbar_Wt), dat.stacked$CS), function(x_stratum) {
            stratum_Z <- as.matrix(x_stratum[, 1:ncol(Z_matrix)])
            stratum_W <- x_stratum$W_new_val
            if (nrow(stratum_Z) == 0 || sum(stratum_W) == 0) return(matrix(0, nrow=ncol(Z_matrix), ncol=ncol(Z_matrix)))
            sum_S2 <- matrix(0, nrow = ncol(stratum_Z), ncol = ncol(stratum_Z))
            for (i in 1:nrow(stratum_Z)) {
               sum_S2 <- sum_S2 + stratum_Z[i,] %*% t(stratum_Z[i,]) * stratum_W[i]
            }
            sum_S2
         })
         S2_flat <- do.call(rbind, lapply(S2_list, function(mat) as.vector(mat)))
         S2 <- matrix(S2_flat, nrow=length(S2_list), byrow=TRUE)
         rownames(S2) <- names(S2_list)
         colnames(S2) <- paste0(rep(colnames(Z_matrix), each = ncol(Z_matrix)), rep(colnames(Z_matrix), times = ncol(Z_matrix)))
         S2[is.na(S2)] <- 0

         mu0k <- tapply(dat.stacked$W.new, as.factor(dat.stacked$CS), sum) / S0
         mu0k[is.na(mu0k)] <- 0

         Sbar <- S1 / S0[match(rownames(S1), names(S0))]
         Sbar[is.na(Sbar)] <- 0

         Sk_denom_adj <- S0[match(rownames(S2), names(S0))]
         Sk_denom_adj[is.na(Sk_denom_adj) | Sk_denom_adj == 0] <- 1e-10
         Sk_vec <- S2 / Sk_denom_adj -
            t(apply(as.matrix(Sbar), 1, tcrossprod))

         Sk_blocks <- lapply(1:nrow(Sk_vec), function(i) matrix(Sk_vec[i,], nrow=ncol(Z_matrix), ncol=ncol(Z_matrix), byrow=TRUE))
         names(Sk_blocks) <- rownames(Sk_vec)

         A_contrib <- matrix(0, nrow=nrow(dat.stacked), ncol=length(betahat)^2)
         for(i in 1:nrow(dat.stacked)){
            cs_val <- as.character(dat.stacked$CS[i])
            if (cs_val %in% names(Sk_blocks)) {
               current_Sk_block <- Sk_blocks[[cs_val]]
               A_contrib[i,] <- as.vector(current_Sk_block * dat.stacked$W.new[i])
            }
         }
         A_agg <- aggregate(A_contrib, by = list(ID = dat.stacked$ID), sum)[,-1]
         A <- matrix(colSums(A_agg) / nrow(dat_full),
                     nrow = length(betahat), ncol = length(betahat), byrow = TRUE)
         if (det(A) == 0) {
            warning("Matrix A is singular for stratified multiplicative method. Cannot estimate variance. Returning NA.")
            return(list("betahat" = betahat, "se" = rep(NA, length(betahat))))
         }

         e_term1_val <- Z_matrix - Sbar[match(dat.stacked$CS, rownames(Sbar)), ]
         e_term1_val[is.na(e_term1_val)] <- 0

         e_term2_val <- (dat.stacked$Yik - mu0k[match(dat.stacked$CS, names(mu0k))] * as.vector(exp(Z_matrix %*% betahat)))
         e_term2_val[is.na(e_term2_val)] <- 0

         e <- dat.stacked$Wt * e_term1_val * as.vector(e_term2_val)
         e[is.na(e)] <- 0

         E <- cbind.data.frame("ID" = dat.stacked$ID, e)
         Ei <- aggregate(. ~ ID, E, sum)
         B <- matrix(0, ncol = length(betahat), nrow = length(betahat))
         for (i in 1:nrow(Ei)){
            ei <- t(as.matrix(Ei[i, -1]))
            B = B + ei %*% t(ei) / nrow(dat_full)
         }
         var <- diag(solve(A) %*% B %*% t(solve(A))) / nrow(dat_full)
         se <- sqrt(var)
      }
   }

   if (!is.null(betahat) && !is.null(se) && length(betahat) == length(se)) {
      names(se) <- names(betahat)
   } else {
      warning("Betahat or SE calculation resulted in mismatch or NAs. Returning NA.")
      return(list("betahat" = NA, "se" = NA))
   }

   return(list("betahat" = betahat, "se" = se))
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
   entry <- rep(0, n_total) # Assuming everyone starts at time 0

   # Za as binary group (0 or 1)
   Za <- rep(c(0, 1), each = n_per_group)

   # Zb as a continuous covariate, duplicated for Zb1...ZbK
   Zb_base <- rnorm(n_total, mean = 0, sd = 1)
   Zb_cols <- paste0("Zb", 1:K)
   Zb_df <- as.data.frame(matrix(Zb_base, ncol = K, nrow = n_total, byrow = FALSE))
   colnames(Zb_df) <- Zb_cols

   # Simulate survival times (X)
   survival_times <- numeric(n_total)
   for (i in 1:n_total) {
      # Hazard influenced by Za (group) and Zb (using Zb_base)
      hazard_survival <- base_rate_survival * exp(group_effect_Za * Za[i] + Zb_effect * Zb_base[i])
      survival_times[i] <- rexp(1, rate = hazard_survival)
   }

   # Simulate censoring times (C)
   censoring_times <- rexp(n_total, rate = censoring_rate)

   # Simulate treatment assignment times (A)
   # This influences when a subject becomes 'active' or 'control'
   assignment_times <- rexp(n_total, rate = treatment_rate) # Time until assignment
   # Make sure assignment_times are less than X or C, or within a reasonable range for dp logic.
   assignment_times[assignment_times > max(L, max(survival_times, censoring_times)) + 5] <- # Cap them
      runif(sum(assignment_times > max(L, max(survival_times, censoring_times)) + 5), 0, max(L, max(survival_times, censoring_times)))

   # Determine observed X, deltaD, deltaT
   # Event hierarchy: Death > Treatment > Censoring > End of study
   # For simplicity, if death happens first, it's death. Else if treatment happens, it's treatment. Else censoring.
   # A is the time of potential treatment. If treatment occurs, deltaT=1.
   deltaD <- rep(0, n_total) # Death indicator
   deltaT <- rep(0, n_total) # Treatment indicator (receipt)

   # Determine X (observed time) and event indicators
   # X is min(Survival, Censoring)
   X <- pmin(survival_times, censoring_times)

   # Decide deltaD (death event)
   deltaD <- as.numeric(survival_times <= censoring_times)

   # Define actual treatment status (`status` for datA) and assignment time (`A` for datA)
   # For datA, 'A' is the assignment time. 'status' determines if they are in active arm or control.
   # Let's say: if assignment_times < X, they were assigned.
   # Half of those assigned go to 'active', half to 'control' (randomly)
   df_datA_status <- sample(c("active", "control"), n_total, replace = TRUE, prob = c(0.5, 0.5))

   # Ensure A (assignment time) in datA is meaningful within the context of X and censoring.
   # For simplicity, let's set A to the actual `assignment_times` and let `dp` filter based on `Sik < A`.
   df_datA <- data.frame(ID = ID, A = assignment_times, status = df_datA_status)

   # Final dat construction
   dat <- data.frame(ID = ID, entry = entry, X = X, deltaD = deltaD, deltaT = deltaT, Za = Za)
   dat <- bind_cols(dat, Zb_df) # Add Zb1...ZbK columns

   # Correct deltaT: if a subject was assigned treatment before their observed X, mark deltaT=1
   # This `deltaT` is for the `Surv(X, deltaT == 1)` in `modT` of `dp` function.
   # It marks if treatment occurred as an event (not just assignment).
   # Assuming `df_datA$A` represents time of treatment *initiation* if `status == "active"`.
   # If `status == "control"`, they effectively never receive treatment for deltaT.
   # This makes deltaT 1 if the assignment time for 'active' is before X, and they are indeed 'active'.
   for (i in 1:n_total) {
      if (df_datA$status[i] == "active" && df_datA$A[i] < dat$X[i]) {
         dat$deltaT[i] <- 1
      } else {
         dat$deltaT[i] <- 0
      }
   }


   return(list(dat = dat, datA = df_datA))
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

   dat <- sim_data$dat
   datA <- sim_data$datA

   # Call the main dp function
   result <- dp(dat = dat, datA = datA, L = L, K = K, interval = interval,
                method = dp_method, link = dp_link, stratified = dp_stratified,
                weights = dp_weights)

   # Check if dp function returned NA (indicating failure)
   if (is.null(result) || any(is.na(result$betahat)) || any(is.na(result$se))) {
      return(NA)
   }

   # Extract the coefficient and SE for the 'Za' (group) variable
   za_beta_name <- "Za"
   if (!za_beta_name %in% names(result$betahat)) {
      warning(paste("Coefficient for", za_beta_name, "not found in dp result. This simulation run failed."))
      return(NA)
   }

   betahat_za <- result$betahat[za_beta_name]
   se_za <- result$se[za_beta_name]

   # Perform Wald test
   if (is.na(se_za) || se_za <= 1e-9 || is.infinite(se_za)) {
      return(NA) # Problematic SE
   }

   test_statistic <- betahat_za / se_za
   p_value <- 2 * pnorm(abs(test_statistic), lower.tail = FALSE)

   return(p_value)
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

   dp_method <- match.arg(dp_method)
   dp_link <- match.arg(dp_link)
   dp_stratified <- match.arg(dp_stratified)
   dp_weights <- match.arg(dp_weights)

   rejections <- 0
   num_failed_sims <- 0

   for (i in 1:n_sim) {
      p_value <- perform_dp_test(n_per_group, L, K, interval,
                                 base_rate_survival, group_effect_Za, Zb_effect,
                                 censoring_rate, treatment_rate,
                                 dp_method, dp_link, dp_stratified, dp_weights)

      if (!is.na(p_value)) {
         if (p_value < alpha) {
            rejections <- rejections + 1
         }
      } else {
         num_failed_sims <- num_failed_sims + 1
      }
   }

   total_successful_sims <- n_sim - num_failed_sims
   if (total_successful_sims == 0) {
      warning("All simulations failed to produce a valid p-value. Power cannot be estimated.")
      power <- NA
   } else {
      power <- rejections / total_successful_sims
   }

   message(paste("Completed", n_sim, "simulations."))
   message(paste(num_failed_sims, "simulations failed to produce a valid p-value."))
   message(paste("Estimated Power:", round(power, 4)))

   return(power)
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

   dp_method <- match.arg(dp_method)
   dp_link <- match.arg(dp_link)
   dp_stratified <- match.arg(dp_stratified)
   dp_weights <- match.arg(dp_weights)

   # Initial guess for n_per_group
   current_n_per_group <- 50
   power_achieved <- 0
   iter <- 0

   message("Starting sample size estimation for dynamic prediction...")

   while (abs(power_achieved - target_power) > tol && iter < max_iter) {
      iter <- iter + 1
      message(paste("Iteration", iter, ": Testing n_per_group =", current_n_per_group))

      current_power <- calculate_dp_power(
         n_per_group = current_n_per_group,
         L = L, K = K, interval = interval,
         base_rate_survival = base_rate_survival,
         group_effect_Za = group_effect_Za,
         Zb_effect = Zb_effect,
         censoring_rate = censoring_rate,
         treatment_rate = treatment_rate,
         dp_method = dp_method, dp_link = dp_link,
         dp_stratified = dp_stratified, dp_weights = dp_weights,
         alpha = alpha,
         n_sim = n_sim_per_iter # Use fewer simulations per iteration for speed
      )

      if (is.na(current_power)) {
         warning("Power calculation failed for current n_per_group. Adjusting n_per_group or parameters.")
         current_n_per_group <- current_n_per_group + 25 # Try a larger N
         if (current_n_per_group > 10000) { # Prevent infinite loop
            warning("Sample size search diverging. Consider adjusting simulation parameters or stopping criteria.")
            return(NA)
         }
         next # Skip to next iteration
      }

      power_achieved <- current_power

      if (power_achieved < target_power) {
         current_n_per_group <- ceiling(current_n_per_group * (target_power / power_achieved)^(1/2)) # Increase N, faster adjustment
      } else {
         current_n_per_group <- floor(current_n_per_group * (target_power / power_achieved)^(1/2)) # Decrease N
      }
      # Ensure n_per_group is at least 10 for meaningful simulation
      if (current_n_per_group < 10) current_n_per_group <- 10

      # Prevent large jumps or negative values
      current_n_per_group <- max(current_n_per_group, 5) # Minimum reasonable n_per_group
   }

   if (iter == max_iter) {
      warning("Max iterations reached without achieving target power. Consider increasing max_iter or adjusting parameters.")
      return(NA)
   }

   message(paste("Estimated sample size per group for", target_power * 100, "% power:", current_n_per_group))
   return(current_n_per_group)
}

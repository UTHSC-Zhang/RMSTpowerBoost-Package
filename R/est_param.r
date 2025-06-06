
##########################################################################################
#
#  Purpose:            Function to estimate model parameters and the corresponding 
#                      asymptotic variance
#
#  Input data files:   a data frame containing baseline covariates, time-to-event outcome, 
#                      and stratum
#
#  Output data files:  a list of a point estimate vector of model parameters, 
#                      an asymptotic variance vector, and a data frame of stratum mean RMST
#
##########################################################################################

library(survival)

estBeta <- function(dat, Xname = "X", deltaXname = "deltaX", 
                    Znames = c("Z1", "Z2"), ZCnames = c("Z1", "Z2"), 
                    strname = "age", L){
  ## ARGUMENTS OF THE FUNCTION
  ## - dat: an object of class "data.frame"
  ## - Xname: a string representing the name of the observation time
  ## - deltaXname: a string representing the name of the death indicator
  ## - Znames: a vector of names of the covariates predictive of survival
  ## - ZCnames: a vector of names of the covariates predictive of censoring
  ## - strname: a string representing the name of the stratum
  ## - L: restriction time
  
  # convert time-to-event outcome to RMST
  dat$Y <- pmin(dat[[Xname]], L)
  dat$deltaY <- ifelse(dat[[deltaXname]] == 1, 1, ifelse(L < dat[[Xname]], 1, 0))
  
  # remove strata with no events
  dat <- dat[!dat[[strname]] %in% names(which(tapply(dat$deltaY, dat[[strname]], sum) == 0)), ]
  
  # estimate IPCWs
  modC <- coxph(as.formula(paste("Surv(", Xname, ", 1-", deltaXname, ") ~ ", 
                                 paste0(ZCnames, collapse = " + "), " + strata(", strname, ")")), 
                data = dat, ties = "breslow")
  dat.new <- dat
  dat.new[[Xname]] <- dat.new$Y
  dat$W <- dat$deltaY/predict(modC, newdata = dat.new, type = "survival")
  
  # estimate betahat
  n = nrow(dat)
  p = length(Znames)
  J = length(unique(dat[[strname]]))
  strvalue <- unique(dat[[strname]])
  Zbar <-  matrix(ncol = p, nrow = J)
  for (j in 1:J){
    Zbar[j, ] <- colSums(as.data.frame(dat[dat[[strname]] == strvalue[j], Znames]) * 
                           dat$W[dat[[strname]] == strvalue[j]]) / 
      sum(dat$W[dat[[strname]] == strvalue[j]])
  }
  Z <- dat[, Znames]
  Zres <- Z - Zbar[match(dat[[strname]], strvalue), ]
  B <- apply(Zres * dat$W * dat$Y, 2, sum)
  A <- matrix(0, ncol = p, nrow = p)
  for (i in 1:n){
    A = A + t(as.matrix(Zres[i, ])) %*% as.matrix(Zres[i, ]) * dat$W[i]
  }
  betahat <- solve(A) %*% B
  
  # estimate stratum mean RMST
  mu0 <- rep(0, J)
  for(j in 1:J){
    dat.temp <- dat[dat[[strname]] == strvalue[j], ]
    mu0[j] <- sum(dat.temp$W * (dat.temp$Y - as.matrix(dat.temp[, Znames]) %*% betahat)) / 
      sum(dat.temp$W)
  }
  
  # estimate asymptotic variance
  Avar <- matrix(0, ncol = p, nrow = p)
  for (i in 1:n){
    Avar = Avar + t(as.matrix(Zres[i, ])) %*% as.matrix(Zres[i, ]) * dat$W[i] / n
  }
  Bvar <- matrix(0, ncol = p, nrow = p)
  for (i in 1:nrow(Zres)){
    temp <- Zres[i, ] * dat$W[i] * 
      (dat$Y[i] - mu0[match(dat[[strname]][i], strvalue)] - as.matrix(Z[i, ]) %*% betahat)
    Bvar = Bvar + t(as.matrix(temp)) %*% as.matrix(temp) / n
  }
  var <- diag(t(solve(Avar)) %*% Bvar %*% solve(Avar)) / n
  
  return(list("betahat" = betahat, "var" = var, 
              "mu0_df" = cbind.data.frame("stratum" = strvalue, "mu0" = mu0)))
}

##########################################################################################
#
#  Function to perform dynamic prediction based on stratified or link models
#
##########################################################################################

dp <- function(dat, datA, L, K, interval,  # K = the number of cross-sections; interval = the time between consecutive cross-sections
               method = c("link", "stratified"), link = c("linear", "log"), stratified = c("add", "multi"), weights = c("stabilized", "unstabilized")){
  
  # estimate the censoring model
  datC <- dat
  modC <- coxph(Surv(X, 1-(deltaD+deltaT)) ~ Za, data = datC, ties = "breslow")
  
  # estimate the treatment model
  datT <- merge(dat, datA, by = "ID")
  datT$weights <- ifelse(datT$status == "active", 1, 1e-08)
  dtimes <- sort(unique(with(datT, X[deltaT == 1])))
  datT.Extended <- survSplit(Surv(X, deltaT == 1) ~ ., datT, cut = dtimes)
  datT.Extended[, "Zt"] <- ifelse(datT.Extended$X > datT.Extended$A, 1, 0)
  modT <- coxph(Surv(tstart, X, event) ~ Za + Zt, data = datT.Extended, weights = weights, 
                ties = "breslow", timefix = FALSE)
  
  # create the stacked dataset
  dat.stacked <- data.frame()
  for (k in 1:K){
    # select subjects who are at-risk and treatment-eligible at Sik
    dat.temp <- dat
    dat.temp$Sik <- interval*k - dat.temp$entry
    dat.temp <- dat.temp %>%
      filter(Sik > 0 & Sik < X & Sik < A)
    if (nrow(dat.temp) > 0){
      dat.temp <- dat.temp %>%
        mutate(Xik = X - Sik) %>%
        rowwise() %>%
        mutate(Yik = min(Xik, L), deltaYik = ifelse(Xik <= L, deltaD, 1)) %>%
        select(-paste("Zb", (1:K)[-k], sep = ""))
      colnames(dat.temp)[grepl(paste("Zb", k, sep = ""), colnames(dat.temp))] <- "Zb"
      dat.stacked <- rbind.data.frame(dat.stacked, cbind.data.frame("CS" = k, dat.temp))
    }
    else{
      dat.stacked <- dat.stacked
    }
  }
  
  # calculate the inverse probability weights
  dat.temp1 = dat.temp2 = dat.temp3 = dat.stacked
  dat.temp1$X <- dat.temp1$Sik
  dat.temp2$X <- dat.temp2$Sik + dat.temp2$Yik
  dat.stacked$WC <- predict(modC, newdata = dat.temp1, type = "survival")/predict(modC, newdata = dat.temp2, type = "survival")
  
  dat.temp1$Zt <- ifelse(dat.temp1$A < dat.temp1$Sik, 1, 0)
  dat.temp1$tstart <- 0
  # dat.temp1$X <- dat.temp1$Sik
  dat.temp1$event <- dat.temp1$deltaT
  dat.temp2$Zt <- ifelse(dat.temp2$A < dat.temp2$Sik + dat.temp2$Yik, 1, 0)
  dat.temp2$tstart <- 0
  # dat.temp2$X <- dat.temp2$Sik + dat.temp2$Yik
  dat.temp2$event <- dat.temp2$deltaT
  dat.temp3$Zt <- ifelse(dat.temp3$A < dat.temp3$Yik, 1, 0)
  dat.temp3$tstart <- 0
  dat.temp3$X <- dat.temp3$Yik
  dat.temp3$event <- dat.temp3$deltaT
  
  if (weights == "unstabilized"){
    dat.stacked$WT <- predict(modT, newdata = dat.temp1, type = "survival")/predict(modT, newdata = dat.temp2, type = "survival")
  }
  else if (weights == "stabilized"){
    dat.stacked$WT <- predict(modT, newdata = dat.temp1, type = "survival")*predict(modT, newdata = dat.temp3, type = "survival")/predict(modT, newdata = dat.temp2, type = "survival")
  }
  dat.stacked$W <- dat.stacked$WC * dat.stacked$WT
  
  # point estimate and standard error of betahat
  if (method == "link"){
    dat.stacked$Wgee <- dat.stacked$W * dat.stacked$deltaYik
    if (link == "linear"){
      family = gaussian(link = "identity")
    }
    else if (link == "log"){
      family = poisson(link = "log")
    }
    f <- as.formula("Yik ~ Za - 1 + Zb : factor(CS) + factor(CS)")
    betahat <- coef(suppressWarnings(glm(f, data = dat.stacked, family = family, weights = Wgee)))
    
    A <- matrix(0, ncol = length(betahat), nrow = length(betahat))
    B <- matrix(0, ncol = length(betahat), nrow = length(betahat))
    CS <- model.matrix(~ factor(CS) - 1, dat.stacked)
    Z <- model.matrix(f, data = dat.stacked)
    if (link == "linear"){
      for (i in 1:nrow(dat.stacked)){
        Zik <- as.matrix(Z[i, ])
        A <- A + dat.stacked$W[i] * dat.stacked$deltaYik[i] * Zik %*% t(Zik)
      }
      e <- Z * as.vector(dat.stacked$W * dat.stacked$deltaYik * (dat.stacked$Yik - as.matrix(Z) %*% betahat))
    }
    else if (link == "log"){
      for (i in 1:nrow(dat.stacked)){
        Zik <- as.matrix(Z[i, ])
        A <- A + dat.stacked$W[i] * dat.stacked$deltaYik[i] * Zik %*% t(Zik) * as.numeric(exp(t(betahat) %*% Zik))
      }
      e <- Z * as.vector(dat.stacked$W * dat.stacked$deltaYik * (dat.stacked$Yik - exp(as.matrix(Z) %*% betahat)))
    }
    E <- cbind.data.frame("ID" = dat.stacked$ID, e)
    Ei <- aggregate(. ~ ID, E, sum)
    for (i in 1:nrow(Ei)){
      ei <- t(Ei[i, -1])
      B <- B + ei %*% t(ei)
    }
    var <- diag(solve(A) %*% B %*% t(solve(A)))
  }
  
  else if (method == "stratified"){
    dat.stacked$Wt <- dat.stacked$W * dat.stacked$deltaYik
    if (stratified == "add"){
      Zbar <- cbind.data.frame("Za" = sapply(split(dat.stacked, dat.stacked$CS), function(x){weighted.mean(x$Za, x$Wt)}), 
                               "Zb" = sapply(split(dat.stacked, dat.stacked$CS), function(x){weighted.mean(x$Zb, x$Wt)}))
      dat.temp <- dat.stacked
      dat.temp$Za = dat.stacked$Za - Zbar[dat.temp$CS, 1]
      dat.temp$Zb = dat.stacked$Zb - Zbar[dat.temp$CS, 2]
      betahat <- coef(lm(Yik ~ Za + Zb - 1, dat = dat.temp, weights = Wt))
      Zres <- dat.temp[, c("Za", "Zb")]
      A <- matrix(0, ncol = length(betahat), nrow = length(betahat))
      for (i in 1:nrow(dat.stacked)){
        A = A + t(as.matrix(Zres[i, ])) %*% as.matrix(Zres[i, ]) * dat.stacked$Wt[i] / nrow(dat)
      }
      dat.temp$Yres <- dat.stacked$Yik - as.matrix(dat.stacked[, c("Za", "Zb")]) %*% betahat
      mu0k <- sapply(split(dat.temp, dat.temp$CS), function(x){weighted.mean(x$Yres, x$Wt)})
      B <- matrix(0, ncol = length(betahat), nrow = length(betahat))
      e <- diag(as.vector(dat.stacked$Wt * (dat.stacked$Yik - mu0k[dat.stacked$CS] - as.matrix(dat.stacked[, c("Za", "Zb")]) %*% betahat))) %*% as.matrix(Zres)
      E <- cbind.data.frame("ID" = dat.stacked$ID, e)
      Ei <- aggregate(. ~ ID, E, sum)
      for (i in 1:nrow(Ei)){
        ei <- t(Ei[i, -1])
        B <- B + ei %*% t(ei) / nrow(dat)
      }
      var <- diag(solve(A) %*% B %*% t(solve(A))) / nrow(dat)
    }
    else if (stratified == "multi"){
      dat.stacked$W.new <- dat.stacked$Wt * dat.stacked$Yik
      dat.temp <- dat.stacked
      dat.temp <- dat.temp[!(is.na(dat.temp$W.new)|dat.temp$W.new == 0), ]
      dat.temp$deltaX.new <- 1
      dat.temp$X.new <- 1
      mod.multi <- coxph(Surv(X.new, deltaX.new) ~ Za + Zb + offset(-log(Yik)) + strata(CS), data = dat.temp, 
                         weights = W.new, ties = "breslow", id = ID)
      betahat <- coefficients(mod.multi)
      Z <- dat.stacked[, c("Za", "Zb")]
      Zbar <- dat.stacked$Wt * exp(as.matrix(Z) %*% betahat)
      S0 <- tapply(Zbar, as.factor(dat.stacked$CS), sum)
      S1 <- aggregate(Z * Zbar, by = list(dat.stacked$CS), sum)[, -1]
      S2 <- aggregate(data.frame(t(apply(as.matrix(Z), 1, tcrossprod))) * Zbar, by = list(dat.stacked$CS), sum)[, -1]
      mu0k <- tapply(dat.stacked$W.new, as.factor(dat.stacked$CS), sum)/S0
      Sbar <- S1/S0
      Sk <- S2/S0 - t(apply(as.matrix(Sbar), 1, tcrossprod))
      A <- matrix(colSums(aggregate(Sk[dat.stacked$CS, ] * dat.stacked$W.new, by = list(dat.stacked$ID), sum)[, -1]) / nrow(dat), 
                  nrow = length(betahat), ncol = length(betahat))
      B <- matrix(0, ncol = length(betahat), nrow = length(betahat))
      e <- dat.stacked$Wt * (Z - Sbar[dat.stacked$CS, ]) * (dat.stacked$Yik - mu0k[dat.stacked$CS] * as.vector(exp(as.matrix(Z) %*% betahat)))
      E <- cbind.data.frame("ID" = dat.stacked$ID, e)
      Ei <- aggregate(. ~ ID, E, sum)
      for (i in 1:nrow(Ei)){
        ei <- t(Ei[i, -1])
        B <- B + ei %*% t(ei) / nrow(dat)
      }
      var <- diag(solve(A) %*% B %*% t(solve(A))) / nrow(dat)
    }
  }
  
  return(list("betahat" = betahat, "se" = sqrt(var)))
}

# rsEfitSplit obtain the estimation of RSICF assuming a single index model on a slited data. The variance estimaor of the coefficiencts are provided.
# splitIndex has three list. $nuisance, on which we fit nuisance parameter. $fit, on which we fit an initial rsfitSplit. $eff, on which we improve by the efficient score
rsEfitSplit <- function(covariate, response, treatment, splitIndex = NULL, propensityModel = 'glmnet', estimatedPropensity = NULL, outcomeModel = 'kernel', estimatedOutcome = NULL, lossType = 'logistic', weights = NULL, tol = 1e-3, propensityFormula=NULL, outcomeFormula=NULL){

  # fit nuisance parameter on 1
  estimatedNuisance <- NULL
  estimatedNuisance$pi <- estimatedPropensity
  if (is.null(estimatedPropensity)){
    data <- list(predictor = covariate, treatment = treatment)
    predictedPropensityAll <- getPropensityModel(data, method = propensityModel, splitIndex = splitIndex, Formula = propensityFormula, predictAll = TRUE)
    estimatedNuisance$pi <- predictedPropensityAll
  }
  estimatedNuisance$S <- estimatedOutcome$control+estimatedOutcome$treatment
  if (is.null(estimatedOutcome)){
    data <- list(predictor = covariate, treatment = treatment, outcome = response)
    predictedOutcomeAll <- getOutcomeModel(data, method = outcomeModel, splitIndex = splitIndex, Formula = outcomeFormula, predictAll = TRUE)
    estimatedNuisance$S <- predictedOutcomeAll$control+predictedOutcomeAll$treatment
  }

  # weights
  if(is.null(weights)){
    weights <- rep(1, times=NROW(covariate))
  }

  # start cross-validation for M_n
  initSolve <- rsfitSolver(covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex, lossType = lossType, weights =weights, tol = tol, m_0=10^3, withM=TRUE)
  m_max <- sum(abs(initSolve$xi))
  m_seq <- m_max*(1-2^(-seq(0.1, 8, length.out = 10)))
  fit <- NULL
  #if (!parallel){
  score <- array(0,c(10,1))
  pi_tmp <- estimatedNuisance$pi[splitIndex$nuisance]
  S_tmp <- estimatedNuisance$S[splitIndex$nuisance]
  response_tmp <- response[splitIndex$nuisance]
  treatment_tmp <- 2*(treatment[splitIndex$nuisance]-0.5)
  weights_tmp <- weights[splitIndex$nuisance]
  for (iter in 1:10){
    fit[[iter]] <- rsfitSolver(covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex, lossType = lossType, weights =weights, tol = tol, m_0=m_seq[iter], withM=TRUE)
    predict_tmp <- predict.rsfitSplit(fit[[iter]], newx = covariate[splitIndex$nuisance,])
    score[iter] <- weighted.mean(response_tmp/(pi_tmp*(treatment_tmp>0)+(1-pi_tmp)*(treatment_tmp<0))*loss(treatment_tmp * predict_tmp, type = lossType)+treatment_tmp/(2*(pi_tmp*(treatment_tmp>0)+(1-pi_tmp)*(treatment_tmp<0)))*(S_tmp-response_tmp)*predict_tmp, weights_tmp)
  }
  score.min <- which.min(score)
  fit0 <- fit[[iter]]
  # predict eta
  eta_trans <- predict.rsfitSplit(fit0, newx = covariate)
  eta<- link(eta_trans, lossType=lossType)

  # reset splitIndex
  splitIndex_eff <- splitIndex
  splitIndex_eff$fit <- splitIndex$eff
  #splitIndex_eff$nuisance <- (splitIndex$nuisance | splitIndex$fit)
  # augmentation
  ratio <- (1-eta)/(1+eta)
  inv_W <- 1/estimatedNuisance$pi * ratio + 1/(1-estimatedNuisance$pi) * 1/ratio
  weights <- inv_W^(-1) * estimatedNuisance$S
  # refit the efficient loss
  initSolve <- rsfitSolver(covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex_eff, lossType = lossType, weights =weights, tol = tol, m_0=10^3, withM=TRUE)
  m_max <- sum(abs(initSolve$xi))
  m_seq <- m_max*(1-2^(-seq(0.1, 8, length.out = 10)))
  fit <- NULL
  #if (!parallel){
  score <- array(0,c(10,1))
  pi_tmp <- estimatedNuisance$pi[splitIndex_eff$nuisance]
  S_tmp <- estimatedNuisance$S[splitIndex_eff$nuisance]
  response_tmp <- response[splitIndex_eff$nuisance]
  treatment_tmp <- 2*(treatment[splitIndex_eff$nuisance]-0.5)
  weights_tmp <- weights[splitIndex_eff$nuisance]
  for (iter in 1:10){
    fit[[iter]] <- rsfitSolver(covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex_eff, lossType = lossType, weights =weights, tol = tol, m_0=m_seq[iter], withM=TRUE)
    predict_tmp <- predict.rsfitSplit(fit[[iter]], newx = covariate[splitIndex_eff$nuisance,])
    score[iter] <- weighted.mean(response_tmp/(pi_tmp*(treatment_tmp>0)+(1-pi_tmp)*(treatment_tmp<0))*loss(treatment_tmp * predict_tmp, type = lossType)+treatment_tmp/(2*(pi_tmp*(treatment_tmp>0)+(1-pi_tmp)*(treatment_tmp<0)))*(S_tmp-response_tmp)*predict_tmp, weights_tmp)
  }
  score.min <- which.min(score)
  fit_hat <- fit[[score.min]]

  # start inference
  predict <- predict.rsfitSplit(fit_hat, newx = covariate)
  w0 <- predict.rsfitSplit(fit_hat, newx = covariate, derivative=TRUE)
  w1 <- response/(estimatedNuisance$pi*treatment+(1-estimatedNuisance$pi)*(1-treatment)) * loss(2*(treatment-0.5) * predict, type = lossType, order = 'hessian') * w0^2 * weights
  w2 <- (2*(treatment-0.5)*response/(estimatedNuisance$pi*treatment+(1-estimatedNuisance$pi)*(1-treatment)) * loss(2*(treatment-0.5) * predict, type = lossType, order = 'derivative')+2*(treatment-0.5)*(estimatedNuisance$S-response)/(2*(estimatedNuisance$pi*treatment+(1-estimatedNuisance$pi)*(1-treatment))))^2 * w0^2 * weights
  sampleIndex <- c(1:NROW(covariate))
  W1 <- matrix(0, NCOL(covariate), NCOL(covariate))
  W2 <- matrix(0, NCOL(covariate), NCOL(covariate))
  for (iter in sampleIndex[splitIndex$fit]){
    W1 <- W1+w1[iter] * (covariate[iter,]) %*% t(covariate[iter,])
    W2 <- W2+w2[iter] * (covariate[iter,]) %*% t(covariate[iter,])
  }
  W1 <- W1/sum(splitIndex$fit)
  W2 <- W2/sum(splitIndex$fit)
  fit_hat$var <- solve(W1) %*% W2 %*% solve(W1)
  fit_hat
}

# rsEfit obtains sysmetric results
rsEfit <- function(covariate, response, treatment, splitIndex = NULL, propensityModel = 'glmnet', estimatedPropensity = NULL, outcomeModel = 'kernel', estimatedOutcome = NULL, lossType = 'logistic', weights = NULL, tol = 1e-3, propensityFormula=NULL, outcomeFormula=NULL, parallel = FALSE){
  n <- NROW(covariate)
  if(is.null(splitIndex)){
    random_index <- sample(c(1,2,3), n, replace = TRUE)
    splitIndex <- NULL
    splitIndex$nuisance <- (random_index == 1)
    splitIndex$fit <- (random_index == 2)
    splitIndex$eff <- (random_index == 3)
  }
  fit <- NULL
  if(!parallel){
    fit[[1]] <- rsEfitSplit(covariate, response, treatment, splitIndex = splitIndex, propensityModel = propensityModel, estimatedPropensity = estimatedPropensity, outcomeModel=outcomeModel, estimatedOutcome = estimatedOutcome, lossType = lossType, weights = weights, tol=tol, propensityFormula = propensityFormula, outcomeFormula = outcomeFormula)
    # reverse the index
    tmp <- splitIndex$eff
    splitIndex$eff <- splitIndex$fit
    splitIndex$fit <- splitIndex$nuisance
    splitIndex$nuisance <- tmp
    fit[[2]] <- rsEfitSplit(covariate, response, treatment, splitIndex = splitIndex, propensityModel = propensityModel, estimatedPropensity = estimatedPropensity, outcomeModel=outcomeModel, estimatedOutcome = estimatedOutcome, lossType = lossType, weights = weights, tol=tol, propensityFormula = propensityFormula, outcomeFormula = outcomeFormula)
    tmp <- splitIndex$eff
    splitIndex$eff <- splitIndex$fit
    splitIndex$fit <- splitIndex$nuisance
    splitIndex$nuisance <- tmp
    fit[[3]] <- rsEfitSplit(covariate, response, treatment, splitIndex = splitIndex, propensityModel = propensityModel, estimatedPropensity = estimatedPropensity, outcomeModel=outcomeModel, estimatedOutcome = estimatedOutcome, lossType = lossType, weights = weights, tol=tol, propensityFormula = propensityFormula, outcomeFormula = outcomeFormula)
  } else {
    aggreIndex <- 0 * splitIndex$nuisance + 1 * splitIndex$fit + 2 * splitIndex$eff
    library(doParallel)
    n_cores <- detectCores(all.tests = FALSE, logical = TRUE)
    cl <- makeCluster(3)
    registerDoParallel(cl)
    fit <- foreach(iter = c(0,1,2))%dopar%{
      splitIndex_local <- NULL
      splitIndex_local$nuisance <- (aggreIndex==iter)
      splitIndex_local$fit <- (aggreIndex==((iter+1)%%3))
      splitIndex_local$eff <- (aggreIndex==((iter+2)%%3))
      rsEfitSplit(covariate, response, treatment, splitIndex = splitIndex_local, propensityModel = propensityModel, estimatedPropensity = estimatedPropensity, outcomeModel=outcomeModel, estimatedOutcome = estimatedOutcome, lossType = lossType,weights = weights, tol=tol, propensityFormula = propensityFormula, outcomeFormula = outcomeFormula)
    }
    stopCluster(cl)
  }
  fit
}

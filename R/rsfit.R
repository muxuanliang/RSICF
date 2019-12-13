# rsfitSolver obtains the optimizer of the propoed optimization given M_n
rsfitSolver <- function(covariate, response, treatment, estimatedNuisance, splitIndex = NULL, lossType = 'logistic', weights = NULL, tol = 10^(-3), m_0=0, withM = TRUE){
  # initialization
  fit_init <- fitLinkLinear(covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex, lossType = lossType, weights = weights, tol = tol)
  iter <- 0

  # iteration start
  fit_last <- fit_init
  fit_last$xi <- 0
  diff <- 1
  while((diff > tol) & (iter < 100)){
    # updateXi
    fit_xi <- updateXi(fit_last, covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex,lossType = lossType, weights = weights, tol = tol, m_0=m_0)
    # updateBeta
    fit_beta <- updateBeta(fit_xi, covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex,lossType = lossType, weights = weights, tol = tol)
    #update
    diff <- max(1-sum(fit_beta$beta*fit_last$beta)/sqrt(sum(fit_beta$beta^2)*sum(fit_last$beta^2)), c(abs(fit_beta$xi-fit_last$xi)))
    iter <- iter+1
    fit_last <- fit_beta
  }
  # final update
  fit <- updateXi(fit_last, covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex,lossType = lossType, weights = weights, tol = tol, m_0=m_0)
  fit
}

# rsfitSplit obtain the estimation of RSICF assuming a single index model on a slited data. The variance estimaor of the coefficiencts are provided.
rsfitSplit <- function(covariate, response, treatment, splitIndex = NULL, propensityModel = 'glmnet', estimatedPropensity = NULL, outcomeModel = 'kernel', estimatedOutcome = NULL, lossType = 'logistic', weights = NULL, tol = 1e-3, propensityFormula=NULL, outcomeFormula=NULL){
  # fit nuisance parameter
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
  # } else {
  #   library(doParallel)
  #   n_cores <- detectCores(all.tests = FALSE, logical = TRUE)
  #   cl <- makeCluster(min(10, n_cores))
  #   registerDoParallel(cl)
  #   score <- foreach(iter = 1:10,.combine = c)%dopar%{
  #     fit <- rsfitSolver(covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex, initLinkType = initLinkType, weights =weights, tol = tol, m_0=m_seq[iter], withM=TRUE)
  #     predict_tmp <- predict.rsfitSplit(fit, newx = covariate[splitIndex$nuisance,])
  #     pi_tmp <- estimatedNuisance$pi[splitIndex$nuisance]
  #     S_tmp <- estimatedNuisance$S[splitIndex$nuisance]
  #     response_tmp <- response[splitIndex$nuisance]
  #     treatment_tmp <- 2*(treatment[splitIndex$nuisance]-0.5)
  #     score <- mean(treatment_tmp*S_tmp/(2*(pi_tmp*(treatment_tmp>0)+(1-pi_tmp)*(treatment_tmp<0))) * link(predict_tmp, type=initLinkType) - response_tmp/(pi_tmp*(treatment_tmp>0)+(1-pi_tmp)*(treatment_tmp<0))*log(1+treatment_tmp*link(predict_tmp, type=initLinkType)))
  #     score
  #   }
  #   stopCluster(cl)
  #}
  score.min <- which.min(score)
  fit_hat <- fit[[score.min]]

  # start inference
  predict <- predict.rsfitSplit(fit_hat, newx = covariate)
  w0 <- predict.rsfitSplit(fit_hat, newx = covariate, derivative=TRUE)
  w1 <- response/(estimatedNuisance$pi*treatment+(1-estimatedNuisance$pi)*(1-treatment)) * loss(2*(treatment-0.5) * predict, type = lossType, order = 'hessian') * w0^2
  w2 <- (2*(treatment-0.5)*response/(estimatedNuisance$pi*treatment+(1-estimatedNuisance$pi)*(1-treatment)) * loss(2*(treatment-0.5) * predict, type = lossType, order = 'derivative')+2*(treatment-0.5)*(estimatedNuisance$S-response)/(2*(estimatedNuisance$pi*treatment+(1-estimatedNuisance$pi)*(1-treatment))))^2 * w0^2
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

# rsfit obtains sysmetric results
rsfit <- function(covariate, response, treatment, splitIndex = NULL, propensityModel = 'glmnet', estimatedPropensity = NULL, outcomeModel = 'kernel', estimatedOutcome = NULL, lossType = 'logistic', weights = NULL, tol = 1e-3, propensityFormula=NULL, outcomeFormula=NULL, parallel = FALSE){
  n <- NROW(covariate)
  if(is.null(splitIndex)){
    random_index <- sample(c(1,2), n, replace = TRUE)
    splitIndex <- NULL
    splitIndex$nuisance <- (random_index == 1)
    splitIndex$fit <- (random_index == 2)
  }
  fit <- NULL
  if(!parallel){
    fit[[1]] <- rsfitSplit(covariate, response, treatment, splitIndex = splitIndex, lossType = lossType, propensityModel = propensityModel, estimatedPropensity = estimatedPropensity, outcomeModel=outcomeModel, estimatedOutcome = estimatedOutcome, weights = weights, tol=tol, propensityFormula = propensityFormula, outcomeFormula = outcomeFormula)
    # reverse the index
    splitIndex$nuisance <- !splitIndex$nuisance
    splitIndex$fit <- !splitIndex$fit
    fit[[2]] <- rsfitSplit(covariate, response, treatment, splitIndex = splitIndex, lossType = lossType, propensityModel = propensityModel, estimatedPropensity = estimatedPropensity, outcomeModel=outcomeModel, estimatedOutcome = estimatedOutcome, weights = weights, tol=tol, propensityFormula = propensityFormula, outcomeFormula = outcomeFormula)
  } else {
    library(doParallel)
    n_cores <- detectCores(all.tests = FALSE, logical = TRUE)
    cl <- makeCluster(2)
    registerDoParallel(cl)
    fit <- foreach(iter = c(TRUE, FALSE))%dopar%{
      splitIndex_local <- NULL
      splitIndex_local$nuisance <- (splitIndex$nuisance==iter)
      splitIndex_local$fit <- (splitIndex$fit==iter)
      rsfitSplit(covariate, response, treatment, splitIndex = splitIndex_local, lossType = lossType, propensityModel = propensityModel, estimatedPropensity = estimatedPropensity, outcomeModel=outcomeModel, estimatedOutcome = estimatedOutcome, weights = weights, tol=tol, propensityFormula = propensityFormula, outcomeFormula = outcomeFormula)
    }
    stopCluster(cl)
  }
  fit
}






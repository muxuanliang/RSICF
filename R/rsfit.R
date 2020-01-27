# rsfitSolver obtains the optimizer of the propoed optimization given M_n
rsfitSolver <- function(covariate, response, treatment, estimatedNuisance, splitIndex = NULL, lossType = 'logistic', weights = NULL, tol = 10^(-3), m_0=0, constraint = TRUE, boundaryPoint=c(-1,1)){
  # initialization
  fit_init <- fitLinkLinear(covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex, lossType = lossType, weights = weights, tol = tol)
  iter <- 0

  # iteration start
  fit_last <- fit_init
  fit_last$xi <- 0
  diff <- 1
  while((diff > tol) & (iter < 100)){
    # updateXi
    fit_xi <- updateXi(fit_last, covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex,lossType = lossType, weights = weights, tol = tol, m_0=m_0, constraint=constraint, boundaryPoint = boundaryPoint)
    # updateBeta
    fit_beta <- updateBeta(fit_xi, covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex,lossType = lossType, weights = weights, tol = tol, withConstraint = TRUE, boundaryPoint = boundaryPoint)
    #update
    diff <- max(abs(fit_beta$beta-fit_last$beta))*max(c(abs(fit_beta$xi-fit_last$xi)))
    iter <- iter+1
    fit_last <- fit_beta
  }
  # final update
  fit <- updateXi(fit_last, covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex,lossType = lossType, weights = weights, tol = tol, m_0=m_0, constraint=constraint, boundaryPoint = boundaryPoint)
  fit_tmp <- updateBeta(fit, covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex,lossType = lossType, weights = weights, tol = tol, withConstraint = FALSE, boundaryPoint = boundaryPoint)
  fit$solution <- fit_tmp$solution$solution
  fit$boundaryPoint <- boundaryPoint
  fit
}

# rsfitSplit obtain the estimation of RSICF assuming a single index model on a slited data. The variance estimaor of the coefficiencts are provided.
rsfitSplit <- function(covariate, response, treatment, splitIndex = NULL, propensityModel = 'glmnet', estimatedPropensity = NULL, outcomeModel = 'kernel', estimatedOutcome = NULL, lossType = 'logistic', weights = NULL, tol = 1e-3, propensityFormula=NULL, outcomeFormula=NULL, constraint=TRUE, boundaryPoint=c(-1,1)){
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

  # solve
  fit_hat <- rsfitSolver(covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex, lossType = lossType, weights =weights, tol = tol, m_0=500, constraint = constraint, boundaryPoint = boundaryPoint)

  # return
  list(fit=fit_hat, covariate=covariate, response=response, treatment=treatment, estimatedNuisance=estimatedNuisance, splitIndex=splitIndex, lossType = lossType, weights =weights)
}

# rsfit obtains sysmetric results
rsfit <- function(covariate, response, treatment, splitIndex = NULL, propensityModel = 'glmnet', estimatedPropensity = NULL, outcomeModel = 'kernel', estimatedOutcome = NULL, lossType = 'logistic', weights = NULL, tol = 1e-3, propensityFormula=NULL, outcomeFormula=NULL, parallel = FALSE, constraint = TRUE, boundaryPoint = c(-1,1)){
  n <- NROW(covariate)
  if(is.null(splitIndex)){
    random_index <- sample(c(1,2,3), n, replace = TRUE)
    splitIndex <- NULL
    splitIndex$nuisance <- (random_index == 1)
    splitIndex$fit <- (random_index == 2)
    splitIndex$infer <- (random_index == 3)
  }
  fit <- NULL
  fit$fit <- NULL
  if(!parallel){
    fit_tmp <- rsfitSplit(covariate, response, treatment, splitIndex = splitIndex, lossType = lossType, propensityModel = propensityModel, estimatedPropensity = estimatedPropensity, outcomeModel=outcomeModel, estimatedOutcome = estimatedOutcome, weights = weights, tol=tol, propensityFormula = propensityFormula, outcomeFormula = outcomeFormula, constraint = constraint, boundaryPoint = boundaryPoint)
    fit$fit[[1]] <- rsInference(fit_tmp, efficient=efficient)
    # rotate the index
    tmp <- splitIndex$nuisance
    splitIndex$nuisance <- splitIndex$infer
    splitIndex$infer <- splitIndex$fit
    splitIndex$fit <- tmp

    fit_tmp <- rsfitSplit(covariate, response, treatment, splitIndex = splitIndex, lossType = lossType, propensityModel = propensityModel, estimatedPropensity = estimatedPropensity, outcomeModel=outcomeModel, estimatedOutcome = estimatedOutcome, weights = weights, tol=tol, propensityFormula = propensityFormula, outcomeFormula = outcomeFormula, constraint = constraint, boundaryPoint = boundaryPoint)
    fit$fit[[2]] <- rsInference(fit_tmp, efficient=efficient)
    # rotate the index
    tmp <- splitIndex$nuisance
    splitIndex$nuisance <- splitIndex$infer
    splitIndex$infer <- splitIndex$fit
    splitIndex$fit <- tmp

    fit_tmp <- rsfitSplit(covariate, response, treatment, splitIndex = splitIndex, lossType = lossType, propensityModel = propensityModel, estimatedPropensity = estimatedPropensity, outcomeModel=outcomeModel, estimatedOutcome = estimatedOutcome, weights = weights, tol=tol, propensityFormula = propensityFormula, outcomeFormula = outcomeFormula, constraint = constraint, boundaryPoint = boundaryPoint)
    fit$fit[[3]] <- rsInference(fit_tmp, efficient=efficient)
  } else {
    library(doParallel)
    n_cores <- detectCores(all.tests = FALSE, logical = TRUE)
    cl <- makeCluster(2)
    registerDoParallel(cl)
    fit <- foreach(iter = c(TRUE, FALSE))%dopar%{
      splitIndex_local <- NULL
      splitIndex_local$nuisance <- (splitIndex$nuisance==iter)
      splitIndex_local$fit <- (splitIndex$fit==iter)
      rsfitSplit(covariate, response, treatment, splitIndex = splitIndex_local, lossType = lossType, propensityModel = propensityModel, estimatedPropensity = estimatedPropensity, outcomeModel=outcomeModel, estimatedOutcome = estimatedOutcome, weights = weights, tol=tol, propensityFormula = propensityFormula, outcomeFormula = outcomeFormula, constraint = constraint, boundaryPoint = boundaryPoint)
    }
    stopCluster(cl)
  }
  fit$betaAN <- (fit$fit[[1]]$betaAN+fit$fit[[2]]$betaAN+fit$fit[[3]]$betaAN)/3
  W1 <- (fit$fit[[1]]$W1+fit$fit[[2]]$W1+fit$fit[[3]]$W1)/3
  W2 <- (fit$fit[[1]]$W2+fit$fit[[2]]$W2+fit$fit[[3]]$W2)/3
  fit$sigmaAN <- solve(W1) %*% W2 %*% solve(W1)
  fit
}






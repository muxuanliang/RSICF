# rsfitSolver obtains the optimizer of the propoed optimization given M_n
rsfitSolver <- function(covariate, response, treatment, estimatedNuisance, splitIndex = NULL, initLinkType = 'tanh', weights = NULL, tol = 10^(-3), m_0=0){
  # initialization
  fit_init <- fitLinkLinear(covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex, linkType = initLinkType, weights = weights, tol = tol)
  iter <- 0

  # iteration start
  fit_last <- fit_init
  diff <- 1
  while((diff > tol) & (iter < 1000)){
    # updateXi
    fit_xi <- updateXi(fit_last, covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex, weights = weights, tol = tol, m_0=m_0)
    # updateBeta
    fit_beta <- updateBeta(fit_xi, covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex, weights = weights, tol = tol)
    #update
    diff <- 1-sum(fit_beta$beta*fit_last$beta)/sqrt(sum(fit_beta$beta^2)*sum(fit_last$beta^2))
    iter <- iter+1
    fit_last <- fit_beta
  }
  # final update
  fit <- updateXi(fit_last, covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex, weights = weights, tol = tol, m_0=m_0)
  fit
}


# rsfitSplit obtain the estimation of RSICF assuming a single index model on a slited data. The variance estimaor of the coefficiencts are provided.
rsfitSplit <- function(covariate, response, treatment, splitIndex = NULL, propensityModel = 'glmnet', estimatedPropensity = NULL, outcomeModel = 'kernel', estimatedOutcome = NULL, initLinkType = 'tanh', weights = NULL, tol = 1e-3, propensityFormula=NULL, outcomeFormula=NULL){
  # fit nuisance parameter
  estimatedNuisance <- NULL
  estimatedNuisance$pi <- estimatedPropensity
  if (is.null(estimatedPropensity)){
    data <- list(predictor = covariate, treatment = treatment)
    predictedPropensityAll <- getPropensityModel(data, method = propensityModel, splitIndex = splitIndex$nuisance, Formula = propensityFormula, predictAll = TRUE)
    estimatedNuisance$pi <- predictedPropensityAll
  }
  estimatedNuisance$S <- estimatedOutcome$control+estimatedOutcome$treatment
  if (is.null(estimatedOutcome)){
    data <- list(predictor = covariate, treatment = treatment, outcome = outcome)
    predictedOutcomeAll <- getOutcomeModel(data, method = outcomeModel, splitIndex = splitIndex$nuisance, Formula = outcomeFormula, predictAll = TRUE)
    estimatedNuisance$S <- predictedOutcomeAll$control+predictedOutcomeAll$treatment
  }

  # start cross-validation for M_n
  initSolve <- rsfitSolver(covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex, initLinkType = initLinkType, weights =weights, tol = tol, withM=FALSE)
  m_max <- sum(abs(initSolve$xi))
  m_seq <- m_max * 10^(-seq(0,3, length.out = 10))
}

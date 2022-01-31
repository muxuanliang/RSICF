#' Estimation and inference for relative contrast function under a monotonic single-index model assumption
#'
#' @param covariate input matrix of dimension nobs x nvars. Each raw is a observation, each column is a covariate
#' @param response numeric response
#' @param treatment is a vector of binary value representing two treatment, 0 or 1.
#' @param estimatedPropensity estimated propensity score p(trt=1|X). This should be used only when the propensity is estimated by a parametric model.
#' @param estimatedOutcome estimated outcome model E(y|trt, X). This should be used only when the outcome is estimated by a parametric model. It should be a list (list(control=..., treatment=...)).
#' @param outcomeModel this selects method to estimate the outcome when estimatedOutcome=NULL. Options include lm, glmnet, kernel, and gam. If lm is used, user also need to input the outcomeFormula like y~x used in lm. By default, kernel regression is selected.
#' @param propensityModel Similar to outcomeModel.
#' @param boundaryPoint please specify appropriate range for the single-index. Default is [-1,1].
#' @return  A list
#' \describe{
#' \item{betaAN}{The one-step updated coefficients estimates}
#' \item{sigmaAN}{The estimated sd of the estimated coefficients}
#' \item{betaAN.renorm}{The one-step updated coefficients estimates after re-normalized such that \|\beta\|_2=1}
#' \item{sigmaAN.renorm}{The estimated sd of the normalized estimates}
#' \item{fit}{A list used for predict.rsfit which predicts the relative contrast}
#' }
#' @references Muxuan Liang, Menggang Yu (2022). Relative Contrast Estimation and Inference for Treatment Recommendation.
#'
#' @author Muxuan Liang
#' @export
#'
rsfit <- function(covariate, response, treatment, splitIndex = NULL, propensityModel = 'glmnet', estimatedPropensity = NULL, outcomeModel = 'kernel', estimatedOutcome = NULL, lossType = 'logistic', weights = NULL, tol = 1e-3, propensityFormula=NULL, outcomeFormula=NULL, parallel = FALSE, constraint = TRUE, boundaryPoint = c(-1,1), efficient = TRUE, local=TRUE){
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
    fit$fit[[1]] <- rsInference(fit_tmp, efficient=efficient, local = local)
    # rotate the index
    tmp <- splitIndex$nuisance
    splitIndex$nuisance <- splitIndex$infer
    splitIndex$infer <- splitIndex$fit
    splitIndex$fit <- tmp

    fit_tmp <- rsfitSplit(covariate, response, treatment, splitIndex = splitIndex, lossType = lossType, propensityModel = propensityModel, estimatedPropensity = estimatedPropensity, outcomeModel=outcomeModel, estimatedOutcome = estimatedOutcome, weights = weights, tol=tol, propensityFormula = propensityFormula, outcomeFormula = outcomeFormula, constraint = constraint, boundaryPoint = boundaryPoint)
    fit$fit[[2]] <- rsInference(fit_tmp, efficient=efficient, local = local)
    # rotate the index
    tmp <- splitIndex$nuisance
    splitIndex$nuisance <- splitIndex$infer
    splitIndex$infer <- splitIndex$fit
    splitIndex$fit <- tmp

    fit_tmp <- rsfitSplit(covariate, response, treatment, splitIndex = splitIndex, lossType = lossType, propensityModel = propensityModel, estimatedPropensity = estimatedPropensity, outcomeModel=outcomeModel, estimatedOutcome = estimatedOutcome, weights = weights, tol=tol, propensityFormula = propensityFormula, outcomeFormula = outcomeFormula, constraint = constraint, boundaryPoint = boundaryPoint)
    fit$fit[[3]] <- rsInference(fit_tmp, efficient=efficient, local = local)
  } else {
    library(doParallel)
    n_cores <- detectCores(all.tests = FALSE, logical = TRUE)
    cl <- makeCluster(3)
    registerDoParallel(cl)
    result <- foreach(iter = c(1,2,3), .packages = c("RSICF","ITRInference", 'glmnet', 'mgcv'))%dopar%{
      splitIndex_local <- NULL
      if (iter == 1){
        splitIndex_local <- splitIndex
      } else if (iter == 2) {
        splitIndex_local$nuisance <- splitIndex$infer
        splitIndex_local$infer <- splitIndex$fit
        splitIndex_local$fit <- splitIndex$nuisance
      } else {
        splitIndex_local$nuisance <- splitIndex$fit
        splitIndex_local$infer <- splitIndex$nuisance
        splitIndex_local$fit <- splitIndex$infer
      }
      fit_tmp <- rsfitSplit(covariate, response, treatment, splitIndex = splitIndex_local, lossType = lossType, propensityModel = propensityModel, estimatedPropensity = estimatedPropensity, outcomeModel=outcomeModel, estimatedOutcome = estimatedOutcome, weights = weights, tol=tol, propensityFormula = propensityFormula, outcomeFormula = outcomeFormula, constraint = constraint, boundaryPoint = boundaryPoint)
      rsInference(fit_tmp, efficient=efficient, local = local)
    }
    stopCluster(cl)
    fit$fit <- result
  }
  fit$betaAN <- (fit$fit[[1]]$betaAN+fit$fit[[2]]$betaAN+fit$fit[[3]]$betaAN)/3
  fit$sigmaAN <- (fit$fit[[1]]$sigmaAN+fit$fit[[2]]$sigmaAN+fit$fit[[3]]$sigmaAN)/3

  # normalize
  fit$betaAN.renorm <- (fit$fit[[1]]$betaAN.renorm+fit$fit[[2]]$betaAN.renorm+fit$fit[[3]]$betaAN.renorm)/3
  fit$sigmaAN.renorm <- (fit$fit[[1]]$sigmaAN.renorm+fit$fit[[2]]$sigmaAN.renorm+fit$fit[[3]]$sigmaAN.renorm)/3
  norm1 <- sqrt(sum(fit$betaAN.renorm^2))
  trans <- matrix(0, NROW(fit$sigmaAN.renorm), NCOL(fit$sigmaAN.renorm))
  for (i in 1: NROW(fit$sigmaAN.renorm)){
    for (j in 1:NCOL(fit$sigmaAN.renorm)){
      if (i==j){
        trans[i,j] <- (norm1^2-(fit$betaAN.renorm[i])^2)/(norm1^3)
      } else {
        trans[i,j] <- -(fit$betaAN.renorm[i]) * (fit$betaAN.renorm[j])/(norm1^3)
      }
    }
  }
  fit$sigmaAN.renorm <- trans %*% fit$sigmaAN.renorm %*% trans
  fit$betaAN.renorm <- fit$betaAN.renorm/norm1

  #return
  fit$sigmaAN <- diag(fit$sigmaAN)
  fit$sigmaAN.renorm <- diag(fit$sigmaAN.renorm)
  fit
}






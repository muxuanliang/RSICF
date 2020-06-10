rsInference <- function(rsfit, efficient = FALSE){
  # set parameter
  estimatedNuisance <- rsfit$estimatedNuisance
  splitIndex <- rsfit$splitIndex
  covariate <- rsfit$covariate
  response <- rsfit$response
  treatment <- rsfit$treatment
  lossType <- rsfit$lossType

  # set nuisance parameter
  pi <- estimatedNuisance$pi
  S <- estimatedNuisance$S

  # set weight
  weights <- 1
  if(efficient){
    eta <- predict.rsfitSplit(rsfit$fit, newx = covariate, type='eta')
    # augmentation
    ratio <- (1-eta)/(1+eta)
    inv_W <- 1/estimatedNuisance$pi * ratio + 1/(1-estimatedNuisance$pi) * 1/ratio
    weights <- inv_W^(-1) * estimatedNuisance$S
  }

  # set bias of x
  link <- covariate %*% rsfit$fit$beta
  exp_gam <- function(x, y){
    shadow_data <- data.frame(predictor = x, response = y)
    shadow_model <- gam(response~predictor, data = shadow_data)
    as.matrix(predict(shadow_model, newdata = data.frame(predictor = link),type='response'))
  }
  tmp_w <- response/(estimatedNuisance$pi*treatment+(1-estimatedNuisance$pi)*(1-treatment))
  ESX <- apply(covariate[splitIndex$fit,], 2, function(t){exp_gam(link[splitIndex$fit], weights[splitIndex$fit]*tmp_w[splitIndex$fit] * t)})
  ES <- exp_gam(link[splitIndex$fit], tmp_w[splitIndex$fit]*weights[splitIndex$fit])

  # start inference
  predict <- predict.rsfitSplit(rsfit$fit, newx = covariate)
  tildeX <- apply(ESX,2,function(t){t/ES})
  w0 <- predict.rsfitSplit(rsfit$fit, newx = covariate, derivative=TRUE)
  w1 <- response/(estimatedNuisance$pi*treatment+(1-estimatedNuisance$pi)*(1-treatment)) * loss(2*(treatment-0.5) * predict, type = lossType, order = 'hessian') * w0^2 * weights
  w2 <- (2*(treatment-0.5)*response/(estimatedNuisance$pi*treatment+(1-estimatedNuisance$pi)*(1-treatment)) * loss(2*(treatment-0.5) * predict, type = lossType, order = 'derivative')+2*(treatment-0.5)*(estimatedNuisance$S-response)/(2*(estimatedNuisance$pi*treatment+(1-estimatedNuisance$pi)*(1-treatment)))) * w0 * weights
  sampleIndex <- c(1:NROW(covariate))
  W1 <- matrix(0, NCOL(covariate), NCOL(covariate))
  W2 <- matrix(0, NCOL(covariate), NCOL(covariate))
  S <- matrix(0, NCOL(covariate), 1)
  for (iter in sampleIndex[splitIndex$fit]){
    W1 <- W1+w1[iter] * (covariate[iter,]) %*% t(covariate[iter,]-tildeX[iter,])
    W2 <- W2+(w2[iter])^2 * (covariate[iter,]-tildeX[iter,]) %*% t(covariate[iter,]-tildeX[iter,])
    S <- S+w2[iter]* (covariate[iter,]-tildeX[iter,])
  }
  W1 <- W1/sum(splitIndex$fit)
  #W1 <- W1 + 10^(-5)*diag(min(abs(diag(W1))),NCOL(covariate),NCOL(covariate))
  W2 <- W2/sum(splitIndex$fit)
  S <- S/sum(splitIndex$fit)
  betaAN <- rsfit$fit$beta-solve(W1) %*% S
  sigmaAN <- solve(W1) %*% W2 %*% solve(W1)
  # return
  list(fit=rsfit$fit, betaAN=betaAN, W1=W1, W2=W2, sigmaAN=sigmaAN)
}

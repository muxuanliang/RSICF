#predict.rsfit get the predicted treatment
predict.rsfit <- function(fit, newx, lossType = 'logistic', type = 'transformed'){
  splitNumber <- length(fit)
  res_all <- sapply(fit, function(t){
    z = newx %*% t$beta
    pre <- cbind(1,splines::bs(z, knots=t$knots[2:(length(t$knots)-1)], intercept = FALSE, Boundary.knots = c(t$knots[1], t$knots[length(t$knots)])))%*%t$xi
    if(sum(abs(t$xi))<1e-5){
      pre <- z
    }
    pre
  })
  if (type=='eta'){
    res_all <- link(res_all, lossType = lossType)
  }
  apply(res_all,1,mean)
}

# predict.rsfit gets the predicted treatment for a splitted fit
predict.rsfitSplit <- function(fit, newx, derivative = FALSE,lossType = 'logistic', type = 'value'){
  t <- fit
  z <- newx %*% t$beta
  pre <- cbind(1,splines::bs(z, knots=t$knots[2:(length(t$knots)-1)], intercept = FALSE, Boundary.knots = c(t$knots[1], t$knots[length(t$knots)])))%*%t$xi
  if(sum(abs(t$xi))<1e-5){
    pre <- z
  }
  if(derivative){
    pre <- cbind(0, nsd(z, knots=t$knots[2:(length(t$knots)-1)], intercept = FALSE, Boundary.knots = c(t$knots[1], t$knots[length(t$knots)])))%*%t$xi
    if(sum(abs(t$xi))<1e-5){
      pre <- 1
    }
  }
  if (type=='eta'){
    pre <- link(pre, lossType = lossType)
  }
  pre
}

#predict.rsfit get the predicted treatment
predict.rsfit <- function(fit, newx){
  splitNumber <- length(fit)
  res_all <- sapply(fit, function(t){
    z = newx %*% t$beta
    pre <- splines::ns(z, knots=t$knots[2:(length(t$knots)-1)], intercept = TRUE, Boundary.knots = c(t$knots[1], t$knots[length(t$knots)]))%*%t$xi
    if(sum(abs(t$xi))<1e-5){
      pre <- z
    }
    pre
  })
  apply(res_all,1,mean)
}

# predict.rsfit gets the predicted treatment for a splitted fit
predict.rsfitSplit <- function(fit, newx, derivative = FALSE){
  t <- fit
  z <- newx %*% t$beta
  pre <- splines::ns(z, knots=t$knots[2:(length(t$knots)-1)], intercept = TRUE, Boundary.knots = c(t$knots[1], t$knots[length(t$knots)]))%*%t$xi
  if(sum(abs(t$xi))<1e-5){
    pre <- z
  }
  if(derivative){
    pre <- nsd(z, knots=t$knots[2:(length(t$knots)-1)], intercept = TRUE, Boundary.knots = c(t$knots[1], t$knots[length(t$knots)]))%*%t$xi
    if(sum(abs(t$xi))<1e-5){
      pre <- 1
    }
  }
  pre
}

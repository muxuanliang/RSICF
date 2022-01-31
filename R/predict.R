#' predict.rsfit get the predicted relative contrast function
#'
#' @param fit the list fit in the output of rsfit.
#' @param newx input matrix of dimension nobs x nvars. Each raw is a new observation, each column is a covariate
#' @type Default is 'transformed' which gives the prediction of the  target relative contrast function in the paper; another option is 'eta', which gives the prediction of E[Y_1-Y_{-1}|X]/E[Y_1+Y_{-1}|X].
#' @return  A Array
#' \describe{
#' A array includes all predictions.
#' }
#' @references Muxuan Liang, Menggang Yu (2022). Relative Contrast Estimation and Inference for Treatment Recommendation.
#'
#' @author Muxuan Liang
#' @export
#'
predict.rsfit <- function(fit, newx, lossType = 'logistic', type = 'transformed'){
  splitNumber <- length(fit)
  res_all <- sapply(fit, function(t){
    z = newx %*% t$fit$beta
    pre <- cbind(1,splines2::bSpline(z, knots=t$fit$knots, intercept = FALSE, Boundary.knots = t$fit$boundaryPoint))%*%t$fit$xi
    if(sum(abs(t$fit$xi))<1e-5){
      pre <- z
    }
    pre
  })
  if (type=='eta'){
    res_all <- link(res_all, lossType = lossType)
  }
  apply(res_all,1,mean)
}


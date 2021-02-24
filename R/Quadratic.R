#' Quadratic
#' 
#' Quadratic dose-response curve, for MLE estimation. Used as an objective function for mle2
#' 
#' @param p 1x4 vector of parameters: b0 = intercept, b1 = shape parameter, b2 = X-value at the max/min, and s = standard deviation
#' @param Y concentration values for each observation
#' @param C response value for each observation (not used)
#' 
#' @return negative log-likelihood of the model fit
#' 
#' @export
#' 

Quadratic = function(p,X,Y){
  b0 = p[1] # maximum (or minimum)
  b1 = p[2] # shape parameter
  b2 = p[3] # X-value at the max/min
  s = exp(p[4])
  P = b0 + b1*(X-b2)^2 # equation of the line
  -sum(dnorm(x=Y,mean=P,sd=s,log=TRUE),na.rm=TRUE)}
parnames(Quadratic) = c('b0','b1','b2','s')
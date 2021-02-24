#' Int
#' 
#' Intercept-only dose-response curve, for MLE estimation. Used as an objective function for mle2
#' 
#' @param p 1x2 vector of parameters, g = mean and s = standard deviation
#' @param Y concentration values for each observation
#' @param C response value for each observation (not used)
#' 
#' @return negative log-likelihood of the model fit
#' 
#' @export


Int <- function(p,Y,C){
  g = p[1]
  s = exp(p[2])
  -sum(dnorm(x=Y,mean=g,sd=s,log=TRUE),na.rm=TRUE)}
parnames(Int) = c('g','s')
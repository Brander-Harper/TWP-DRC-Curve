#' SigmoidDR
#'
#' 4-parameter sigmoidal dose-response curve, for MLE estimation. Used as an objective function for mle2
#'
#' @param p 1x5 vector of parameters: g = max Y value, h = min Y value, a = exponential intercept term, b = exponential rate term, and s = standard deviation
#' @param Y concentration values for each observation
#' @param C response value for each observation (not used)
#'
#' @return negative log-likelihood of the model fit
#'
#' @export
#'

SigmoidDR = function(p,X,Y){
  g = p[1] # maximum
  h = p[2] # minimum
  a = p[3] # exponential intercept term
  b = p[4] # exponential rate term
  s = exp(p[5]) # standard deviation of response; use exp to avoid negative values
  P = (h-g)/(1+exp(a + b*X))+g # the equation for the curve
  # Assume a normal distribution for the likelihood (probably OK if data are log transformed), this is defining the actual function that you are fitting to the data

  -sum(dnorm(x=Y,mean=P,sd=s,log=TRUE),na.rm=TRUE)}

parnames(SigmoidDR) = c('g','h','a','b','s')

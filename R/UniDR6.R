#' UniDR6
#'
#' 6-parameter unimodal dose-response curve, for MLE estimation. Used as an objective function for mle2
#'
#' @param p 1x7 vector of parameters: g = max Y value, h = min Y value, a1 = exponential intercept term, b1 = exponential rate term (left side), a2 = exponential intercept term, b2 = exponential rate term (right side), and s = standard deviation
#' @param Y concentration values for each observation
#' @param C response value for each observation (not used)
#'
#' @return negative log-likelihood of the model fit
#'
#' @export
#'

uniDR6 = function(p,X,Y){

  g = p[1] # maximum
  h = p[2] # minimum
  a1 = p[3] # exponential intercept term
  b1 = abs(p[4]) # exponential rate term (rate of increase before the hump)
  a2 = p[5] # exponential intercept term
  b2 = -abs(p[6]) # exponential rate term (rate of decrease after the hump)
  s = exp(p[7]) # standard deviation of response (use exp to avoid negative values)
  P = g+h*(1 + exp(-(a1+b1*X)))/(1+exp(-(a2 + b2*X))) # the equation for the curve
  # Assume a normal distribution for the likelihood (probably OK if data are log transformed), this is defining the actual function that you are fitting to the data

  -sum(dnorm(x=Y,mean=P,sd=s,log=TRUE),na.rm=TRUE)}
parnames(uniDR6) = c('g','h','a1','b1','a2','b2','s')

#' Linear
#' 
#' Linear dose-response curve, for MLE estimation. Used as an objective function for mle2
#' 
#' @param p 1x3 vector of parameters, b0 = intercept, b1 = slope, and s = standard deviation
#' @param Y concentration values for each observation
#' @param C response value for each observation (not used)
#' 
#' @return negative log-likelihood of the model fit
#' 
#' @export
#' 

Linear = function(p,X,Y){
b0 = p[1]
b1 = p[2]
s = exp(p[3])
P = b0 + b1*X # equation of the line
-sum(dnorm(x=Y,mean=P,sd=s,log=TRUE),na.rm=TRUE)}
parnames(Linear) = c('b0','b1','s')
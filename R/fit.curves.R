#' fit.curves
#'
#' Fit dose-response curves using MLE
#'
#' @param X concentrations
#' @param Y response values
#' @param Fname filename for saving figures
#' @param Ylab Ylabel for figures
#' @param Xlab Xlabel for figures
#' @param axis.font fontsize for figures
#' @param log.factor value added to log-transformed data to avoid log(0)
#' @param FigW Figure width
#' @param FigH Figure height
#' @param SigInit 2-element vector with start values for sigmoidal function
#' @param Uni6Init 6-element vector with start values for 6-param unimodal
#' @param Uni6Init 5-element vector with start values for 5-param unimodal
#' @param reltol numerical tolerance for mle2
#'
#' @import bbmle
#'
#' @return baseline mle2 object for null model
#' @return sigmoidal mle2 object for sigmoid DR curve
#' @return linear mle2 object for linear DR curve
#' @return unimodal5 mle2 object for 5-param unimodal DR curve
#' @return unimodal6 mle2 object for 6-param unimodal DR curve
#' @return quadratic mle2 object for quadratic DR curve
#' @return AIC.table Table with AICc results for each model
#' @return ANOVA Table with likelihood ratio test results for each model
#' @return FN File name used for plotting
#' @return B Identity of 'best' model as judged by likelihood ratio test
#' @return Pval Table of likelihood ratio test p-values
#'
#' @export


fit.curves = function(X,Y,Fname = '',Ylab='Y',Xlab='X',axis.font = 0.75,log.factor = 0,
                      FigW=7,FigH=7,SigInit=c(0,0),Uni6Init=c(0,-1,2,2,2,2),
                      Uni5Init=NA,reltol=1e-4){

  # Fit a suite of models using mle2 numerical routine, using best guesses for initial conditions

  # Intercept-only
  mBase = mle2(Int,start=c(g = mean(Y),s=1),data=list(Y=Y))

  # Linear
  # Initial guess at slope using lm (yes this is silly because MLE = LS in this case, but this ensures all the results are a MLE2 object)
  mLin = mle2(Linear,start=c(b0 = summary(lm(Y~X))[[4]][[1]],
                             b1=summary(lm(Y~X))[[4]][[2]],
                             s=1),data=list(X=X,Y=Y))

  # Sigmoidal
  mSigm = mle2(SigmoidDR,start=c(g=mean(Y),h=mean(Y),a=SigInit[1],b=SigInit[2],s=1),data=list(X=X,Y=Y))

  # Unimodal (6-parameter)
  # Find good starting values:
  mUni6 = mle2(uniDR6,start=c(g=Uni6Init[1],h=Uni6Init[2],a1=Uni6Init[3],b1=Uni6Init[4],
                              a2=Uni6Init[5],b2=Uni6Init[6],s=1),data=list(X=X,Y=Y),
                              control=list(maxit=5000, reltol=reltol))

  # Unimodal (5-parameter)
  if (is.na(Uni5Init)){ # if it was not user-specified
    G = max(Y)-min(Y) # 0 # Difference between maximum and minimum values (i.e., y-range of data
    H = min(Y) #-1 # Minimum value
    # parameters for the increasing portion:
    B1 = 1 #2 # Rate of increase (must be > 1)
    # parameters for the decreasing portion
    A2 = 1 # Intercept (larger positive numbers move this to the right)
    B2 = -1 # Rate of decrease (must be < 1)
    }else{ # if user-specified
    G = Uni5Init[1]
    H= Uni5Init[2]
    B1= Uni5Init[3]
    A2== Uni5Init[4]
    B2== Uni5Init[5]}
   mUni5 = mle2(uniDR5,start=c(g=G,h=H,b1=B1,a2=A2,b2=B2,s=1),data=list(X=X,Y=Y),
               control=list(maxit=5000, reltol=reltol))

  # Unimodal (4-parameter, deprecated)
  #G = 0 # Difference between maximum and minimum values (i.e., y-range of data
  #H = -1 # Minimum value
  # parameters for the increasing portion:
  # B1 = 2 # Rate of increase (must be > 1)
  # parameters for the decreasing portion
  #  B2 = -4 # Rate of decrease (must be < 1)
  # mUni4 = mle2(uniDR4,start=c(g=G,h=H,b1=B1,b2=B2,s=1),data=list(X=X,Y=Y))

  # Quadratic
  mQuad = mle2(Quadratic,start=c(b0 = mean(Y),b1=1,b2=mean(X),s=1),data=list(X=X,Y=Y),
               control=list(maxit=5000, reltol=reltol))

  # Calculate AICs from mle objects
  AIC.table = AICctab(mBase,mLin,mSigm,mUni6,mUni5,mQuad,nobs=length(X),sort=F)

  # Likelihood ratio tests.
  # Sometimes higher-order models fail to converge, violating the assumptions of the LRT (that the null model always has lower likelihood)
  # A constraint is applied to avoid misleading p-values in that case
  ANOVA = list()
  ANOVA[[1]] = if (logLik(mBase) < logLik(mLin)){anova(mBase,mLin)}else{1}
  ANOVA[[2]] = if (logLik(mBase) < logLik(mQuad)){anova(mBase,mQuad)}else{1}
  ANOVA[[3]] = if (logLik(mBase) < logLik(mSigm)){anova(mBase,mSigm)}else{1}
  ANOVA[[4]] = if (logLik(mBase) < logLik(mUni5)){anova(mBase,mUni5)}else{1}
  ANOVA[[5]] = if (logLik(mBase) < logLik(mUni6)){anova(mBase,mUni6)}else{1}
  # ANOVA.6 = anova(mBase,mUni3)


  # Tally up p-values and deviances
  Pval = rep(1,length(ANOVA)) # number of replicates
  Dev = rep(NA,length(ANOVA)) # number of replicates
  for (a in 1:length(ANOVA)){
    if (length(ANOVA[[a]])>1){
      Pval[a]=ANOVA[[a]][10] # pvalue
      Dev[a]=ANOVA[[a]][4] # deviance (-2 * log likelihood)
    }}

  # Pick the best one and make a plot, based on Pvals
  # Note: it might be better to do this based on AIC...
  Best = match(min(Pval),Pval)
  # Choose best model based on AIC:
  #Best = match(min(AIC.table$dAICc),AIC.table$dAICc)

  Title = Fname #paste(Cycle,' cycle')
  Xticks = unique(Db$Treatment)
  Xticks.at = log10(Xticks + log.factor)
  # some plotting labeling settings suppressed for now, to allow easy plotting of
  # wider data range

  # This opens a new graphic window of a specific size, to be saved as a pdf file
  # You can adjust these
  switch(Sys.info()[['sysname']],
         Windows= {windows(width=FigW,height=FigH)},
         Darwin = {quartz(width=FigW,height=FigH)},
         Linux = {x11(width=FigW,height=FigH)})

  plot(X,Y, #xaxp=c(-3,3,12),
       ylab=Ylab, xlab = Xlab, main = Title, # X & Y labels and Title
       cex = 1, pch = 1, # size & symbol type for the markers
       cex.lab = 1, cex.axis = axis.font, # size for axes labels and tick marks
       xaxt = 'n') # turns off default xtick labels

  axis(side = 1,at=Xticks.at,labels=Xticks,cex.axis=axis.font) # plot xtick lables where the actual concentrations are

  X_temp = seq(-5,5,length.out=1000) # dummy x-values for plotting the curve

  # Choose the correct line type
  # the code below is appropriate if basing results on an AIC table.
  #if (Best > 1){ # Best == 1 corresponds to intercept-only model which should always be dashed
  #if (Pval[Best-1]<0.05){Lty = 1}else{Lty=2} }
  # the code below is appropriate if basing results on the LRT ('ANOVA') table
  if (Pval[Best]<0.05){Lty = 1}else{Lty=2}


  # Note that the numbers corresponding to each model depend on whether you use AIC or LRT to determine which ones to plot.
  # Currently configured for LRT
  # if (Best==1){ # Baseline
  #  lines(c(X_temp[1],X_temp[1000]),c(mean(Y),mean(Y)),lty=2)
  #  Resid = Y-mean(Y)
  #  File.name = paste("dose_response_figures/",Fname,"_",Cycle,"_baseline_dose_response.pdf",sep="")
  #  dev.print(device=pdf,file=File.name,useDingbats=FALSE)
  #}
  if (Best==1){ # Linear
    B = coef(mLin)
    P = B[1]+B[2]*X_temp
    Pr = B[1]+B[2]*X
    Resid = Y-X
    lines(X_temp,P,lty=Lty)
    File.name = paste("dose_response_figures/",Fname,"_linear_dose_response.pdf",sep="")
    dev.print(device=pdf,file=File.name,useDingbats=FALSE)
  }
  if (Best==2){ # Quadratic
    B = coef(mQuad)
    P = B[1]+B[2]*(X_temp-B[3])^2
    Pr = B[1]+B[2]*(X-B[3])^2
    Resid = Y-X
    lines(X_temp,P,lty=Lty)
    File.name = paste("dose_response_figures/",Fname,"_quadratic_dose_response.pdf",sep="")
    dev.print(device=pdf,file=File.name,useDingbats=FALSE)
  }
  if (Best==3){ # Sigmoidal
    B = coef(mSigm)
    P = (B[2]-B[1])/(1+exp(B[3] + B[4]*X_temp))+B[1]
    Pr = (B[2]-B[1])/(1+exp(B[3] + B[4]*X))+B[1]
    Resid = Y-X
    lines(X_temp,P,lty=Lty)
    File.name = paste("dose_response_figures/",Fname,"_sigmoidal_dose_response.pdf",sep="")
    dev.print(device=pdf,file=File.name,useDingbats=FALSE)
  }
  if (Best==4){ # Uni6
    B = coef(mUni6)
    P = B[1]+B[2]*(1 + exp(-(B[3]+B[4]*X_temp)))/(1+exp(-(B[5] + B[6]*X_temp)))
    Pr = B[1]+B[2]*(1 + exp(-(B[3]+B[4]*X)))/(1+exp(-(B[5] + B[6]*X)))
    Resid = Y-X
    lines(X_temp,P,lty=Lty)
    File.name = paste("dose_response_figures/",Fname,"_unimodal6_dose_response.pdf",sep="")
    dev.print(device=pdf,file=File.name,useDingbats=FALSE)
  }
  if (Best==5){ # Uni5
    B = coef(mUni5)
    P = B[1]+B[2]*(1 + B[3]*X_temp)/(1+exp(-(B[4] + B[5]*X_temp)))
    Pr = B[1]+B[2]*(1 + B[3]*X)/(1+exp(-(B[4] + B[5]*X)))
    Resid = Y-X
    lines(X_temp,P,lty=Lty)
    File.name = paste("dose_response_figures/",Fname,"_unimodal5_dose_response.pdf",sep="")
    dev.print(device=pdf,file=File.name,useDingbats=FALSE)
    dev.off()
  }



  # Also create normal QQ plot
  # Collect residuals...
  # This opens a new graphic window of a specific size, to be saved as a pdf file
  # You can adjust these
  switch(Sys.info()[['sysname']],
         Windows= {windows(width=FigW,height=FigH)},
         Darwin = {quartz(width=FigW,height=FigH)},
         Linux = {x11(width=FigW,height=FigH)})
  qqnorm(y=Resid)
  qqline(y=Resid)
  File.name = paste("dose_response_figures/",Fname,"_bestmodel_qqplot.pdf",sep="")
  dev.print(device=pdf,file=File.name,useDingbats=FALSE)
  dev.off()

  return(list("baseline"=mBase,"sigmoidal"=mSigm,"linear"=mLin,
              "unimodal5"=mUni5,"unimodal6"=mUni6,"quadratic"=mQuad,
              "AIC"=AIC.table,"ANOVA"=ANOVA,'FN'=File.name,'B'=Best,
              "Pval"=Pval))
}
#parnames(fit.curves) = c('Db','Cycle')

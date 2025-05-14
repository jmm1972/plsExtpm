#########################################################################
# plsExtcheck()
# Runs gam.check(mgcv function to check on the basis dimensions, which 
# can help to flag up terms for which k is too low. Grossly too small k will
# also be visible from partial residuals available with plot.gam.
#########################################################################
#########################################################################
plsExtcheck <- function(pls.object, plot="n"){
  #
  if(!require("mgcv")) stop("package 'mgcv' not available")
  #
  if(class(pls.object) != "plsExtpm")
    stop(pls.object, " is not of class 'plsExtpm'.")
  strucmod = pls.object$model$strucmod
  innerW.rlshp = pls.object$model$innerW.rlshp
  factor_scores = pls.object$factor_scores
  endog = pls.object$model$endogenous
  ##
  for(i in endog){
    if(all(strucmod[,3]=="ln")) stop("No nonlinear relationshipt in inner model.\n")
      if(any(innerW.rlshp[[i]][["type"]]=="nln")){
      #if(strucmod[which(strucmod[,2]==i),3]=="nln"){
        cat("Relastionship: ",innerW.rlshp[[i]][["formula"]],"\n")
        indpnt <- factor_scores[,innerW.rlshp[[i]][["predictors"]]]
        dpnt <- factor_scores[,i]
        mod.data = data.frame(dpnt,indpnt); names(mod.data) = c(i,innerW.rlshp[[i]][["predictors"]])
        gam.check(gam(formula(innerW.rlshp[[i]][["formula"]]),data=mod.data))
        cat ("Next relationship: Press [enter] to continue")
        line <- readline()
      }
  }
}
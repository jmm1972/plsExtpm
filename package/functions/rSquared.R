#########################################################################
# rSquared
# Method for accessing the R-squares for the submodels of a path model fitted by 'sempls'
#########################################################################
rSquared <- function(object){
  latent <- object$model$latent             # names of latent variables
  fscores <- object$factor_scores           # matrix of estimated factor scores
  strucmod <- object$model$strucmod           # matrix of estimated factor scores
  
  R.squared<- matrix(NA, nrow=length(latent), ncol=1)
  rownames(R.squared) <- latent
  colnames(R.squared) <- c("R-Squared")
  
  for(i in latent){
    if (i %in% strucmod[,2]){
      R.squared[i,1] = cor(object$path_coefficients[["predicted"]][,i],fscores[,i])
    }
  }
  return(R.squared)
}
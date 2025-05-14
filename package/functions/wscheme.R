#########################################################################
# centroid()
#' A function for the centroid fitting scheme for the inner weights within sempls. The centroid scheme only uses the sign of the correlation between latent variables connected with each other.
centroid <- function(model, fscores){
  D <- model$D
  C <- D + t(D)
  innerW <- C
  innerW[C == 1] <- sign(cor(fscores))[C == 1]
  return(innerW)
}
#########################################################################
# factorial()
# factorial weighting scheme
factorial <- function(model, fscores){
  D <- model$D
  C <- D + t(D)
  innerW <- C
  innerW[C == 1] <- cor(fscores, use="everything")[C == 1]
  return(innerW)
}
#########################################################################
# pathWeighting()
# path weighting scheme
pathWeighting <- function(model, fscores){
  D <- model$D
  latent = model$latent
  E <- D - t(D)
  pred <- predecessors(model)
  # calculating the inner weights
  innerW <- E
  for (i in latent)
  {
    if(length(pred[[i]])==0) next
    else if (length(pred[[i]])==1){
      innerW[pred[[i]], i] <- cor(fscores[,pred[[i]]], fscores[,i])    }
      innerW[pred[[i]], i] <- solve(cor(as.matrix(fscores[,pred[[i]]]))) %*% cor(fscores.sc[,pred[[i]]], fscores[,i])
    
  }
  
  innerW[E == 0] <- 0
  innerW[E == -1] <- cor(as.matrix(fscores[, latent]))[E == -1]
  # return the matrix of inner weights
  return(innerW)
}
#########################################################################
# smoothWeighting()
# smooth weighting scheme for smoothing PLS
smoothWeighting <- function(model, fscores){
  innerW.rnames = model$innerW.rlshp[["innerW.rnames"]]
  innerW.formula = model$innerW.rlshp[["formula"]]
  innerW.rlshp = model$innerW.rlshp
  latent = model$latent
  endog = model$endogenous
  innerW.exosucc = model$innerW.exosucc
  innerW = matrix(0,ncol=length(model$latent),nrow=length(innerW.rnames))
  colnames(innerW) = latent
  rownames(innerW) = innerW.rnames
  #
  for(i in endog){
    indpnt <- fscores[,innerW.rlshp[[i]][["predictors"]]]
    dpnt <- fscores[,i]
    mod.data=data.frame(dpnt,indpnt); names(mod.data) = c(i,innerW.rlshp[[i]][["predictors"]])
    cf = coef(gam(formula(innerW.rlshp[[i]][["formula"]]),data=mod.data))
    innerW[names(cf[-1]),i] <- cf[-1]
  }
  for (i in names(innerW.exosucc)){
    succ <- fscores[,innerW.exosucc[[i]][["successors"]]]
    exog <- matrix(scale(fscores[,i]), ncol=1); colnames(exog)=i
    mod.data = data.frame(scale(data.frame(exog,succ))); names(mod.data) = c(i,innerW.exosucc[[i]][["successors"]])
    mod.data = model.matrix(gam(formula(innerW.exosucc[[i]][["formula"]]),data=mod.data))[,-1]
    #if(length(setdiff(innerW.exosucc[[i]][["innerW.rnames"]],colnames(mod.data))!=0){mod.data = cbind(mod.data)
    innerW[innerW.exosucc[[i]][["innerW.rnames"]],i] <- cor(mod.data[,innerW.exosucc[[i]][["innerW.rnames"]]], exog[,i])
  }
return(innerW)
}

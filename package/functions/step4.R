#########################################################################
# step4()
# Outer estimation of the factor scores
step4 <- function(data, outerW, model){
  latent= model$latent
  innerW.rnames = model$innerW.rlshp[["innerW.rnames"]]
  innerW.formula = model$innerW.rlshp[["formula"]]
  innerW.rlshp = model$innerW.rlshp
  
  src.Latent <- as.matrix(data) %*% outerW
  colnames(src.Latent) <- model$latent
  #print(colnames(data))
  #print(colnames(Latent))
  #cat ("Step 4: Press [enter] to continue")
  #line <- readline()
  mod.scores = src.Latent
  for(i in latent){
    if(i %in% names(innerW.rlshp)){
      indpnt <- mod.scores[,innerW.rlshp[[i]][["predictors"]]]
      dpnt <- mod.scores[,i]
      mod.data = data.frame(dpnt,indpnt); names(mod.data) = c(i,innerW.rlshp[[i]][["predictors"]])
      #mod.scores = cbind(mod.scores,model.matrix(gam(formula(innerW.rlshp[[i]][["formula"]]),data=mod.data)))
      mod.scores = cbind(mod.scores,model.matrix(gam(formula(innerW.rlshp[[i]][["formula"]]),data=mod.data))[,-1])
    }
  }
  Latent = scale(mod.scores[,innerW.rnames])
  
  return(Latent=Latent)
}

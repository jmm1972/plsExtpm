###################################################################################################  
# step1 - initialisation
# Initial estimation of outer proxies
###################################################################################################  
step1 <- function(model, data, sum1){
  M <- model$M
  latent= model$latent
  innerW.rnames = model$innerW.rlshp[["innerW.rnames"]]
  innerW.formula = model$innerW.rlshp[["formula"]]
  innerW.rlshp = model$innerW.rlshp
  ##
  if(sum1){M <- apply(model$M, 2, function(x){x <- x/sum(x)})}
  ##
  src.latent = as.matrix(data) %*% M
  ##
  mod.scores = src.latent
  for(i in latent){
    if(i %in% names(innerW.rlshp)){
      indpnt <- mod.scores[,innerW.rlshp[[i]][["predictors"]]]
      dpnt <- mod.scores[,i]
      mod.data = data.frame(dpnt,indpnt); names(mod.data) = c(i,innerW.rlshp[[i]][["predictors"]])
      mod.scores = cbind(mod.scores,model.matrix(gam(formula(innerW.rlshp[[i]][["formula"]]),data=mod.data)))
    }
  }
  Latent = scale(mod.scores[,innerW.rnames])
  return(list(Latent=Latent,outerW=M))
}
###################################################################################################  
# END OF step1()
###################################################################################################  
#########################################################################
# Calculates the total effects
# total effects for nonlinear modelling is not correct
totEffects <- function(model,pCoeff,factor_scores){

  latent = model$latent
  strucmod <- model$strucmod                        # names of manifest variables
  fscores <- factor_scores           # matrix of estimated factor scores
  
  ret <- pCoeff
  
  step <- pCoeff
  class(pCoeff)
  
  for (i in 2:ncol(pCoeff)){
    step <- step %*% pCoeff
    ret <- step + ret
  }
  totEff<-ret
  
  indEff=totEff-pCoeff
  
  total.effects = data.frame(relationships = 0,
                             directPLS = 0,
                             indirectPLS = 0,
                             totalPLS = 0)
  #
  r=1
  #
  for (i in latent){
    for (j in setdiff(latent,i)){
      total.effects[r,"relationships"]=paste0(i, " -> ", j)
      total.effects[r,"directPLS"]=pCoeff[i,j]
      total.effects[r,"indirectPLS"]=indEff[i,j]
      total.effects[r,"totalPLS"] = totEff[i,j]
      #print(paste0(i, " -> ", j))
      r = r + 1
    }
  } 
  
  total.effects = total.effects[which(rowSums(total.effects[,-1])!=0),]
  
  return(list(totEff=totEff,totEff.summary=total.effects))
}

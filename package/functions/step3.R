#########################################################################
# step3()
# Outer Approximation (step 3 within PLS-Algorithm)
step3 <- function(Latent, data, model, sum1){
  blocks <- model$blocks
  W <- model$M
  latent=model$latent
  for (i in latent){
    if(length(blocks[[i]])==1) next
    mf <- as.matrix(subset(data, select=blocks[[i]]))
    fscores <- as.matrix(Latent[,i])

    # only for Mode "B": transform the MVs of a block
    if (attr(blocks[[i]], "mode") == "B") {
      S <- cor(mf, mf, use="everything")
      T <- solve(chol(S))
      mf <-  t(t(T) %*% t(mf)) # ((X'X)^{-1})X'
    }
    
    # the same for mode "A" and "B"Y'X
    W[blocks[[i]],i] <- cor(fscores, mf, use="everything") # 
    
    # only for Mode "B": retransform the weights according to the MVs
    if (attr(blocks[[i]], "mode") == "B") {
      W[blocks[[i]],i] <- T %*% W[blocks[[i]],i] #((X'X)^{-1})X'(X'Y)
    }
  }
  
  ## Normalize weights to colwise sum up to 1?
  if(sum1){W <- apply(W, 2, function(x){x <- x/sum(x)})}
  return(W=W)
}

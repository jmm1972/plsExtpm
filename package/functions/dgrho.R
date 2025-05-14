#########################################################################
# dgrho()
dgrho <- function(object){
  dgr <- matrix(NA, nrow=length(object$model$latent), ncol=2)
  rownames(dgr) <- object$model$latent
  colnames(dgr) <- c("Dillon-Goldstein's rho", "reflective MVs")
  for(i in object$model$latent){
    if(attr(object$model$blocks[[i]], "mode")=="B"){
      next
    }
    x <- object$outer_loadings[, i]
    ind <- which(x!=0)
    if(length(ind)==1){
      dgr[i,2] <- 1
      next
    }
    else {
      x <- x[ind]
      dgr[i,1] <- sum(x)^2 / (sum(x)^2 + sum(1-x^2))
      dgr[i,2] <- length(ind)
    }
  }
  return(dgr)
}

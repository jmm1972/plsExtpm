###########################################################################################
# gof()
# requires: rSquared, communality
gof <- function(object){
  rSq <- object$R.Squared
  ## aveRsq <- nrow(rSq)^-1 * sum(rSq[,1], na.rm=TRUE)
  aveRsq <- mean(rSq[,1], na.rm=TRUE)
  com <- communality(object)
  aveCom <- sum(com[!is.na(com[,1]), 2], na.rm=TRUE)^-1 *
    sum(com[,1] * com[,2], na.rm=TRUE)
  gof <- sqrt(aveCom * aveRsq)
  gof <- matrix(c(aveRsq, aveCom, gof), nrow=3, ncol=1)
  rownames(gof) <- c("Average R-squared", "Average Communality", "GoF")
  colnames(gof) <- c("Value")
  class(gof) <- c("gof", class(gof))
  return(gof)
} 

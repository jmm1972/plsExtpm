#########################################################################
# step2()
# Inner estimation of factor scores
step2 <- function(Lat, innerW, model){
  Latent<-Lat
  Latent <- scale(Latent %*% innerW)
  # the attributes for the scale are meaningless
  attributes(Latent)[c(3,4)] <- NULL
  return(Latent=Latent)
}

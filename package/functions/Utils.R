#########################################################################
# miscellabeous functions
# for object of class: plsm

str.mode = function(object){
  mode=c()
  for (i in object$model$latent)
  {
    mode[i]=attr(object$model$blocks[[i]], "mode")
  }
  return(mode)
}

str.type = function(object){
  type=c()
  for (i in object$model$latent)
  {
    mode[i]=attr(object$model$blocks[[i]], "mode")
  }
  return(mode)
}

exogenous <- function(object){
    ret <- setdiff(object$model$strucmod[,1],object$model$strucmod[,2])
    return(ret)
}

endogenous <- function(object){
    ret <- unique(object$model$strucmod[,2])
    return(ret)
}

formative <- function(object){
    ret <- names(which(lapply(object$model$blocks, function(x){attr(x, "mode")})=="B"))
    return(ret)
}

reflective <- function(object){
    ret <- names(which(lapply(object$model$blocks, function(x){attr(x, "mode")})=="A"))
    return(ret)
}

indicators <- function(object, LV){
    if(!LV %in% object$model$latent) stop("The LV must be contained in the model!")
    ret <- object$model$blocks[[LV]]
    return(ret)
}

noindicators <- function(object){
  noind=c()
  for (i in object$model$latent)
  {
    noind[i]=length(object$model$blocks[[i]])
  }
  return(noind)
}

# used in 'pathWeighting'
predecessors <- function(model){
    D <- model$D
    foo <- function(x) names(which(x==1))
    pred <- apply(D, 2, foo)
    return(pred)
}

successors <- function(model){
    D <- model$D
    foo <- function(x) names(which(x==1))
    succ <- apply(D, 1, foo)
    return(succ)
}
#########################################################################
# requires: outer loadings (factor scores), model
communality <- function(object, ...){
    com <- matrix(NA, nrow=length(object$model$latent), ncol=2)
    rownames(com) <- object$model$latent
    colnames(com) <- c("communality", "reflective MVs")
    for(i in object$model$latent){
        if(attr(object$model$blocks[[i]], "mode")=="B"){
            next
        }
        x <- object$outer_loadings[, i]
        ind <- which(x!=0)
        if(length(ind)==1){
            com[i,2] <- 1
            next
        }
        else {
            x <- x[ind]
            com[i,1] <- 1/length(x)*sum(x^2)
            com[i,2] <- length(ind)
        }
    }
    return(com)
}
#########################################################################
# function to evaluate different values for
# the omission distance 'd'
ommissionTest <- function(object, drange){
  options(digits=2)
  omt <- list()
  omt$PLS<-matrix(NA, nrow=length(endogenous(object$model)),ncol=length(drange))
  omt$PLSc<-matrix(NA, nrow=length(endogenous(object$model)),ncol=length(drange))
  omt$PLSs<-matrix(NA, nrow=length(endogenous(object$model)),ncol=length(drange))
  rownames(omt$PLS) <- endogenous(object$model)
  colnames(omt$PLS) <- as.character(drange)
  rownames(omt$PLSc) <- endogenous(object$model)
  colnames(omt$PLSc) <- as.character(drange)
  rownames(omt$PLSs) <- endogenous(object$model)
  colnames(omt$PLSs) <- as.character(drange)
  for(i in 1:length(drange)){
    qSqrd <- try(my.qSquared(object, d=drange[i]),silent=TRUE)
    if(!inherits(qSqrd, "try-error")){
      omt$PLS[,i] <- qSqrd[,1]
      omt$PLSc[,i] <- qSqrd[,2]
      omt$PLSs[,i] <- qSqrd[,3]
    }
  }
  return(omt)
}


#### see
#
sse <- function(model, data, dblind, i, impfun){
  #plsm <- sempls(model, impfun(dblind), pairwise=TRUE, verbose=FALSE, ...)
  res <- plsextpm(model, impfun(dblind), verbose=FALSE)
  #model=model;data=impfun(dblind);pairwise=TRUE;verbose=FALSE
  
  m <- attr(res$data, "scaled:center")[model$blocks[[i]]]
  s <- attr(res$data, "scaled:scale")[model$blocks[[i]]]
  sse <- vector("numeric", length=length(m))
  ssec <- vector("numeric", length=length(m))
  sses <- vector("numeric", length=length(m))
  for(l in 1:length(m))
  {
    # Estimation
    yind <- which(i==model$latent)
    ind <- is.na(dblind[, model$blocks[[i]][l]])
    indf <- complete.cases(res$factor_scores[, -yind])
    ind <- as.logical(ind*indf)
    e <- res$factor_scores[ind, -yind, drop=FALSE] %*%
      res$path_coefficients$pC[-yind, i, drop=FALSE] *
      res$outer_loadings[model$blocks[[i]][l],i]
    ec <- res$factor_scores[ind, -yind, drop=FALSE] %*%
      res$path_coefficients$pCc[-yind, i, drop=FALSE] *
      res$PLSc_outer_loadings[model$blocks[[i]][l],i]
    es <- res$path_coefficients$spl.pC[[i]]$f[ind,"Sum"] *
      res$outer_loadings.spl[model$blocks[[i]][l],i]
    # Rescaling
    e <- e * s[l] + m[l]
    ec <- ec * s[l] + m[l]
    es <- es * s[l] + m[l]
    sse[l] <- sum((data[ind, model$blocks[[i]][l]] - e)^2, na.rm=TRUE)
    ssec[l] <- sum((data[ind, model$blocks[[i]][l]] - ec)^2, na.rm=TRUE)
    sses[l] <- sum((data[ind, model$blocks[[i]][l]] - es)^2, na.rm=TRUE)
  }
  return(list(sum(sse),sum(ssec),sum(sses)))
}
#
#
sso <- function(model, data, dblind, i){
  m <- try(apply(dblind[, model$blocks[[i]]], 2, mean, na.rm=TRUE), silent=TRUE)
  if(inherits(m, "try-error")) m <- mean(dblind[, model$blocks[[i]]], na.rm=TRUE)
  sso <- vector("numeric", length=length(m))
  for(k in 1:length(m)){
    ind <- is.na(dblind[, model$blocks[[i]][k]])
    sso[k] <- sum((data[ind, model$blocks[[i]][k]] - m[k])^2, na.rm=TRUE)
  }
  return(sum(sso))
}

#######################################################
################ RECYCLE BIN ##########################
#######################################################
#########################################################################
# rhoA()
# Calculates the vector of attenuation factors.
rhoA <- function(model, data, weights){
  latent <- model$latent             # names of latent variables
  strucmod <- model$strucmod         # 
  measuremod <- model$measuremod     #
  data<-data				 # data: manifest variables
  weights<-weights			 # outer weights matrix
  
  rhoA<- numeric(length(latent))
  names(rhoA) <- latent
  
  for(i in latent){
    index <- which(i == measuremod[,1])
    S=cor(data[,index], use="everything")
    dS=matrix(0,nrow=nrow(S),ncol=ncol(S))
    diag(dS)=diag(S)
    w=weights[index,i]
    dw = matrix(0,nrow=length(w),ncol=length(w))
    diag(dw) = diag(w%*%t(w))
    rhoA[i]=(t(w)%*%w)^2*(t(w)%*%(S-dS)%*%w)/(t(w)%*%(w%*%t(w)-dw)%*%w)
    if(attr(model$blocks[[1]], "mode")=="B"){rhoA[i,1]=1}
  }
  return(rhoA)
}
###########################################################################################################################
### deriv(): computes derivatives of inner relationships for PLS
###########################################################################################################################
plsderiv <- function(x, y) {
  if (length(x) != length(y)) {
    stop('x and y vectors must have equal length')
  }
  
  n <- length(x)
  
  # Initialize a vector of length n to enter the derivative approximations
  fdx <- vector(length = n)
  
  # Iterate through the values using the forward differencing method
  for (i in 2:n) {
    fdx[i-1] <- (y[i-1] - y[i]) / (x[i-1] - x[i])
  }
  
  # For the last value, since we are unable to perform the forward differencing method 
  # as only the first n values are known, we use the backward differencing approach
  # instead. Note this will essentially give the same value as the last iteration 
  # in the forward differencing method, but it is used as an approximation as we 
  # don't have any more information
  fdx[n] <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])
  
  return(fdx)
}
###########################################################################################################################
### END OF plsderiv(): computes derivatives of inner relationships for PLS
###########################################################################################################################

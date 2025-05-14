###############################################################################
# plsExtboot()
plsExtboot <- function(pls.object, nboot=200, start=c("ones", "old"), conf.level = 0.95, verbose=TRUE){
    result <- list()
    #
    if(!(start %in% c("old","ones"))) stop("Starting outer weights must be 'old' or 'ones'.")
    #
    if(class(pls.object) != "plsExtpm")
      stop(pls.object, " is not of class 'plsExtpm'.")
    #
    class(result) = "plsExtboot"
    #
    path.rlshp = list()
    #
    latent = pls.object$model$latent
    scaled = pls.object$scaled
    nrlshp = sum(sapply(pls.object$path_coefficients$impact,dim)[2,])
    ngrid = pls.object$ngrid
    impact = pls.object$path_coefficients$impact
    endogenous = pls.object$model$endogenous
    data <- pls.object$model$src.data
    fscores = pls.object$factor_scores
    attr(fscores,"scaled:center") = attr(pls.object$factor_scores,"scaled:center")
    attr(fscores,"scaled:scale") = attr(pls.object$factor_scores,"scaled:scale")
    #
    if (!require("boot")) stop("package boot not available")
    refit <- function(){
        data <- data[indices,]
        refitted_model <- replsExtpm(plsextpm=pls.object, data = data, start=start)
        refitted_model
    }
    # the following 2 lines borrowed from boot in package boot
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    warn <- options(warn=-2)
    on.exit(options(warn)) # insure restore even in event of error
    nErrors <- 0
    N <- nrow(data)
    tryErrorIndices <- NULL
    bootIndices <- matrix(NA, ncol=nboot, nrow=nrow(data))
    for (b in 1:nboot){bootIndices[,b] <- sample(N, N, replace=TRUE)}
    cat("Replicate: ")
    for (b in 1:nboot){
        indices <- bootIndices[,b]
        treil <- paste(rep(" ", floor(log10(nboot)) - floor(log10(b))), collapse="")
        ndel <- paste(rep("\b", floor(log10(nboot)) + 1), collapse="")
        if(b==1) cat(paste(treil, b, sep=""))
        if(b!=nboot) cat(paste(ndel, treil, b, sep="")) else cat(paste(ndel, b, " Done.\n", sep=""))
        for (try in 1:11){
        	if (try > 10){
			      res <- replsExtpm(plsextpm=pls.object, data = data, start)
			      cat(paste(" Replicate : ", b, " failed. Replaced by original data.\n",sep=""))
    			  cat("Replicate: ")
		  	    result[[b]] = list(inner_weights=res$inner_weights,
		  	                       factor_scores=res$factor_scores,
		  	                       path_coefficients=res$path_coefficients,
		  	                       cross_loadings=res$cross_loadings,
		  	                       outer_loadings=res$outer_loadings,
		  	                       outer_weights=res$outer_weights,
		  	                       iterations=res$iterations,
		  	                       bootIndices=1:N,
		  	                       seed=seed,
		  	                       nboot=nboot,
		  	                       data=data[indices,])
			      break()
		      }
        	res <- try(refit(), silent=TRUE)
        	if(!inherits(res, "try-error")){
		  	    result[[b]] = list(inner_weights=res$inner_weights,
		  	                       factor_scores=res$factor_scores,
		  	                       path_coefficients=res$path_coefficients,
		  	                       cross_loadings=res$cross_loadings,
		  	                       outer_loadings=res$outer_loadings,
		  	                       outer_weights=res$outer_weights,
		  	                       iterations=res$iterations,
		  	                       bootIndices=1:N,
		  	                       seed=seed,
		  	                       nboot=nboot,
		  	                       data=data[indices,])
      			break()				  
		      }
		    indices <- sample(N, N, replace=TRUE)
		    bootIndices[,b] <- indices
       }
    }
    options(warn)
    if (nErrors > 0) warning("There were ", nErrors,
                             " apparent convergence failures;\n",
                             "  these are discarded from the ",
                              nboot, " bootstrap replications returned.")
    #
    cat("Computing bootstrap credible intervals.\n")
    #
    tbi = 0
    pb <- txtProgressBar(min=1,max=nrlshp*nboot,style=3)
    for(i in endogenous){
        path.rlshp[[i]] = list()
        for (p in colnames(impact[[i]])){
          path.rlshp[[i]][[p]] = list()
          path.rlshp[[i]][[p]][["replicates"]] = matrix(0,nrow = ngrid, ncol = nboot)
          for (r in 1:nboot){
            tbi = tbi + 1
            path.rlshp[[i]][[p]][["replicates"]][,r] = result[[r]]$path_coefficients$impact[[i]][,p]
            setTxtProgressBar(pb, tbi)
          }
          path.rlshp[[i]][[p]][["mean"]] = rowSums(path.rlshp[[i]][[p]][["replicates"]])/ncol(path.rlshp[[i]][[p]][["replicates"]])
          path.rlshp[[i]][[p]][["lower"]] = apply(path.rlshp[[i]][[p]][["replicates"]],1,quantile,probs=(1-conf.level)/2)
          path.rlshp[[i]][[p]][["upper"]] = apply(path.rlshp[[i]][[p]][["replicates"]],1,quantile,probs=conf.level+(1-conf.level)/2)
        }
      }
    close(pb)
    #result[[1]]$factor_scores*rep(attr(factor_scores,"scaled:scale"),each=nrow(factor_scores))+rep(attr(factor_scores,"scaled:center"),each=nrow(factor_scores))
    #result[[1]]$path_coefficients$impact$Value
    #?rowSums
    #dim(path.rlshp$Value$Quality$replicates)
    #head(path.rlshp$Value$Quality$replicates)
    #path.rlshp$Value$Quality$mean
    #path.rlshp$Value$Quality$lower
    #for (r in 1:nboot){
    #  #r = sample(1:30,1)
    #  print(r)
    #  plot(result[[r]]$path_coefficients$grid[,"Satisfaction"],result[[r]]$path_coefficients$impact[["Loyalty"]][,"Satisfaction"], type="l")
    #}
    #plot(pls.object$path_coefficients$grid[,"Satisfaction"],path.rlshp[["Loyalty"]][["Satisfaction"]][["mean"]], type="l")
    #lines(pls.object$path_coefficients$grid[,"Satisfaction"],path.rlshp[["Loyalty"]][["Satisfaction"]][["lower"]], type="l")
    #lines(pls.object$path_coefficients$grid[,"Satisfaction"],path.rlshp[["Loyalty"]][["Satisfaction"]][["upper"]], type="l")
    #
    result = list(Replicates = nboot, boot.result = result, boot.ci = path.rlshp)
    class(result) = "plsExtboot"
    return(result)
}
#########################################################################################
# replsextpm ()
# Refitting the plsextpm object for the plsextpm.boot method.
replsExtpm <- function(plsextpm, data, start=c("ones", "old")){
  start = match.arg(start)
  sum1 = plsextpm$sum1
  scaled = plsextpm$scaled
  tol = plsextpm$tolerance
  convCrit = plsextpm$convCrit
  maxit = plsextpm$maxit
  wscheme = plsextpm$model$weighting_scheme
  model = plsextpm$model
  gridLV = apply(data.frame(plsextpm$factor_scores),2,function(x){seq(min(x), max(x), length.out=plsextpm$ngrid)})
  #
  result = list()

  ## scale data?
  # Note: scale() changes class(data) to 'matrix'
  if(plsextpm$scaled) data <- scale(data)
  
  #############################################
  # step 1
  # initialization with results from fitted model
  if (!require("mgcv")) stop("package mgcv not available")
  if(start=="old"){
    Wold<-plsextpm$outer_weights

    latent = plsextpm$model$latent
    innerW.rnames = plsextpm$model$innerW.rlshp[["innerW.rnames"]]
    innerW.formula = plsextpm$model$innerW.rlshp[["formula"]]
    innerW.rlshp = plsextpm$model$innerW.rlshp
    ##
    src.latent = as.matrix(data) %*% Wold
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
    #
    factor_scores <- fscores <- Latent
    #
    i <- plsextpm$iterations
    #
    index <- plsextpm$weights_evolution$iteration==i
  }
  if(start=="ones"){
    # Weights not adding up to 1 (14.08.2009)
    # changed (30.03.2010)
    stp1 <- step1(model, data, sum1=sum1)
    Wold <- list()
    factor_scores <- stp1$Latent
    if(!sum1){
      # to ensure: w'Sw=1
      sdYs <- rep(attr(factor_scores, "scaled:scale")[model$latent],each=length(model$manifest))
      Wold <- stp1$outerW / sdYs
    } else{  
      Wold <- stp1$outerW
    }
  }
  #############################################
  # weighting scheme
  if(wscheme == "centroid") {
    innerWe <- centroid
  }
  if(wscheme == "factorial") {
    innerWe <- factorial
  }
  if(wscheme == "pathWeighting") {
    innerWe <- pathWeighting
  }
  if(wscheme == "smoothWeighting"){
    innerWe <- smoothWeighting
  }
  
  #############################################
  # Iterate over step 2 to 5
  i <- c()
  innerWeights <- c()
  eval(plsLoop.boot)

  #######################################################################
  result$factor_scores = factor_scores[,model$latent]
  attr(result$factor_scores,"scaled:center") = attr(factor_scores,"scaled:center")[model$latent]
  attr(result$factor_scores,"scaled:scale") = attr(factor_scores,"scaled:scale")[model$latent]
  #
  result$path_coefficients <- pathCoeffs(model=model,factor_scores=result$factor_scores,grid=gridLV,scaled=scaled)
  #
  if(wscheme %in% c("centroid","factorial","pathWeighting")){result$total_effects <- totEffects(model=model,pCoeff=result$path_coefficients$pC, factor_scores=factor_scores)}
  #
  result$R.Squared = rSquared(object=result)
  #
  result$cross_loadings <- cor(data, factor_scores)
  #
  result$outer_loadings <- result$cross_loadings
  #
  result$outer_loadings[Wnew==0] <- 0
  #
  result$inner_weights <- innerWeights
  #
  result$outer_weights <- Wnew
  #
  result$iterations <- (i-1)
  #
  return(result)
}

plsLoop.boot <- expression({
  #######################################################################
  # Iterate over step 2 to step 5
  i <- 1
  converged <- FALSE
  while(!converged){
    #############################################
    # step 2
    innerWeights <- innerWe(model, fscores=factor_scores)
    factor_scores <- step2(Lat=factor_scores, innerW=innerWeights, model)
    #############################################
    # step 3
    Wnew <-  step3(Latent=factor_scores, data, model,sum1=sum1)
    #############################################
    # step 4
    factor_scores <- step4(data, outerW = Wnew, model)
    if(!sum1){
      # to ensure: w'Sw=1
      sdYs <- rep(attr(factor_scores, "scaled:scale")[model$latent],each=length(model$manifest))
      Wnew <- Wnew / sdYs
    }
    #############################################
    # step 5
    stp5 <- step5(Wold, Wnew, tol, converged, convCrit)
    Wold <- stp5$Wold
    converged <- stp5$converged
    
    #############################################
    
    if(i == maxit && !converged){
      # 'try-error' especially for replsExtpm()
      class(result) <- c(class(result), "try-error")
      i <- i+1
      break
    }
    i <- i+1
  }
})

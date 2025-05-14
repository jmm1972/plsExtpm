##########################################################################################
# Estimates factor scores and parameters for PLS path models
##########################################################################################
plsExtpm <- function(data, strucmod, measuremod,...){
  UseMethod("plsextpm", model, strucmod, measuremod)
}
plsExtpm <-
function(data, strucmod, measuremod, maxit=200, tol=1e-7, scaled=TRUE, imputed = FALSE, wscheme = "centroid",
         sum1=FALSE, ngrid=300,verbose=TRUE,
         convCrit="relative"){
  result <- list(model = NULL,
                 inner_weights = NULL,
                 factor_scores = NULL,
                 ngrid = ngrid,
                 path_coefficients = NULL,
                 grid = grid,
                 R.Squared = NULL,
                 cross_loadings = NULL,
                 outer_loadings = NULL,
                 outer_weights = NULL,
                 weights_evolution = NULL,
                 scaled = scaled,
                 sum1 = sum1,
                 iterations = NULL,
                 tolerance = tol,
                 converged = NULL,
                 maxit = maxit,
                 incomplete = NULL, 
                 convCrit = NULL,
                 verbose = verbose)
  #
  class(result) <- "plsExtpm"
  if((convCrit %in% c("absolute","relative"))==FALSE){stop("Covergence criterion must be either 'relative' or 'absolute'.")}
  #
  result$convCrit = convCrit
  #
  if(!require("mgcv")) stop("package 'mgcv' not available")
  #
  # building model
  model=plsmodel(data=data, strucmod=strucmod, measuremod=measuremod, imputed, wscheme, verbose)
  #
  ## scale data?
  # Note: scale() changes class(data) to 'matrix'
  if(scaled){data <- scale(model$src.data)}
  #cat ("Init scale data: Press [enter] to continue")
  #line <- readline()
  

  #############################################
  # step 1: Initialisation
  stp1 <- step1(model,data,sum1=sum1)
  Wold <- list()
  factor_scores <- stp1$Latent
  if(!sum1){
    # to ensure: w'Sw=1
    sdYs <- rep(attr(factor_scores, "scaled:scale")[model$latent],each=length(model$manifest))
    Wold <- stp1$outerW / sdYs
  } else{  
    Wold <- stp1$outerW
  }
  weights_evolution <- reshape(as.data.frame(Wold),
                               v.names="weights",
                               ids=rownames(Wold),
                               idvar="MVs",
                               times=colnames(Wold),
                               timevar="LVs",
                               varying=list(colnames(Wold)),
                               direction="long")
  ##
  weights_evolution <- cbind(weights_evolution, iteration=0)
  #############################################
  # Select the function according to the weighting scheme
  if(model$wscheme == "centroid") {
    innerWe <- centroid
  }
  if(model$wscheme == "factorial") {
    innerWe <- factorial
  }
  if(model$wscheme == "pathWeighting") {
    innerWe <- pathWeighting
  }
  if(model$wscheme == "smoothWeighting"){
    innerWe <- smoothWeighting
  }

  converged <- c()
  i <- c()
  innerWeights <- c()
  eval(plsLoop)

  ## print
  if(converged & verbose){
      cat(paste("\nConverged after ", (i-1), " iterations.\n",
                "Tolerance: ", tol ,"\n", sep=""))
      if (model$wscheme == "centroid") cat("Scheme: centroid\n")
      if (model$wscheme == "factorial") cat("Scheme: factorial\n")
      if (model$wscheme == "pathWeighting") cat("Scheme: path weighting\n")
      if (model$wscheme == "smoothWeighting") cat("Scheme: smooth weighting\n")
  }
  if(!converged){cat(paste("\nResult did not converge after ", result$maxit, " iterations.\n",
           "\nIncrease 'maxit' and rerun.\n", sep=""))
  }
  #
  weights_evolution <- weights_evolution[weights_evolution!=0,]
  weights_evolution$LVs <- factor(weights_evolution$LVs,  levels=model$latent)
  # create result list
  #
  result$model$src.data = model$src.data
  #
  result$model$data = data
  #
  result$model$D = model$D
  #
  result$model$blocks <- model$blocks
  #
  result$model$N <- model$N
  #
  result$model$incomplete <- model$missings
  #
  result$model$strucmod = model$strucmod
  #
  result$model$measurmod = model$measuremod
  #
  result$model$latent = model$latent
  #
  result$model$endogenous = model$endogenous
  #
  result$model$exogenous = model$exogenous
  #
  result$model$manifest = model$manifest
  #
  result$model$innerW.exosucc = model$innerW.exosucc
  #
  result$model$innerW.rlshp = model$innerW.rlshp
  #
  result$model$M <- model$M
  #
  result$model$weighting_scheme = model$wscheme
  #
  result$inner_weights <- innerWeights
  #
  result$factor_scores = factor_scores[,model$latent]
  attr(result$factor_scores,"scaled:center") = attr(factor_scores,"scaled:center")[model$latent]
  attr(result$factor_scores,"scaled:scale") = attr(factor_scores,"scaled:scale")[model$latent]
  #
  result$ngrid = ngrid
  #
  gridLV = apply(data.frame(result$factor_scores),2,function(x){seq(min(x), max(x), length.out=ngrid)})
  result$path_coefficients <- pathCoeffs(model=model,factor_scores=result$factor_scores,grid=gridLV,scaled=result$scaled)
  #
  if(model$wscheme %in% c("centroid","factorial","pathWeighting")){result$total_effects <- totEffects(model=model,pCoeff=result$path_coefficients$pC, factor_scores=factor_scores)}
  #
  result$R.Squared = rSquared(object=result)
  #
  result$path_coefficients[["predicted"]] = result$path_coefficients[["predicted"]][,result$model$endogenous]
  #
  result$cross_loadings <- cor(data, result$factor_scores)
  #
  result$outer_loadings <- result$cross_loadings
  result$outer_loadings[Wnew==0] <- 0
  #
  result$outer_weights <- Wnew
  #
  result$weights_evolution <- weights_evolution
  #
  result$sum1 <- sum1
  #
  result$iterations <- (i-1)
  #
  result$tolerance <- tol
  #
  result$converged <- converged
  #
  if(converged) result$gof = gof(result)
  #
  class(result) <- "plsExtpm"
  #
  return(result)
}

plsLoop <- expression({
  #######################################################################
  # Iterate over step 2 to step 5
  i <- 1
  converged <- FALSE
  cat("Iteration:     ")
  while(!converged){
    treil <- paste(rep(" ", floor(log10(maxit)) - floor(log10(i))), collapse="")
    ndel <- paste(rep("\b", floor(log10(maxit)) + 1), collapse="")
    cat(paste(ndel, treil, i, sep=""))
    #############################################
    # step 2
    innerWeights = innerWe(model, fscores=factor_scores)
    factor_scores <- step2(Lat=factor_scores, innerW=innerWeights, model)
    #############################################
    # step 3
    Wnew <-  step3(Latent=factor_scores, data, model,sum1=sum1)
    #############################################
    # step 4
    factor_scores <- step4(data, outerW=Wnew, model)
    if(!sum1){
      # to ensure: w'Sw=1
      sdYs <- rep(attr(factor_scores, "scaled:scale")[model$latent],each=length(model$manifest))
      Wnew <- Wnew / sdYs
    }
    weights_evolution_tmp <- reshape(as.data.frame(Wnew),
                                     v.names="weights",
                                     ids=rownames(Wnew),
                                     idvar="MVs",
                                     times=colnames(Wnew),
                                     timevar="LVs",
                                     varying=list(colnames(Wnew)),
                                     direction="long")
    weights_evolution_tmp <- cbind(weights_evolution_tmp, iteration=i)
    weights_evolution <- rbind(weights_evolution, weights_evolution_tmp)
    #############################################
    # step 5
    stp5 <- step5(Wold, Wnew, tol, converged, convCrit)
    Wold <- stp5$Wold
    converged <- stp5$converged

    #############################################

    if(i == maxit && !converged){
      # 'try-error' especially for replsextpm.R
      class(result) <- c(class(result), "try-error")
      i <- i+1
      break
    }

    i <- i+1
  }
})

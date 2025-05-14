#########################################################################
# pathCoeffs()
# Calculates the matrix of path coefficients.
pathCoeffs <- function(model,factor_scores,grid,scaled){
  #
  innerW.rnames = model$innerW.rlshp[["innerW.rnames"]]
  innerW.formula = model$innerW.rlshp[["formula"]]
  innerW.rlshp = model$innerW.rlshp
  latent = model$latent
  endog = model$endogenous
  ##
  path.rlshp = list()
  path.rlshp[["pC"]] = matrix(0,ncol=length(latent),nrow=(length(innerW.rnames)))
  colnames(path.rlshp[["pC"]]) <- latent
  rownames(path.rlshp[["pC"]]) <- c(innerW.rnames)
  path.rlshp[["predicted"]] = matrix(0,ncol=length(latent),nrow=nrow(factor_scores))
  colnames(path.rlshp[["predicted"]]) <- latent
  path.rlshp[["impact"]] = list()
  #
  for(i in endog){
    indpnt <- factor_scores[,innerW.rlshp[[i]][["predictors"]]]
    dpnt <- factor_scores[,i]
    new.data <- data.frame(grid[,innerW.rlshp[[i]][["predictors"]]]); names(new.data) = innerW.rlshp[[i]][["predictors"]]
    #grid.indpnt <- grid[,innerW.rlshp[[i]][["predictors"]]]
    #grid.dpnt <- grid[,i]
    mod.data = data.frame(dpnt,indpnt); names(mod.data) = c(i,innerW.rlshp[[i]][["predictors"]])
    #new.data = data.frame(grid.dpnt,grid.indpnt); names(new.data) = c(i,innerW.rlshp[[i]][["predictors"]])
    #new.data = data.frame(grid.dpnt,grid.indpnt); names(new.data) = c(i,innerW.rlshp[[i]][["predictors"]])
    cf = coef(gam(formula(innerW.rlshp[[i]][["formula"]]),data=mod.data))
    mod.X = predict(gam(formula(innerW.rlshp[[i]][["formula"]]),data=mod.data),newdata=new.data,type="lpmatrix")
    path.rlshp[["predicted"]][,i] = predict(gam(formula(innerW.rlshp[[i]][["formula"]]),data=mod.data))
    #path.rlshp[["pC"]][names(cf),i] <- cf[-1]
    path.rlshp[["pC"]][names(cf[-1]),i] <- cf[-1]
    #path.rlshp[["pC"]]["aic",i] <- AIC(gam(formula(innerW.rlshp[[i]][["formula"]]),data=mod.data))
    path.rlshp[["impact"]][[i]] = matrix(1,nrow=nrow(grid),ncol=1)
    if(length(innerW.rlshp[[i]][["predictors"]])!=1 & length(cf)!=1){
      for (p in innerW.rlshp[[i]][["predictors"]]){
        cf.p = cf
        cf.p[!grepl(p,names(cf))] = 0
        if(!scaled){
          #path.rlshp[["impact"]][[i]] = cbind(path.rlshp[["impact"]][[i]],predict(gam(formula(innerW.rlshp[[i]][["formula"]]),data=mod.data),newdata=new.data,exclude=setdiff(innerW.rlshp[[i]][["predictors"]],p))*attr(factor_scores,"scaled:scale")[i] + attr(factor_scores,"scaled:center")[i])
          path.rlshp[["impact"]][[i]] = cbind(path.rlshp[["impact"]][[i]],mod.X%*%cf.p*attr(factor_scores,"scaled:scale")[i] + attr(factor_scores,"scaled:center")[i])
        } else{
          #path.rlshp[["impact"]][[i]] = cbind(path.rlshp[["impact"]][[i]],predict(gam(formula(innerW.rlshp[[i]][["formula"]]),data=mod.data),newdata=new.data,exclude=setdiff(innerW.rlshp[[i]][["predictors"]],p)))
          path.rlshp[["impact"]][[i]] = cbind(path.rlshp[["impact"]][[i]],mod.X%*%cf.p)
        }
      }
      path.rlshp[["impact"]][[i]] = as.matrix(path.rlshp[["impact"]][[i]][,-1])
      colnames(path.rlshp[["impact"]][[i]]) = innerW.rlshp[[i]][["predictors"]]
    }
    if(length(innerW.rlshp[[i]][["predictors"]])==1 & length(cf)!=1){
      if(!scaled){
        #path.rlshp[["impact"]][[i]] = cbind(path.rlshp[["impact"]][[i]],predict(gam(formula(innerW.rlshp[[i]][["formula"]]),data=mod.data),newdata=new.data,exclude=setdiff(innerW.rlshp[[i]][["predictors"]],p))*attr(factor_scores,"scaled:scale")[i] + attr(factor_scores,"scaled:center")[i])
        path.rlshp[["impact"]][[i]] = cbind(path.rlshp[["impact"]][[i]],mod.X%*%cf*attr(factor_scores,"scaled:scale")[i] + attr(factor_scores,"scaled:center")[i])
      } else{
        #path.rlshp[["impact"]][[i]] = cbind(path.rlshp[["impact"]][[i]],predict(gam(formula(innerW.rlshp[[i]][["formula"]]),data=mod.data),newdata=new.data,exclude=setdiff(innerW.rlshp[[i]][["predictors"]],p)))
        path.rlshp[["impact"]][[i]] = cbind(path.rlshp[["impact"]][[i]],mod.X%*%cf)
      }
      path.rlshp[["impact"]][[i]] = as.matrix(path.rlshp[["impact"]][[i]][,-1])
      colnames(path.rlshp[["impact"]][[i]]) = innerW.rlshp[[i]][["predictors"]]
    }
    if(length(innerW.rlshp[[i]][["predictors"]])==1 & length(cf)==1){
      if(!scaled){
        #path.rlshp[["impact"]][[i]] = cbind(path.rlshp[["impact"]][[i]],predict(gam(formula(innerW.rlshp[[i]][["formula"]]),data=mod.data),newdata=new.data)*attr(factor_scores,"scaled:scale")[i] + attr(factor_scores,"scaled:center")[i])
        path.rlshp[["impact"]][[i]] = cbind(path.rlshp[["impact"]][[i]],mod.X*cf*attr(factor_scores,"scaled:scale")[i] + attr(factor_scores,"scaled:center")[i])
      } else{
        #path.rlshp[["impact"]][[i]] = cbind(path.rlshp[["impact"]][[i]],predict(gam(formula(innerW.rlshp[[i]][["formula"]]),data=mod.data),newdata=new.data))
        #path.rlshp[["impact"]][[i]] = cbind(path.rlshp[["impact"]][[i]],mod.X*cf)
        path.rlshp[["impact"]][[i]] = cbind(path.rlshp[["impact"]][[i]],mod.X%*%cf)
      }
      path.rlshp[["impact"]][[i]] = as.matrix(path.rlshp[["impact"]][[i]][,-1])
      colnames(path.rlshp[["impact"]][[i]]) = innerW.rlshp[[i]][["predictors"]]
    }
  }
  if(!scaled){
    path.rlshp[["grid"]] = grid*rep(attr(factor_scores,"scaled:scale")[i],each=nrow(grid)) + rep(attr(factor_scores,"scaled:center")[i],each=nrow(grid))
  } else{
    path.rlshp[["grid"]] = grid
  }
  return(path.rlshp)
}

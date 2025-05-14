###################################################################################################
plsmodel <- function(data, strucmod, measuremod, imputed, wscheme, verbose)
{
  result <- list()
  #########################################################################################################################
  # checking structural and measurements models ###########################################################################
  #########################################################################################################################
  #
  if(!is.data.frame(strucmod)){
    strucmod <- as.data.frame(strucmod)
  }
  if(!is.matrix(measuremod)){
    measuremod <- as.matrix(measuremod)
  }
  #
  if(any(ncol(strucmod)!=4 | names(strucmod) != c("source","target","type","K") | class(strucmod)!="data.frame"))
    stop("The argument 'strucmod' must be a four-column data frame.\n",
         "Columns should be named 'source', 'target', 'type' and 'K'.\n",
         "Columns 'source' and 'target' determine the direction of the partial relationship.\n")
  #
  if(any(any(is.na(strucmod[,3])) | !match(unique(strucmod[,3]),c("ln","nln")))){
    stop("Column 'type' indicate the nature of the partial relationship, either 'linear' ('ln') or nonlinear ('nln').\n")
  }
  #
  if(!all(strucmod[,3]=="ln")){
    if(!is.numeric(strucmod[which(strucmod[,3]=="nln"),4]))
      stop("Column 'K' indicate the basis (integer) dimension of the smooth terms of the nonliner relationship.\n
           It is ignored for the linear ones.\n
           It should be an integer value.")
    strucmod$K = as.integer(strucmod$K)
  } 
  #
  if(all(ncol(measuremod)!=2 | colnames(measuremod)!=c("source","target") | mode(measuremod)!="character" | class(measuremod)!="matrix"))
    stop("The argument 'measuremod' must be a two column character matrix.\n",
         "Columns should be named 'source' and 'target'.\n",
         "Position of MV or LV as 'source' or 'target' determine mode 'A' or 'B'.")
  #
  result$strucmod <- strucmod
  #
  if(!(wscheme %in% c("centroid","factorial","pathWeighting","smoothWeighting"))){
    stop("Inner weighting scheme must be a valid option: 'centroid', or 'factorial', or 'pathWeighting', or 'smoothWeighting'.\n")
  }
  if(length(unique(strucmod[,3]))==1){
    if(unique(strucmod[,3])=="ln" & wscheme=="smoothWeighting" & verbose){
      wscheme = "centroid"
    }
  }
  #
  if("nln" %in% strucmod[,3] & wscheme!="smoothWeighting" & verbose){
      wscheme = "smoothWeighting"
      cat("All or some structural relationships are nonlinear and inner weighting scheme choosen is not 'smoothWeighting'.\n",
          "Switching to inner weighting scheme: 'smoothWeighting'.\n")
  }
  # latent variables
  latent = unique(as.vector(as.matrix(strucmod[,1:2])))
  result$latent <- latent
  #
  # manifest variables
  manifest <- sort(setdiff(as.vector(measuremod), latent))
  if(!all(manifest %in% colnames(data)))
    stop("The manifest variables must be contained in the data.")
  #
  #########################################################################################################################
  # checking data #########################################################################################################
  #########################################################################################################################
  #
  if(class(data)!="data.frame")
    stop("The argument 'data' must be of class 'data.frame'.")
  
  data <- data[, manifest]
  N <- nrow(data)
  missings <- which(complete.cases(data)==FALSE)
  if(length(missings)==0 & verbose){
    cat("All", N ,"observations are valid.\n")
  }

  if(length(missings)!=0 & imputed==FALSE & verbose){
    # Just keeping the observations, that are complete.
    data <- na.omit(data)
    cat("Data rows:", paste(missings, collapse=", "),
        "\nare not taken into acount, due to missings in the manifest variables.\n",
        "Total number of complete cases:", N-length(missings), "\n")
  }
  if(length(missings)!=0 & imputed==TRUE & verbose){
    # Imputing incomplete observations.
    data <- apply(data,2,function(x){x[is.na(x)]=mean(x,na.rm=TRUE)})
    cat("Data rows:", paste(missings, collapse=", "),
        "\nwith missing data values imputed by the mean of the respective MV.\n",
        "Total number of complete cases:", N, "\n")
    missings <- which(complete.cases(data)==FALSE)
  }
  
  ## check the variances of the data
  if(!all(apply(data, 2, sd, na.rm=TRUE)!=0)){
    stop("The MVs: ",
         paste(colnames(data)[which(apply(data, 2, sd)==0)], collapse=", "),
         "\n  have standard deviation equal to 0.\n",
         "Recheck model!\n")
  }
  #
  #########################################################################################################################
  # latent variables and structural relationships #########################################################################
  #########################################################################################################################
  #
  if(any(result$latent %in% colnames(data)))
    stop("The latent variables are not allowed to coincide with names of observed variables.\n",
         "Recheck model!\n")
  #
  # Adjacency matrix D for the structural model
  D <- innerW(strucmod, latent)
  result$endogenous = unique(strucmod[,2])
  #
  # inner relationships, predictors, types, basis dimension and number of parameters to estimate
  #
  exogenous = setdiff(strucmod[,1],strucmod[,2])
  result$exogenous = exogenous
  innerW.exosucc = list()
  for (i in exogenous){
    innerW.exosucc[[i]][["successors"]] = as.vector(strucmod[which(strucmod[,1]==i),2])
    innerW.exosucc[[i]][["successors.f"]] = character()
    innerW.exosucc[[i]][["innerW.rnames"]] = character()
    for (l in 1:nrow(strucmod)){
      if(strucmod[l,1]==i){
        if(nrow(strucmod[which(strucmod[,1] == strucmod[l,2]),c(1,3,4)])!=0){
          innerW.exosucc[[i]][["successors.f"]] = c(innerW.exosucc[[i]][["successors.f"]],unlist(apply(strucmod[which(strucmod[,1] == strucmod[l,2]),c(1,3,4)],1,function(x){ifelse(x[2]=="nln",paste0("s(",x[1],",bs='cr',k=",as.integer(x[3]),")"),x[1])}),use.names=FALSE))
          x=strucmod[which(strucmod[,1] == strucmod[l,2]),c(1,3,4)]
          for (k in 1: nrow(x)){
            if(x[k,2]=="nln"){
              innerW.exosucc[[i]][["innerW.rnames"]] = append(innerW.exosucc[[i]][["innerW.rnames"]],paste0("s(",x[k,1],").",1:(as.integer(x[k,3])-1)))
            } else if(x[k,2]=="ln"){
              innerW.exosucc[[i]][["innerW.rnames"]] = append(innerW.exosucc[[i]][["innerW.rnames"]],x[k,1])
            }
          }
        } 
        else{
            innerW.exosucc[[i]][["successors.f"]] = c(innerW.exosucc[[i]][["successors.f"]],strucmod[l,2])
            innerW.exosucc[[i]][["innerW.rnames"]] = append(innerW.exosucc[[i]][["innerW.rnames"]],strucmod[l,2])
        }
      }
    }
    names(innerW.exosucc[[i]][["successors.f"]]) = NULL
    names(innerW.exosucc[[i]][["innerW.rnames"]]) = NULL
    innerW.exosucc[[i]][["innerW.rnames"]] = unlist(innerW.exosucc[[i]][["innerW.rnames"]])
    #innerW.exosucc[[i]][["formula"]] = paste0(i, " ~ -1 + ",paste(innerW.exosucc[[i]][["successors.f"]],collapse=" + "))  
    innerW.exosucc[[i]][["formula"]] = paste0(i, " ~ ",paste(innerW.exosucc[[i]][["successors.f"]],collapse=" + "))  
  }
  ##
  innerW.rlshp = list()
  innerW.rlshp[["innerW.rnames"]] = c(latent,unlist(apply(strucmod[which(strucmod[,3]=="nln"),c(1,4)],1,function(x){paste0("s(",x[1],").",1:(as.integer(x[2])-1))}),use.names=FALSE))
  innerW.rlshp[["innerW.formula"]] = c(latent,unlist(apply(strucmod[which(strucmod[,3]=="nln"),c(1,4)],1,function(x){paste0("s(",x[1],",bs='cr',k=",as.integer(x[2]),")")}),use.names=FALSE))
  names(innerW.rlshp[["innerW.formula"]]) = NULL
  df=0
  for (i in latent){
    if(length(as.vector(strucmod[which(strucmod[,2]==i),1]))!=0){
      innerW.rlshp[[i]][["predictors"]] = as.vector(strucmod[which(strucmod[,2]==i),1])
      innerW.rlshp[[i]][["type"]] = as.vector(strucmod[which(strucmod[,2]==i),3])
      innerW.rlshp[[i]][["K"]] = as.vector(strucmod[which(strucmod[,2]==i),4])
      if(length(which(!is.na(innerW.rlshp[[i]][["K"]]) & innerW.rlshp[[i]][["K"]] < 3))!=0){
        cat("Nonlinear relationships with basis dimension, K, fewer than allowed. Increased to the minimum possible.")
        innerW.rlshp[[i]][["K"]][which(!is.na(innerW.rlshp[[i]][["K"]]) & innerW.rlshp[[i]][["K"]] < 3)] = 3
      }
      innerW.rlshp[[i]][["predictors.f"]] = unlist(apply(strucmod[which(strucmod[,2]==i),c(1,3,4)],1,function(x){ifelse(x[2]=="nln",paste0("s(",x[1],",bs='cr',k=",as.integer(x[3]),")"),x[1])}),use.names=FALSE)
      #innerW.rlshp[[i]][["formula"]] = paste0(i, " ~ -1 + ",paste(innerW.rlshp[[i]][["predictors.f"]],collapse=" + "))
      innerW.rlshp[[i]][["formula"]] = paste0(i, " ~ ",paste(innerW.rlshp[[i]][["predictors.f"]],collapse=" + "))
      innerW.rlshp[[i]][["predictors.f"]] = NULL
      innerW.rlshp[[i]][["param"]] = as.vector(strucmod[which(strucmod[,2]==i),4])
      innerW.rlshp[[i]][["param"]][is.na(innerW.rlshp[[i]][["param"]])] = 1
      innerW.rlshp[[i]][["param"]][which(innerW.rlshp[[i]][["type"]]=="nln")] = innerW.rlshp[[i]][["param"]][which(innerW.rlshp[[i]][["type"]]=="nln")]-1
      innerW.rlshp[[i]][["param.tot"]] = sum(innerW.rlshp[[i]][["param"]])
      #innerW.rlshp[[i]][["predictors.f"]] = character()
      #for (p in 1:length(innerW.rlshp[[i]][["predictors"]])){
      #  if(innerW.rlshp[[i]][["type"]][p] == "ln"){innerW.rlshp[[i]][["predictors.f"]] = append(innerW.rlshp[[i]][["predictors.f"]],innerW.rlshp[[i]][["predictors"]][p])}
      #  if(innerW.rlshp[[i]][["type"]][p] == "nln"){innerW.rlshp[[i]][["predictors.f"]] = append(innerW.rlshp[[i]][["predictors.f"]],paste0("s(",innerW.rlshp[[i]][["predictors"]][p],").",1:(as.integer(innerW.rlshp[[i]]["K"][p])-1)))}
      #}
      df=max(df,innerW.rlshp[[i]][["param.tot"]], na.rm=TRUE)
    }
  }
  #
  # Checking identifiability issues
  if((N-length(missings))<df)
    stop("Identifiability issue: the number of coefficients is greater or equal than effective sample size. Increase sample size")
  #
  # building blocks of manifest variables (including 'measurement mode')
  blocks <- block(latent, manifest, measuremod)
  #
  # Ordering of MVs and measuremod
  MVs <- NULL
  mm <- NULL
  for(i in names(blocks)){
    MVs <- append(MVs, blocks[[i]])
    if(attr(blocks[[i]], "mode") == "A"){
      mm <- rbind(mm, (cbind(i, blocks[[i]])))
    }
    if(attr(blocks[[i]], "mode") == "B"){
      mm <- rbind(mm, (cbind(blocks[[i]], i)))
    }
  }
  #
  dimnames(mm) <- dimnames(measuremod)
  measuremod <- mm
  
  #
  result$manifest <- MVs # manifest variables
  result$measuremod <- measuremod
  result$D = D
  result$innerW.exosucc = innerW.exosucc
  result$innerW.rlshp = innerW.rlshp
  result$wscheme = wscheme
  result$M = initM1(result)
  result$blocks <- blocks
  result$src.data = data[,MVs]
  result$N <- N
  result$incomplete <- missings
  return(result)
}

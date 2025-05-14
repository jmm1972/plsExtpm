plsExpplot <- function(pls.object, boot.object = NULL, pages = 0, x.LV = NULL, y.LV = NULL, xlim = NULL, ylim = NULL, legend = FALSE, verbose=TRUE,...){
  UseMethod("plsExpplot")
}
#
plsExtplot = function(pls.object, boot.object = NULL, x.LV = NULL, y.LV = NULL, xlim = NULL, ylim = NULL, legend = FALSE, verbose = TRUE){
  # check inputs
  if(!is.null(boot.object) & class(boot.object) != "plsExtboot")
    stop(boot.object, " is not of class 'plsExtboot'.")
  #    
  if(class(pls.object) != "plsExtpm")
    stop(pls.object, " is not of class 'plsExtpm'.")
  #
  latent = pls.object$model$latent
  innerW.rlshp = pls.object$model$innerW.rlshp
  scaled = pls.object$scaled
  nplots = sum(sapply(pls.object$path_coefficients$impact,dim)[2,])
  impact = pls.object$path_coefficients$impact
  grid = pls.object$path_coefficients$grid
  endog = pls.object$model$endogenous
  #
  if(!scaled){
    fscores = pls.object$factor_scores*rep(attr(pls.object$factor_scores,"scaled:scale"),each=nrow(pls.object$factor_scores)) + rep(attr(pls.object$factor_scores,"scaled:center"),each=nrow(pls.object$factor_scores))
  } else {
    fscores = pls.object$factor_scores
  }
  #
  if(length(x.LV) != 1 || length(y.LV) != 1) stop("Relationship requested must be between two Latent Variables only.")
  #
  if(!(y.LV %in% endog)) stop("Target Latent Variable is not endogenous")
  #
  if(!(x.LV %in% colnames(impact[[y.LV]]))) stop("Relationship requested does not exist in the inner model.")
  #
  if(!is.null(x.LV) & !is.null(y.LV)){
    if(innerW.rlshp[[y.LV]]$type[which(innerW.rlshp[[y.LV]][["predictors"]] == x.LV)] == "nln"){
      sub = paste("(Type: ", "nonlinear)")
    } else{
      sub = paste("(Type: ", "linear)")
    }
    xlab = x.LV
    ylab = paste0("f(",x.LV,")")
    main = paste(y.LV, " = ", ylab)
    x.fscores = fscores[,x.LV]
    y.fscores = fscores[,y.LV]
    y.line = impact[[y.LV]][,x.LV]
    x.line = grid[,x.LV]
    if(is.null(xlim)){
      if(!is.null(boot.object)){
        xlim = c(min(x.fscores,boot.object[["boot.ci"]][[y.LV]][[x.LV]][["lower"]])*0.8,max(x.fscores,boot.object[["boot.ci"]][[y.LV]][[x.LV]][["upper"]])*1.2)
      } else{
        xlim = c(min(x.fscores)*0.8,max(x.fscores)*1.2)
      }
    }
    if(is.null(ylim)){
      if(!is.null(boot.object)){
        ylim = c(min(y.fscores,boot.object[["boot.ci"]][[y.LV]][[x.LV]][["lower"]])*0.8,max(y.fscores,boot.object[["boot.ci"]][[y.LV]][[x.LV]][["upper"]])*1.2)
      } else{
        ylim = c(min(y.fscores)*0.8,max(y.fscores)*1.2)
      }
    }
    plot(x.fscores,y.fscores,
         pch=19,
         cex=0.7,
         xlab=xlab,
         ylab=ylab,
         main=main,
         sub=sub,
         cex.main=1,
         cex.lab=1,
         cex.axis=1,
         col="blue",
         xlim = xlim,
         ylim = ylim,
         cex.main=1,
         bty="n",
         font=2,
         cex.sub = 0.8)
    # Linear PLS estimated relationship
    lines(x.line,y.line, lty="solid", col="blue", lwd=2)
    lines(x.line,plsderiv(x.line,y.line), lwd=2, lty=2, col="blue")
    # Linear PLS estimated credible intervals
    #plot(c(0,10),c(0,10), type="n", xaxt="n",yaxt="n",bty="n")
    #lines(x.line,boot.object[["boot.ci"]][[y.LV]][[x.LV]][["upper"]])
    #lines(x.line,boot.object[["boot.ci"]][[y.LV]][[x.LV]][["lower"]])
    #lines()
    polygon(c(x.line, rev(x.line)),
            c(boot.object[["boot.ci"]][[y.LV]][[x.LV]][["upper"]], rev(boot.object[["boot.ci"]][[y.LV]][[x.LV]][["lower"]])),
            border = NA, col = adjustcolor("blue", alpha.f = 0.3))
    if(legend){
      legend("topleft",
             c("Estimated factor scores",
               "Estimated relationship",
               "Relationship first derivative",
               "95% credibility band"
             ),
             pch=c(19,NA,
                   NA,NA
             ),
             lty=c(NA,1,
                   2,1
             ),
             col=c("blue","blue",
                   "blue",adjustcolor("blue", alpha.f = 0.3)
             ),
             lwd=c(NA,2,
                   2,10
             ),bty="n",cex=0.7, text.font=2)
    }
  } else{
  #
    for (i in endog){
      for (p in colnames(impact[[i]])){
        if(innerW.rlshp[[i]]$type[which(innerW.rlshp[[i]][["predictors"]] == p)] == "nln"){
          sub = paste("Type: ", "nonlinear")
        } else{
          sub = paste("Type: ", "linear")
        }
        xlab = p
        ylab = paste0("f(",p,")")
        main = paste(i, " = ", ylab)
        x.fscores = fscores[,p]
        y.fscores = fscores[,i]
        y.line = impact[[i]][,p]
        x.line = grid[,i]
        if(is.null(xlim)){
          if(!is.null(boot.object)){
            xlim = c(min(x.fscores,boot.object[["boot.ci"]][[y.LV]][[x.LV]][["lower"]])*0.8,max(x.fscores,boot.object[["boot.ci"]][[y.LV]][[x.LV]][["upper"]])*1.2)
          } else{
            xlim = c(min(x.fscores)*0.8,max(x.fscores)*1.2)
          }
        }
        if(is.null(ylim)){
          if(!is.null(boot.object)){
            ylim = c(min(y.fscores,boot.object[["boot.ci"]][[y.LV]][[x.LV]][["lower"]])*0.8,max(y.fscores,boot.object[["boot.ci"]][[y.LV]][[x.LV]][["upper"]])*1.2)
          } else{
            ylim = c(min(y.fscores)*0.8,max(y.fscores)*1.2)
          }
        }
        plot(x.fscores,y.fscores,
             pch=19,
             cex=0.7,
             xlab=xlab,
             ylab=ylab,
             main=main,
             sub=sub,
             cex.main=1,
             cex.lab=1,
             cex.axis=1,
             col="blue",
             xlim = xlim,
             ylim = ylim,
             cex.main=1,
             bty="n",
             font=2,
             cex.sub = 0.8)
        # Linear PLS estimated relationship
        lines(x.line,y.line, lty="solid", col="blue", lwd=2)
        lines(x.line,plsderiv(x.line,y.line), lwd=2, lty=2, col="blue")
        # Linear PLS estimated credible intervals
        polygon(c(x.line, rev(x.line)),
                c(boot.object[["boot.ci"]][[i]][[p]][["upper"]], rev(boot.object[["boot.ci"]][[i]][[p]][["lower"]])),
                border = NA, col = adjustcolor("blue", alpha.f = 0.3))
        if(legend){
          legend("topleft",
                 c("Estimated factor scores",
                   "Estimated relationship",
                  "Relationship first derivative",
                  "95% credibility band"
                ),
                pch=c(19,NA,
                      NA,NA
                ),
                lty=c(NA,1,
                      2,1
                ),
                col=c("blue","blue",
                      "blue",adjustcolor("blue", alpha.f = 0.3)
                ),
                lwd=c(NA,2,
                      2,10
                ),bty="n",cex=0.7, text.font=2)
        }
        cat ("Next plot: Press [enter] to continue")
        line <- readline()
      }
    }
  }    
}


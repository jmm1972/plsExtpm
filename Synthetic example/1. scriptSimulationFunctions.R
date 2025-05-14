###########################################################################################################################
### This script contains the functions for running the synthetic Monte Carlo process
### It also contains the functions that produce synthetic example figures and tables in the article
###
###########################################################################################################################
### run.results(): run results for PLS and PLSs
###########################################################################################################################
run.results = function(h,pop,sample.idx,samples=1000){
  PLS.results = list()
  PLSs.results = list()
  for (ss in c("n75","n100","n150","n250","n300","n500","n600","n750","n900"))
  {
    PLS.results[[ss]] = list()
    PLSs.results[[ss]] = list()
    for (s in 1:samples)
    {
      if(h=="h1"){data = sim.pop.h1[sample.index[[h]][[ss]][[s]], c(paste0("xi.",1:5),paste0("eta1.",1:5),paste0("eta2.",1:5),paste0("eta3.",1:5),paste0("eta4.",1:5),paste0("eta5.",1:5))]}
      if(h=="h2"){data = sim.pop.h2[sample.index[[h]][[ss]][[s]], c(paste0("xi.",1:5),paste0("eta1.",1:5),paste0("eta2.",1:5),paste0("eta3.",1:5),paste0("eta4.",1:5),paste0("eta5.",1:5))]}
      if(h=="h3"){data = sim.pop.h3[sample.index[[h]][[ss]][[s]], c(paste0("xi.",1:5),paste0("eta1.",1:5),paste0("eta2.",1:5),paste0("eta3.",1:5),paste0("eta4.",1:5),paste0("eta5.",1:5))]}
      cat(paste(h," - ", ss," - Sample: ", s,"\n",sep=""))
      PLS.results[[ss]][[s]] <- my.sempls(data = data, strucmod = sm.m, measuremod = mm.m, mode="PLS", maxit=200, scaled=FALSE, sum1=TRUE)
      PLS.results[[ss]][[s]]$weights_evolution = NULL
      PLS.results[[ss]][[s]]$cross_loadings = NULL
      PLS.results[[ss]][[s]]$inner_weights = NULL
      PLS.results[[ss]][[s]]$blocks = NULL
      PLS.results[[ss]][[s]]$data = NULL
      PLS.results[[ss]][[s]]$weighting_scheme = NULL
      PLS.results[[ss]][[s]]$sum1 = NULL
      PLS.results[[ss]][[s]]$tolerance = NULL
      PLS.results[[ss]][[s]]$maxit = NULL
      PLS.results[[ss]][[s]]$N = NULL
      PLS.results[[ss]][[s]]$verbose = NULL
      PLS.results[[ss]][[s]]$total_effects = NULL
      PLS.results[[ss]][[s]]$tol = NULL
      PLS.results[[ss]][[s]]$model$src.data = NULL
      PLS.results[[ss]][[s]]$model$blocks = NULL
      PLS.results[[ss]][[s]]$model$D = NULL
      PLS.results[[ss]][[s]]$model$C = NULL
      PLS.results[[ss]][[s]]$model$M = NULL
      PLS.results[[ss]][[s]]$model$N = NULL
      PLS.results[[ss]][[s]]$model$incomplete = NULL
      PLS.results[[ss]][[s]]$model$manifest = NULL
      PLS.results[[ss]][[s]]$model$strucmod = NULL
      PLS.results[[ss]][[s]]$model$measuremod = NULL
      ##
      PLSs.results[[ss]][[s]] <- my.sempls(data = data, strucmod = sm.m, measuremod = mm.m, mode="PLSs",maxit=200, scaled=FALSE, sum1=TRUE)
      PLSs.results[[ss]][[s]]$weights_evolution = NULL
      PLSs.results[[ss]][[s]]$cross_loadings = NULL
      PLSs.results[[ss]][[s]]$inner_weights = NULL
      PLSs.results[[ss]][[s]]$blocks = NULL
      PLSs.results[[ss]][[s]]$data = NULL
      PLSs.results[[ss]][[s]]$weighting_scheme = NULL
      PLSs.results[[ss]][[s]]$sum1 = NULL
      PLSs.results[[ss]][[s]]$tolerance = NULL
      PLSs.results[[ss]][[s]]$maxit = NULL
      PLSs.results[[ss]][[s]]$N = NULL
      PLSs.results[[ss]][[s]]$verbose = NULL
      PLSs.results[[ss]][[s]]$total_effects = NULL
      PLSs.results[[ss]][[s]]$tol = NULL
      PLSs.results[[ss]][[s]]$model$src.data = NULL
      PLSs.results[[ss]][[s]]$model$blocks = NULL
      PLSs.results[[ss]][[s]]$model$D = NULL
      PLSs.results[[ss]][[s]]$model$C = NULL
      PLSs.results[[ss]][[s]]$model$M = NULL
      PLSs.results[[ss]][[s]]$model$N = NULL
      PLSs.results[[ss]][[s]]$model$incomplete = NULL
      PLSs.results[[ss]][[s]]$model$manifest = NULL
      PLSs.results[[ss]][[s]]$model$strucmod = NULL
      PLSs.results[[ss]][[s]]$model$measuremod = NULL
    }
  }
  return(list(PLS=PLS.results,PLSs=PLSs.results))
}
###########################################################################################################################
### END OF run.results(): run results for PLS and PLSs
###########################################################################################################################
##
##
###########################################################################################################################
### run.converg.anal(): run convergence analysis for PLS and PLSs
###########################################################################################################################
run.converg.anal = function(h,data,samples=1000){
  cM = list()
  iter = list()
  weig.neg = list()
  PLS = list()
  PLSs = list()
  for (ss in c("n75","n100","n150","n250","n300","n500","n600","n750","n900"))
  {
    cM[[ss]] = matrix(0, 2, 2);rownames(cM[[ss]]) = c("PLS_converged","PLS_non_converged");colnames(cM[[ss]])=c("PLSs_converged","PLSs_non_converged")
    iter[[ss]] = matrix(0,2,200);colnames(iter[[ss]])=c(paste0("iter",1:200));rownames(iter[[ss]]) = c("PLS","PLSs")
    weig.neg[[ss]] = c(0,0);names(weig.neg[[ss]])=c("PLS_neg_weights","PLSs_neg_weights")
    PLS[[ss]] = numeric()
    PLSs[[ss]] = numeric()
    h.PLS=0
    h.PLSs=0
    for (s in 1:samples)
    {
      if(data[["PLS"]][[ss]][[s]]$converged==TRUE  & data[["PLSs"]][[ss]][[s]]$converged==TRUE){cM[[ss]][1,1] = cM[[ss]][1,1] + 1}
      if(data[["PLS"]][[ss]][[s]]$converged==FALSE & data[["PLSs"]][[ss]][[s]]$converged==FALSE){cM[[ss]][2,2] = cM[[ss]][2,2] + 1}
      if(data[["PLS"]][[ss]][[s]]$converged==TRUE  & data[["PLSs"]][[ss]][[s]]$converged==FALSE){cM[[ss]][1,2] = cM[[ss]][1,2] + 1}
      if(data[["PLS"]][[ss]][[s]]$converged==FALSE & data[["PLSs"]][[ss]][[s]]$converged==TRUE){cM[[ss]][2,1] = cM[[ss]][2,1] + 1}
      if(data[["PLS"]][[ss]][[s]]$converged==TRUE){iter[[ss]][1,data[["PLS"]][[ss]][[s]]$iterations] = iter[[ss]][1,data[["PLS"]][[ss]][[s]]$iterations] + 1}
      if(data[["PLSs"]][[ss]][[s]]$converged==TRUE){iter[[ss]][2,data[["PLSs"]][[ss]][[s]]$iterations] = iter[[ss]][1,data[["PLSs"]][[ss]][[s]]$iterations] + 1}
      if(data[["PLS"]][[ss]][[s]]$converged==TRUE  & sum(rowSums(data[["PLS"]][[ss]][[s]]$outer_weights)<0)>0){weig.neg[[ss]][1] = weig.neg[[ss]][1] + 1}
      if(data[["PLSs"]][[ss]][[s]]$converged==TRUE  & sum(rowSums(data[["PLSs"]][[ss]][[s]]$outer_weights)<0)>0){weig.neg[[ss]][2] = weig.neg[[ss]][2] + 1}
      if(data[["PLS"]][[ss]][[s]]$converged==TRUE){PLS[[ss]][h.PLS]=s;h.PLS=h.PLS+1}
      if(data[["PLSs"]][[ss]][[s]]$converged==TRUE){PLSs[[ss]][h.PLSs]=s;h.PLSs=h.PLSs+1}
    }
  }
  return(list(cM=cM,iter.freq=iter,weig.neg.freq=weig.neg,index=list(PLS=PLS,PLSs=PLSs)))
}
###########################################################################################################################
### END OF run.converg.anal(): run convergence analysis for PLS and PLSs
###########################################################################################################################
##
##
###########################################################################################################################
### run.converg.summary(): run convergence analysis for PLS and PLSs
###########################################################################################################################
run.converg.summary = function(data){
  converg.sum=matrix(0,9,3)
  rownames(converg.sum)=c("n75","n100","n150","n250","n300","n500","n600","n750","n900")
  colnames(converg.sum)=c("conv.PLS_conv.PLSs","ncon.PLS_conv.PLSs","conv.PLS_ncon.PLSs")
  for (ss in c("n75","n100","n150","n250","n300","n500","n600","n750","n900")){
    converg.sum[ss,1] = sum(diag(data[["cM"]][[ss]]))/sum(data[["cM"]][[ss]])
    converg.sum[ss,2] = data[["cM"]][[ss]][2,1]/sum(data[["cM"]][[ss]])
    converg.sum[ss,3] = data[["cM"]][[ss]][1,2]/sum(data[["cM"]][[ss]])
  }
  rownames(converg.sum)=c(75,100,150,250,300,500,600,750,900)
  return(converg.sum)
}
###########################################################################################################################
### END OF run.converg.summary(): run convergence analysis for PLS and PLSs
###########################################################################################################################
##
##
###########################################################################################################################
### run.wght.neg(): computes the samples with negative outer weights for PLS and PLSs
###########################################################################################################################
run.wght.neg = function(data){
  weig.neg=matrix(0,9,2)
  rownames(weig.neg)=c("n75","n100","n150","n250","n300","n500","n600","n750","n900")
  colnames(weig.neg)=c("PLS", "PLSs")
  for (ss in c("n75","n100","n150","n250","n300","n500","n600","n750","n900")){
    weig.neg[ss,]=data[["weig.neg.freq"]][[ss]]/c(data[["cM"]][[ss]][1,1]+data[["cM"]][[ss]][1,2],data[["cM"]][[ss]][1,1]+data[["cM"]][[ss]][2,1])
  }
  rownames(weig.neg)=c(75,100,150,250,300,500,600,750,900)
  return(weig.neg)
}
###########################################################################################################################
### END OF run.wght.neg(): computes the samples with negative outer weights for PLS and PLSs
###########################################################################################################################
##
##
###########################################################################################################################
### run.mean.iter.neg(): computes mean number of iterations to obtain converge convergence analysis for PLS and PLSs
###########################################################################################################################
run.mean.iter = function(data){
  iter.med=matrix(0,9,2)
  rownames(iter.med)=c("n75","n100","n150","n250","n300","n500","n600","n750","n900")
  colnames(iter.med)=c("PLS", "PLSs")
  for (ss in c("n75","n100","n150","n250","n300","n500","n600","n750","n900")){
    iter.med[ss,]=c(sum(data[["iter.freq"]][[ss]][1,]*1:200)/(data[["cM"]][[ss]][1,1]+data[["cM"]][[ss]][1,2]),
                    sum(data[["iter.freq"]][[ss]][2,]*1:200)/(data[["cM"]][[ss]][1,1]+data[["cM"]][[ss]][2,1]))
  }
  rownames(iter.med)=c(75,100,150,250,300,500,600,750,900)
  return(iter.med)
}
###########################################################################################################################
### END OF run.mean.iter.neg(): computes mean number of iterations to obtain converge convergence analysis for PLS and PLSs
###########################################################################################################################
##
##
###########################################################################################################################
### run.sample.relshp(): computes the relationships between latent variables for PLS and PLSs
###########################################################################################################################
run.sample.relshp<-function(h,data,grid,index,center,scale,var=c("eta1","eta2","eta3","eta4","eta5"), samples=1000)
{
  #
  relshp = list()

  for (m in c("PLS","PLSs"))
  {
    relshp[[m]] = list()
    for (ss in c("n75","n100","n150","n250","n300","n500","n600","n750","n900"))
    {
      relshp[[m]][[ss]] = list()
      for (v in var)
      {
        relshp[[m]][[ss]][[v]] = list()
        relshp[[m]][[ss]][[v]][["rep"]] =  matrix(0,300,length(index[[m]][[ss]]))
        for (s in 1:samples)
        {
          cat(paste("Communality: ", h, " - Method: ", m, " - Sample size: ", ss," - Variable: ", v, " - Sample: ", s,"\n",sep=""))
          PLS.factor_scores = data.frame(data[[m]][[ss]][[s]]$factor_scores)
          PLSs.factor_scores = data.frame(data[[m]][[ss]][[s]]$factor_scores)
          form.PLS=paste(v,"~ xi")
          form.PLSs=paste(v, "~ s(xi, bs='cr',k=10)")
          if (m=="PLS"){relshp[[m]][[ss]][[v]][["rep"]][,which(index[[m]][[ss]]==s)]=predict(lm(formula(form.PLS),data=PLS.factor_scores),newdata=data.frame(xi=grid))*rep(scale[v],length(grid))+rep(center[v],length(grid))}
          if (m=="PLSs"){relshp[[m]][[ss]][[v]][["rep"]][,which(index[[m]][[ss]]==s)]=predict(gam(formula(form.PLSs),data=PLSs.factor_scores),newdata=data.frame(xi=grid))*rep(scale[v],length(grid))+rep(center[v],length(grid))}
        }
        relshp[[m]][[ss]][[v]][["mean"]]=rowSums(relshp[[m]][[ss]][[v]][["rep"]]/length(index[[m]][[ss]]))
        #
        relshp[[m]][[ss]][[v]][["lower"]]=apply(relshp[[m]][[ss]][[v]][["rep"]],1,quantile,probs=0.05/2)
        #
        relshp[[m]][[ss]][[v]][["upper"]]=apply(relshp[[m]][[ss]][[v]][["rep"]],1,quantile,probs=1-0.05/2)
      }
    }
  }
  return(relshp)
}
###########################################################################################################################
### END OF run.sample.relshp(): computes the relationships between latent variables for PLS and PLSs
###########################################################################################################################
##
##
###########################################################################################################################
### plot.impact(): plots relationships between latent variables for PLS and PLSs and respective 95% credible intervals
###########################################################################################################################
plot.impact.fn = function(h,n,pop,main,var,data, grid){
  ymin=min(pop[,var],data[[h]][["PLS"]][[n]][[var]][["mean"]],data[[h]][["PLS"]][[n]][[var]][["lower"]],data[[h]][["PLS"]][[n]][[var]][["upper"]],
           data[[h]][["PLSs"]][[n]][[var]][["mean"]],data[[h]][["PLSs"]][[n]][[var]][["lower"]],data[[h]][["PLSs"]][[n]][[var]][["upper"]])
  ymax=max(pop[,var],data[[h]][["PLS"]][[n]][[var]][["mean"]],data[[h]][["PLS"]][[n]][[var]][["lower"]],data[[h]][["PLS"]][[n]][[var]][["upper"]],
           data[[h]][["PLSs"]][[n]][[var]][["mean"]],data[[h]][["PLSs"]][[n]][[var]][["lower"]],data[[h]][["PLSs"]][[n]][[var]][["upper"]])
  if(var=="eta1"){var.plot=expression(eta[1])}
  if(var=="eta2"){var.plot=expression(eta[2])}
  if(var=="eta3"){var.plot=expression(eta[3])}
  if(var=="eta4"){var.plot=expression(eta[4])}
  if(var=="eta5"){var.plot=expression(eta[5])}
  plot(pop[,"xi"],pop[,var],
       ylim=c(ymin,ymax),
       main=NULL,
       pch=19,type="p",
       xlab=expression(xi),
       ylab=expression("f("*xi*")"), col="lightgray",cex=0.5,font=2,cex.main=1,bty="n")
  lines(pop[,"xi"],pop[,paste0(var,"f")],lwd=2)
  lines(grid,data[[h]][["PLSs"]][[n]][[var]][["mean"]],type="l",lwd=2, col="red")
  # Plot of 2.5% and 97.5% quantiles of 1000 PLSs
  polygon(c(grid, rev(grid)),
          c(data[[h]][["PLSs"]][[n]][[var]][["upper"]], rev(data[[h]][["PLSs"]][[n]][[var]][["lower"]])),
          border = NA, col = adjustcolor("red", alpha.f = 0.2))
  # Linear PLS line
  lines(grid,data[[h]][["PLS"]][[n]][[var]][["mean"]],type="l",lwd=2,col="blue")
  # Plot of 2.5% and 97.5% quantiles of 1000 PLS
  polygon(c(grid, rev(grid)),
          c(data[[h]][["PLS"]][[n]][[var]][["upper"]], rev(data[[h]][["PLS"]][[n]][[var]][["lower"]])),
          border = NA, col = adjustcolor("blue", alpha.f = 0.2))
}
###########################################################################################################################
### END OF plot.impact(): plots relationships between latent variables for PLS and PLSs and respective 95% credible intervals
###########################################################################################################################
##
##
###########################################################################################################################
### relshp.accuracy: computes RMSE and Bias for PLS and PLSs ##############################################################
###########################################################################################################################
relshp.accuracy<-function(data, mu, index,grid,var=c("eta1","eta2","eta3","eta4","eta5")){
  #
  RMSE=list()
  B=list()
  #
  for (h in c("h1","h2","h3"))
  {
    RMSE[[h]] = list()
    B[[h]] = list()
    for (m in c("PLS","PLSs"))
    {
      RMSE[[h]][[m]] = list()
      B[[h]][[m]] = list()
      for (ss in c("n75","n100","n150","n250","n300","n500","n600","n750","n900"))
        {
          RMSE[[h]][[m]][[ss]] = list()
          B[[h]][[m]][[ss]] = list()
          for (v in var)
          {
            cat(paste("Communality: ", h, " - Method: ", m, " - Sample size: ", ss," - Variable: ", v,"\n",sep=""))
            RMSE[[h]][[m]][[ss]][[v]] = sqrt(sum((data[[h]][[m]][[ss]][[v]][["rep"]]-mu[,v])^2)/(length(grid)*length(index[[h]][["index"]][[m]][[ss]])))
            B[[h]][[m]][[ss]][[v]] = sum(abs(data[[h]][[m]][[ss]][[v]][["mean"]]-mu[,v]))/length(grid)
          }
        }
      }
    }
  return(list(RMSE=RMSE,B=B))
}
###########################################################################################################################
### END OF relshp.accuracy: computes RMSE and Bias for PLS and PLSs #######################################################
###########################################################################################################################
##
##
###########################################################################################################################
### plot.converg: plots a matrix of cofsusion matrices ####################################################################
###########################################################################################################################
plot.confusion.matrix = function(data)
  {
  layout(matrix(c(1,2,3,4,5,6,
                  7,8,9,10,11,12,
                  13,14,15,16,17,18,
                  19,20,21,22,23,24,
                  19,25,26,27,28,29,
                  19,30,31,32,33,34,
                  19,35,36,37,38,39,
                  19,40,41,42,43,44,
                  19,45,46,47,48,49,
                  19,50,51,52,53,54,
                  19,55,56,57,58,59),
                nrow=11,byrow=TRUE),
         heights=c(1.5,1.5,1.5,3,3,3,3,3,3,3,3,3),
         widths=c(0.5,0.5,1.5,3,3,3))
  #layout.show(n=59)
  par(mar=c(0,0,1,0))
  plot.new() #1
  plot.new() #2
  plot.new() #3
  plot.new() #4
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  text((90+300)/2, (350+450)/2.1, 'PLSs-PM', cex=1.2, adj=0.5, font=2)
  plot.new() #6
  plot.new() #7
  plot.new() #8
  plot.new() #9
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  text((90+300)/2+30, (350+450)/2.1, 'Communality = 25%', cex=1, adj=0.5, font=2)
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  text((90+300)/2+30, (350+450)/2.1, 'Communality = 50%', cex=1, adj=0.5, font=2)
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  text((90+300)/2+30, (350+450)/2.1, 'Communality = 75%', cex=1, adj=0.5, font=2)
  plot.new() #11
  plot.new() #13
  plot.new() #14
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  text((100+215)/2, (380+440)/2.1, 'Converged', cex=0.7, adj=0.5, font=2)
  text((225+340)/2, (380+440)/2.1, 'Non-converged', cex=0.7, adj=0.5, font=2)
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  text((100+215)/2, (380+440)/2.1, 'Converged', cex=0.7, adj=0.5, font=2)
  text((225+340)/2, (380+440)/2.1, 'Non-converged', cex=0.7, adj=0.5, font=2)
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  text((100+215)/2, (380+440)/2.1, 'Converged', cex=0.7, adj=0.5, font=2)
  text((225+340)/2, (380+440)/2.1, 'Non-converged', cex=0.7, adj=0.5, font=2)
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  text((90+300)/2, (350+450)/2.1, 'PLS-PM', cex=1.2, adj=0.5,srt=90, font=2)
  ##
  plot.matrix.cM.l1(CM=data[["line1"]])
  plot.matrix.cM(data[["line2"]])
  plot.matrix.cM(data[["line3"]])
  plot.matrix.cM(data[["line4"]])
  plot.matrix.cM(data[["line5"]])
  plot.matrix.cM(data[["line6"]])
  plot.matrix.cM(data[["line7"]])
  plot.matrix.cM(data[["line8"]])
  par(mfrow=c(1,1))
}
##
plot.matrix.cM.l1 = function(CM){
  #
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  text(195, 375, paste0("n = ",CM$n), cex=0.7, srt=90, adj=0.5, font=2)
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  text(157.5, 400, 'Converged', cex=0.7, adj=0, font=2)
  text(157.5, 335, 'Non-converged', cex=0.7,adj=0, font=2)
  #
  plot(c(90, 350), c(300, 450), type = "n", main = "",
     xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  rect(100, 380, 215, 440, col='#3F97D0')
  rect(225, 380, 340, 440, col='#F7AD50')
  rect(100, 310, 215, 370, col='#F7AD50')
  rect(225, 310, 340, 370, col='#3F97D0')

  # add in the cm results 
  res <- as.numeric(CM[[1]])
  text((100+215)/2, (380+440)/2, paste0(res[1],"%"), cex=1, font=2, col='white')
  text((225+340)/2, (380+440)/2, paste0(res[2],"%"), cex=1, font=2, col='white')
  text((100+215)/2, (310+370)/2, paste0(res[3],"%"), cex=1, font=2, col='white')
  text((225+340)/2, (310+370)/2, paste0(res[4],"%"), cex=1, font=2, col='white')
  ##
  ##
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  rect(100, 380, 215, 440, col='#3F97D0')
  rect(225, 380, 340, 440, col='#F7AD50')
  rect(100, 310, 215, 370, col='#F7AD50')
  rect(225, 310, 340, 370, col='#3F97D0')
  
  # add in the cm results 
  res <- as.numeric(CM[[2]])
  text((100+215)/2, (380+440)/2, paste0(res[1],"%"), cex=1, font=2, col='white')
  text((225+340)/2, (380+440)/2, paste0(res[2],"%"), cex=1, font=2, col='white')
  text((100+215)/2, (310+370)/2, paste0(res[3],"%"), cex=1, font=2, col='white')
  text((225+340)/2, (310+370)/2, paste0(res[4],"%"), cex=1, font=2, col='white')
  ##
  ##
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  rect(100, 380, 215, 440, col='#3F97D0')
  rect(225, 380, 340, 440, col='#F7AD50')
  rect(100, 310, 215, 370, col='#F7AD50')
  rect(225, 310, 340, 370, col='#3F97D0')
  
  # add in the cm results 
  res <- as.numeric(CM[[3]])
  text((100+215)/2, (380+440)/2, paste0(res[1],"%"), cex=1, font=2, col='white')
  text((225+340)/2, (380+440)/2, paste0(res[2],"%"), cex=1, font=2, col='white')
  text((100+215)/2, (310+370)/2, paste0(res[3],"%"), cex=1, font=2, col='white')
  text((225+340)/2, (310+370)/2, paste0(res[4],"%"), cex=1, font=2, col='white')
}
##################### LINE 2 ##############################
plot.matrix.cM = function(CM) {
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  text(195, 375, paste0("n = ",CM$n), cex=0.7, srt=90, adj=0.5, font=2)
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  text(157.5, 400, 'Converged', cex=0.7, adj=0, font=2)
  text(150, 335, 'Non-converged', cex=0.7,adj=0, font=2)
  #
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  rect(100, 380, 215, 440, col='#3F97D0')
  rect(225, 380, 340, 440, col='#F7AD50')
  rect(100, 310, 215, 370, col='#F7AD50')
  rect(225, 310, 340, 370, col='#3F97D0')
  
  # add in the cm results 
  res <- as.numeric(CM[[1]])
  text((100+215)/2, (380+440)/2, paste0(res[1],"%"), cex=1, font=2, col='white')
  text((225+340)/2, (380+440)/2, paste0(res[2],"%"), cex=1, font=2, col='white')
  text((100+215)/2, (310+370)/2, paste0(res[3],"%"), cex=1, font=2, col='white')
  text((225+340)/2, (310+370)/2, paste0(res[4],"%"), cex=1, font=2, col='white')
  ##
  ##
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  rect(100, 380, 215, 440, col='#3F97D0')
  rect(225, 380, 340, 440, col='#F7AD50')
  rect(100, 310, 215, 370, col='#F7AD50')
  rect(225, 310, 340, 370, col='#3F97D0')
  
  # add in the cm results 
  res <- as.numeric(CM[[2]])
  text((100+215)/2, (380+440)/2, paste0(res[1],"%"), cex=1, font=2, col='white')
  text((225+340)/2, (380+440)/2, paste0(res[2],"%"), cex=1, font=2, col='white')
  text((100+215)/2, (310+370)/2, paste0(res[3],"%"), cex=1, font=2, col='white')
  text((225+340)/2, (310+370)/2, paste0(res[4],"%"), cex=1, font=2, col='white')
  ##
  ##
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  rect(100, 380, 215, 440, col='#3F97D0')
  rect(225, 380, 340, 440, col='#F7AD50')
  rect(100, 310, 215, 370, col='#F7AD50')
  rect(225, 310, 340, 370, col='#3F97D0')
  
  # add in the cm results 
  res <- as.numeric(CM[[3]])
  text((100+215)/2, (380+440)/2, paste0(res[1],"%"), cex=1, font=2, col='white')
  text((225+340)/2, (380+440)/2, paste0(res[2],"%"), cex=1, font=2, col='white')
  text((100+215)/2, (310+370)/2, paste0(res[3],"%"), cex=1, font=2, col='white')
  text((225+340)/2, (310+370)/2, paste0(res[4],"%"), cex=1, font=2, col='white')
}
###########################################################################################################################
### plot.converg: plots a matrix of consusion matrices ####################################################################
###########################################################################################################################
##
##
###########################################################################################################################
### plot.accuracy: plots RMSE and absolute bias for PLS and PLSs ##########################################################
###########################################################################################################################
plot.accuracy = function(data,measure,var=c("eta1","eta2","eta3","eta4","eta5"))
{
  layout(matrix(c(1,2,2,2,2,2,3,
                  4,5,6,7,8,9,3,
                  10,11,12,13,14,15,3,
                  10,16,16,16,16,16,3,
                  10,17,18,19,20,21,3,
                  10,22,22,22,22,22,3,
                  10,23,24,25,26,27,3,
                  10,28,28,28,28,28,3),
                nrow=8,byrow=TRUE),
         heights=c(0.8,0.5,3,0.8,3,0.8,3,1),
         widths=c(2,3,3,3,3,3,1))
  #layout.show(n=28)
  par(mar=c(0,0,0,0))
  plot.new() #1
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  text((90+300)/2+25, (350+450)/2-25, 'Communality = 25%', cex=1.2, adj=0.5, font=2)
  legend("right", c("PLS-PM","PLSs-PM"), col=c("blue","red"), bty="n", pch=19,text.font=2,cex=1.2)
  plot.new() #3
  plot.new() #4
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  text((90+300)/2+5, (350+450)/2.1, expression(eta[1]), cex=1.2, adj=0.5, font=2)
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  text((90+300)/2+5, (350+450)/2.1, expression(eta[2]), cex=1.2, adj=0.5, font=2)
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  text((90+300)/2+5, (350+450)/2.1, expression(eta[3]), cex=1.2, adj=0.5, font=2)
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  text((90+300)/2, (350+450)/2.1, expression(eta[4]), cex=1.2, adj=0.5, font=2)
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  text((90+300)/2+5, (350+450)/2.1, expression(eta[5]), cex=1.2, adj=0.5, font=2)
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  if(measure=="RMSE"){text((90+300)/2, (350+450)/2.1, "Root Mean Square Error", cex=1.2, adj=0.5, font=2,srt=90)}
  if(measure=="Bias"){text((90+300)/2, (350+450)/2.1, "Absolute bias", cex=1.2, adj=0.5, font=2,srt=90)}
  #
  x=c("75","100","150","250","300","500","750","900")
  y=data[which(data$LV==var[1] & data$h=="h1"),]
  plot(y[which(y$Method=="PLS"),measure], type = "p", main = "",pch=19, ylim=c(0,2.5),
       xlab="", ylab="", xaxt='n', yaxt="n", col="blue",cex=1.2)
  points(y[which(y$Method=="PLSs"),measure], type = "p", pch=19,col="red",cex=1.2)
  axis(2,at=c(0,0.5,1,1.5,2,2.5),cex.axis=0.8,font.axis=2, font.lab=2)
  abline(v=c(1,2,3,4,5,6,7,8), col="gray", lty=3)
  for(v in var[-1])
  {
    y=data[which(data$LV==v & data$h=="h1"),]
    plot(y[which(y$Method=="PLS"),measure], type = "p", main = "",pch=19, ylim=c(0,2.5),
         xlab="", ylab="", xaxt='n', yaxt="n", col="blue",cex=1.2)
    points(y[which(y$Method=="PLSs"),measure], type = "p", pch=19,col="red",cex=1.2)
    abline(v=c(1,2,3,4,5,6,7,8), col="gray", lty=3)
  }  
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  text((90+300)/2+25, (350+450)/2-25, 'Communality = 50%', cex=1.2, adj=0.5, font=2)
  #
  y=data[which(data$LV==var[1] & data$h=="h2"),]
  for(v in var[-length(var)])
  {
    y=data[which(data$LV==v & data$h=="h2"),]
    plot(y[which(y$Method=="PLS"),measure], type = "p", main = "",pch=19, ylim=c(0,2.5),
         xlab="", ylab="", xaxt='n', yaxt="n", col="blue",cex=1.2)
    points(y[which(y$Method=="PLSs"),measure], type = "p", pch=19,col="red",cex=1.2)
    abline(v=c(1,2,3,4,5,6,7,8), col="gray", lty=3)
  }
  y=data[which(data$LV==var[length(var)] & data$h=="h2"),]
  plot(y[which(y$Method=="PLS"),measure], type = "p", main = "",pch=15, ylim=c(0,2.5),
       xlab="", ylab="", xaxt='n', yaxt="n",col="blue",cex=1.2)
  points(y[which(y$Method=="PLSs"),measure], type = "p", pch=19,col="red",cex=1.2)
  abline(v=c(1,2,3,4,5,6,7,8), col="gray", lty=3)
  axis(4,at=c(0,0.5,1,1.5,2,2.5),cex.axis=0.8,font.axis=2, font.lab=2)
  plot(c(90, 350), c(300, 450), type = "n", main = "",
       xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
  text((90+300)/2+25, (350+450)/2-25, 'Communality = 75%', cex=1.2, adj=0.5, font=2)
  #
  y=data[which(data$LV==var[1] & data$h=="h3"),]
  plot(y[which(y$Method=="PLS"),measure], type = "p", main = "",pch=19, ylim=c(0,2.5),
       xlab="", ylab="", xaxt="n", yaxt="n",col="blue",cex=1.2)
  points(y[which(y$Method=="PLSs"),measure], type = "p", pch=19,col="red",cex=1.2)
  abline(v=c(1,2,3,4,5,6,7,8), col="gray", lty=3)
  axis(1,at=c(1,2,3,4,5,6,7,8),labels=x,las=2,cex.axis=0.8,font.axis=2, font.lab=2)
  axis(2,at=c(0,0.5,1,1.5,2.0,2.5),cex.axis=0.8,font.axis=2, font.lab=2)
  for(v in var[-1])
  {
    y=data[which(data$LV==v & data$h=="h3"),]
    plot(y[which(y$Method=="PLS"),measure], type = "p", main = "",pch=19, ylim=c(0,2.5),
         xlab="", ylab="", xaxt="n", yaxt="n", col="blue",cex=1.2)
    points(y[which(y$Method=="PLSs"),measure], type = "p", pch=19,col="red",cex=1.2)
    abline(v=c(1,2,3,4,5,6,7,8), col="gray", lty=3)
    axis(1,at=c(1,2,3,4,5,6,7,8),labels=x,las=2,cex.axis=0.8,font.axis=2, font.lab=2)
  }
  plot.new()
  par(mfrow=c(1,1))
}
###########################################################################################################################
### END OF plot.accuracy: plots RMSE and absolute bias for PLS and PLSs ###################################################
###########################################################################################################################



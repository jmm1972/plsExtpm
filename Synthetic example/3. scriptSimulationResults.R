###############################################################################################################################################
### This script processes the synthetic example Monte Carlo results producing article's plots and tables
###############################################################################################################################################
rm(list=ls())
#install.packages("mgcv")
library(mgcv)
#install.packages("lattice")
library(lattice)
#intall.packages ("reshape")
library(reshape)
######################################################################################################
## ALL LISTS ARE ORGANISED HIERARCHICALLY: Communality | Method | Sample size | Variable | Other object
######################################################################################################
source("your path/Synthetic example/scriptSimulationFunctions.R")
###########################################################################################################
## Convergence assessment and relationships estimation ####################################################
###########################################################################################################
results.converg = list()
samp.relshp = list(h1=list(),h2=list(),h3=list())
##
load("your_path/Synthetic example/Synthetic data generation/sim.pop.RData")
center=sim.pop.scale$center
scale=sim.pop.scale$scale
grid=sim.pop.scale$grid.xi
##
########################## Communality = 25% ##############################################################
load("your path/Synthetic example/Results/results.h1.RData")
results.converg[["h1"]] = run.converg.anal(h="h1",data=results.h1)
samp.relshp[["h1"]]=run.sample.relshp(h="h1",data=results.h1,grid=grid,index=results.converg[["h1"]]$index,center=center,scale=scale)
rm(list=c("results.h1"))
########################## Communality = 50% ##############################################################
load("your path/Synthetic example/Results/results.h2.RData")
results.converg[["h2"]] = run.converg.anal(h="h2",data=results.h2)
samp.relshp[["h2"]]=run.sample.relshp(h="h2",data=results.h2,grid=grid,index=results.converg[["h2"]]$index,center=center,scale=scale)
rm(list=c("results.h2"))
########################## Communality = 75% ##############################################################
load("your path/Synthetic example/Results/results.h3.RData")
results.converg[["h3"]] = run.converg.anal(h="h3",data=results.h3)
samp.relshp[["h3"]]=run.sample.relshp(h="h3",data=results.h3,grid=grid,index=results.converg[["h3"]]$index,center=center,scale=scale)
rm(list=c("results.h3"))
###
save(results.converg, file="your path/Synthetic example/Results/Convergence.results.RData")
save(samp.relshp,file = "your path/Synthetic example/Results/samp.relshp.RData")
###########################################################################################################
## END OF Convergence assessment and relationships estimation #############################################
###########################################################################################################
##
load("your path/Synthetic example/Results/Convergence.results.RData")
load("your path/Synthetic example/Results/samp.relshp.RData")
##
###########################################################################################################
## Convergence assessment plots ###########################################################################
###########################################################################################################
#
conv.summ.h1 = run.converg.summary(data=results.converg[["h1"]])
conv.summ.h2 = run.converg.summary(data=results.converg[["h2"]])
conv.summ.h3 = run.converg.summary(data=results.converg[["h3"]])
##
weig.neg.fn = function(data,h){
  weig.neg=matrix(0,9,2)
  rownames(weig.neg)=c("n75","n100","n150","n250","n300","n500","n600","n750","n900")
  colnames(weig.neg)=c("PLS", "PLSs")
  for (ss in c("n75","n100","n150","n250","n300","n500","n600","n750","n900")){
    weig.neg[ss,]=data[["weig.neg.freq"]][[h]][[ss]]/c(data[["cM"]][[h]][[ss]][1,1]+data[["cM"]][[h]][[ss]][1,2],data[["cM"]][[h]][[ss]][1,1]+data[["cM"]][[h]][[ss]][2,1])
  }
  rownames(weig.neg)=c(75,100,150,250,300,500,600,750,900)
  return(weig.neg)
}
Weig.neg.h1 = run.wght.neg(data=results.converg[["h1"]])
Weig.neg.h1 = Weig.neg.h1[which(row.names(Weig.neg.h1)!=600),]
Weig.neg.h2 = run.wght.neg(data=results.converg[["h2"]])
Weig.neg.h2 = Weig.neg.h2[which(row.names(Weig.neg.h2)!=600),]
Weig.neg.h3 = run.wght.neg(data=results.converg[["h3"]])
Weig.neg.h3 = Weig.neg.h3[which(row.names(Weig.neg.h3)!=600),]
##
old.par=par(no.readonly = TRUE)
layout(matrix(c(1,2,3),nrow=1))
par(mar=c(5.1,4,4.1,0))
barplot(t(Weig.neg.h1)*100,beside=T,
        xlab="", main="Communality = 25%",
        ylab="(%)",
        col=c("blue","red"))
#
barplot(t(Weig.neg.h2)*100,beside=T, cex.names = 1, cex.lab=0.75, ylim=c(0,100),
        xlab="Sample size", main="Communality = 50%", 
        ylab="",
        col=c("blue","red"), yaxt="n")
#
barplot(t(Weig.neg.h3)*100,beside=T, cex.names = 1, cex.lab=0.75, ylim=c(0,100),
        xlab="",  main="Communality = 75%",
        ylab="",
        col=c("blue","red"), yaxt="n")
legend("topright", c("PLS","PLSs"), col=c("blue","red"), bty="n", pch=c(15,15), cex=1.5)
par(old.par)
##
iter.med.h1 = run.mean.iter(data=results.converg[["h1"]])
iter.med.h1 = iter.med.h1[which(rownames(iter.med.h1)!=600),]
iter.med.h2 = run.mean.iter(data=results.converg[["h2"]])
iter.med.h2 = iter.med.h2[which(rownames(iter.med.h2)!=600),]
iter.med.h3 = run.mean.iter(data=results.converg[["h3"]])
iter.med.h3 = iter.med.h3[which(rownames(iter.med.h3)!=600),]
##
############
##
old.par=par(no.readonly = TRUE)
layout(matrix(c(1,2,3),nrow=1))
par(mar=c(5.1,4,4.1,0))
plot(rownames(iter.med.h1),iter.med.h1[,1], type="b", col="blue", 
     lty="dotdash", lwd=2, pch=15,
     ylim=c(min(iter.med.h1[,1],iter.med.h1[,2]),
            max(c(iter.med.h1[,1],iter.med.h1[,2]))),
     bty="n", xlab="", ylab="Iterations (no.)", main="Communality = 25%")
lines(rownames(iter.med.h1),iter.med.h1[,2], type="b", col="red", lty="dotdash", lwd=2, pch=16)
#
plot(rownames(iter.med.h2),iter.med.h2[,1], type="b", col="blue",
     lty="twodash", lwd=2, pch=15,
     ylim=c(min(iter.med.h2[,1],iter.med.h2[,2]),
            max(c(iter.med.h2[,1],iter.med.h2[,2]))),
     bty="n", xlab="Sample size", ylab="", main="Communality = 50%")
lines(rownames(iter.med.h2),iter.med.h2[,2], type="b", col="red", lty="twodash", lwd=2, pch=16)
#
plot(rownames(iter.med.h3),iter.med.h3[,1], type="b", col="blue",
     lty="twodash", lwd=2, pch=15,
     ylim=c(min(iter.med.h3[,1],iter.med.h3[,2]),
            max(c(iter.med.h3[,1],iter.med.h3[,2]))),
     bty="n", xlab="Sample size", ylab="", main="Communality = 75%")
lines(rownames(iter.med.h3),iter.med.h3[,2], type="b", col="red", lty=1, lwd=2, pch=16)
legend("topright", c("PLS","PLSs"), col=c("blue","red"), bty="n", pch=c(15,16), lty=c("dotdash","dotdash","twodash","twodash",1,1))
layout(matrix(1,nrow=1))
par(old.par)
##
###########################################################################################################
## END OF Convergence assessment plots ####################################################################
###########################################################################################################
##
##
#######################################################################################################
## PLOTS of impact ####################################################################################
#######################################################################################################
##
##
###################### Communality = 75% ##############################################################
plot.impact.fn(h="h3",n="n750",pop=sim.pop.latent,var="eta1",data=samp.relshp,grid=grid)
plot.impact.fn(h="h3",n="n750",pop=sim.pop.latent,var="eta2",data=samp.relshp,grid=grid)
plot.impact.fn(h="h3",n="n750",pop=sim.pop.latent,var="eta3",data=samp.relshp,grid=grid)
plot.impact.fn(h="h3",n="n750",pop=sim.pop.latent,var="eta4",data=samp.relshp,grid=grid)
plot.impact.fn(h="h3",n="n750",pop=sim.pop.latent,var="eta5",data=samp.relshp,grid=grid)
#
plot.impact.fn(h="h3",n="n600",pop=sim.pop.latent,var="eta1",data=samp.relshp,grid=grid)
plot.impact.fn(h="h3",n="n600",pop=sim.pop.latent,var="eta2",data=samp.relshp,grid=grid)
plot.impact.fn(h="h3",n="n600",pop=sim.pop.latent,var="eta3",data=samp.relshp,grid=grid)
plot.impact.fn(h="h3",n="n600",pop=sim.pop.latent,var="eta4",data=samp.relshp,grid=grid)
plot.impact.fn(h="h3",n="n600",pop=sim.pop.latent,var="eta5",data=samp.relshp,grid=grid)
######################################################################################################
############## ACCURACY ##############################################################################
######################################################################################################
##
mu=matrix(0,nrow=length(grid),ncol=5)
colnames(mu)=c("eta1","eta2","eta3","eta4","eta5")
mu[,1]<-4.85 + 0.5*grid-0.35*grid^2
mu[,2]<-5 - 0.5*grid^2
mu[,3]<-5 + 0.04*(1-exp(-1.5*grid))
mu[,4]<-5 + (0.1*grid - 0.3*grid^2)*(grid<0) + (0.1*grid - 0.1*grid^2)*(grid>0)
mu[,5]<-3.45 + 0.9*grid
##
accuracy=relshp.accuracy(data=samp.relshp, mu=mu, index=results.converg, grid=grid)
##
###################################################################################
acc=data.frame(Sample=rep(NA,3*2*5*9),
               LV=rep(NA,3*2*5*9),
               Method=rep(NA,3*2*5*9),
               h=rep(NA,3*2*5*9),
               RMSE=rep(NA,3*2*5*9),
               Bias=rep(NA,3*2*5*9))
count=0
for (h in c("h1","h2","h3"))
{
  for (m in c("PLS","PLSs"))
  {
    for (v in c("eta1","eta2","eta3","eta4","eta5"))
    {
      for (ss in c("n75","n100","n150","n250","n300","n500","n600","n750","n900"))
      {
        count=count+1
        cat(paste0("Row: ",count, " - h: ", h," - Method: ", m, " - var: ", v, " - ss: ", ss,"\n"))
        acc$h[count]=h
        acc$Method[count]=m
        acc$LV[count]=v
        #acc$Sample[count]=as.integer(substr(ss,2,1000000L))
        acc$Sample[count]=ss
        acc$RMSE[count]=accuracy[["RMSE"]][[h]][[m]][[ss]][[v]]
        acc$Bias[count]=accuracy[["B"]][[h]][[m]][[ss]][[v]]
      }
    }
  }
}
#
write.table(acc,"your path/Synthetic example/Synthetic data generation/Results/accuracy.txt")
##
acc=read.table("your path/Synthetic example/Synthetic data generation/Results/accuracy.txt",header=T)
acc = acc[acc$Sample!="n600",]
##
dotplot(Bias~Sample|LV+h,data=acc,groups = Method,subscripts=T,
        xlab = "Sample size",ylab="Absolute bias",scales=list(rot=c(90,90),cex=0.6), pch=c(16,16),col=c("red","blue"))
##
##
dotplot(RMSE~Sample|LV+h,data=acc,groups = Method,subscripts=T,
        xlab = "Sample size",ylab="Absolute bias",scales=list(rot=c(90,90),cex=0.6), pch=c(16,16),col=c("red","blue"))

######################################################################################################################
#####################################
# Index number for Bias and RMSE
index.B = data.frame(Sample=subset(acc,acc$Method=="PLS")$Sample,LV=subset(acc,acc$Method=="PLS")$LV, h=subset(acc,acc$Method=="PLS")$h,
                     Index=subset(acc,acc$Method=="PLSs")$Bias/subset(acc,acc$Method=="PLS")$Bias*100)
#
index.RMSE<-data.frame(Sample=subset(acc,acc$Method=="PLS")$Sample,LV=subset(acc,acc$Method=="PLS")$LV, h=subset(acc,acc$Method=="PLS")$h,
                       Index=subset(acc,acc$Method=="PLSs")$RMSE/subset(acc,acc$Method=="PLS")$RMSE*100)
#
index<-rbind(index.B,index.RMSE)
index$Measure=c(rep("Bias",nrow(index.B)),rep("RMSE",nrow(index.RMSE)))
index$Measure<-factor(index$Measure)
dotplot(Index~Sample|LV+h,data=index,groups = Measure,xlab = "Sample size",
        ylab="Absolute bias/Root mean square error",
        scales=list(rot=c(90,90),cex=0.6),
        pch=c(16,16),col=c("green","yellow"))
#
###############################################################################
################################ ANOVA ########################################
###############################################################################
#
acc1=acc
acc1$Method<-as.factor(acc1$Method)
acc1$Sample<-as.factor(acc1$Sample)
levels(acc$Sample)=c("75","100","150","250","300","500","750","900")
acc1$LV=as.factor(acc1$LV)
acc1$h<-as.factor(acc1$h)
levels(acc1$h)<-c("25%","50%","75%")

parameter=c("$n=100$","$n=150$","$n=250$","$n=300$","$n=500$","$n=750$","$n=900$",
            "PLSs",
            "h=50\\%","h=75\\%",
            "PLSs|h=50\\%",
            "PLSs|h=75\\%",
            "PLSs|$n=100$","PLSs|$n=150$","PLSs|$n=250$","PLSs|$n=300$","PLSs|$n=500$","PLSs|$n=750$","PLSs|$n=900$")
B.eta1.coeff=round(summary(lm(Bias~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta1")))$coefficients[-1,1],3)
B.eta2.coeff=round(summary(lm(Bias~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta2")))$coefficients[-1,1],3)
B.eta3.coeff=round(summary(lm(Bias~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta3")))$coefficients[-1,1],3)
B.eta4.coeff=round(summary(lm(Bias~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta4")))$coefficients[-1,1],3)
B.eta5.coeff=round(summary(lm(Bias~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta5")))$coefficients[-1,1],3)
B.eta1.pv=round(summary(lm(Bias~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta1")))$coefficients[-1,1],3)
B.eta2.pv=round(summary(lm(Bias~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta2")))$coefficients[-1,1],3)
B.eta3.pv=round(summary(lm(Bias~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta3")))$coefficients[-1,1],3)
B.eta4.pv=round(summary(lm(Bias~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta4")))$coefficients[-1,1],3)
B.eta5.pv=round(summary(lm(Bias~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta5")))$coefficients[-1,1],3)

ANOVA.Bias=data.frame(parameter=parameter,eta1=rep(NA,19),eta2=rep(NA,19),eta3=rep(NA,19),eta4=rep(NA,19),eta5=rep(NA,19))
for (i in 1:length(summary(lm(Bias~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta1")))$coefficients[-1,1])){
  ANOVA.Bias[i,"eta1"]=paste0(B.eta1.coeff[i], "(",B.eta1.pv[i],")")    
  ANOVA.Bias[i,"eta2"]=paste0(B.eta2.coeff[i], "(",B.eta2.pv[i],")")    
  ANOVA.Bias[i,"eta3"]=paste0(B.eta3.coeff[i], "(",B.eta3.pv[i],")")    
  ANOVA.Bias[i,"eta4"]=paste0(B.eta4.coeff[i], "(",B.eta4.pv[i],")")    
  ANOVA.Bias[i,"eta5"]=paste0(B.eta5.coeff[i], "(",B.eta5.pv[i],")")    
}

RMSE.eta1.coeff=round(summary(lm(RMSE~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta1")))$coefficients[-1,1],3)
RMSE.eta2.coeff=round(summary(lm(RMSE~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta2")))$coefficients[-1,1],3)
RMSE.eta3.coeff=round(summary(lm(RMSE~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta3")))$coefficients[-1,1],3)
RMSE.eta4.coeff=round(summary(lm(RMSE~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta4")))$coefficients[-1,1],3)
RMSE.eta5.coeff=round(summary(lm(RMSE~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta5")))$coefficients[-1,1],3)
RMSE.eta1.pv=round(summary(lm(RMSE~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta1")))$coefficients[-1,1],3)
RMSE.eta2.pv=round(summary(lm(RMSE~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta2")))$coefficients[-1,1],3)
RMSE.eta3.pv=round(summary(lm(RMSE~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta3")))$coefficients[-1,1],3)
RMSE.eta4.pv=round(summary(lm(RMSE~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta4")))$coefficients[-1,1],3)
RMSE.eta5.pv=round(summary(lm(RMSE~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta5")))$coefficients[-1,1],3)

ANOVA.RMSE=data.frame(parameter=parameter,eta1=rep(NA,19),eta2=rep(NA,19),eta3=rep(NA,19),eta4=rep(NA,19),eta5=rep(NA,19))
for (i in 1:length(summary(lm(RMSE~Sample+Method+h+Method*h+Method*Sample,data=subset(acc1,LV=="eta1")))$coefficients[-1,1])){
  ANOVA.RMSE[i,"eta1"]=paste0(RMSE.eta1.coeff[i], "(",RMSE.eta1.pv[i],")")    
  ANOVA.RMSE[i,"eta2"]=paste0(RMSE.eta2.coeff[i], "(",RMSE.eta2.pv[i],")")    
  ANOVA.RMSE[i,"eta3"]=paste0(RMSE.eta3.coeff[i], "(",RMSE.eta3.pv[i],")")    
  ANOVA.RMSE[i,"eta4"]=paste0(RMSE.eta4.coeff[i], "(",RMSE.eta4.pv[i],")")    
  ANOVA.RMSE[i,"eta5"]=paste0(RMSE.eta5.coeff[i], "(",RMSE.eta5.pv[i],")")    
}


write.table(cbind(summary(lm(RMSE~Sample+Method+h+Method*h+Method*Sample,data=subset(acc,LV=="eta1")))$coefficients[,c(1,4)],
summary(lm(RMSE~Sample+Method+h+Method*h+Method*Sample,data=subset(acc,LV=="eta2")))$coefficients[,c(1,4)],
summary(lm(RMSE~Sample+Method+h+Method*h+Method*Sample,data=subset(acc,LV=="eta3")))$coefficients[,c(1,4)],
summary(lm(RMSE~Sample+Method+h+Method*h+Method*Sample,data=subset(acc,LV=="eta4")))$coefficients[,c(1,4)],
summary(lm(RMSE~Sample+Method+h+Method*h+Method*Sample,data=subset(acc,LV=="eta5")))$coefficients[,c(1,4)]),
"your path/Synthetic example/Synthetic data generation/Results/ANOVA_RMSE.csv",sep=";")






















###################################################################################################################################################################
### This scripts allows running plspm , seminR and and plsExtpm to compare the results using scaled manifest variables
### At the end of the scrpt the code using the non-scaled manifest variables that was used in the artice is reproduced
###################################################################################################################################################################
rm(list=ls())
##
### Reading the necessary libraries
##
#install.packages("plspm")
library(plspm)
#install.packages("seminr")
library(seminr)
#install.packages("ltm")
library(ltm)
##
###################################################################################################################################################################
### Reading bank data
ECSIbank_data<-read.table("your_path/ECSI example/ECSI.txt",header=TRUE)
##
###################################################################################################################################################################
### plspm data
# structural model
##
Quality = c(0,0,0,0)
Value = c(1,0,0,0)
Satisfaction = c(1,1,0,0)
Loyalty = c(0,0,1,0)
sat_path = rbind(Quality, Value, Satisfaction, Loyalty)
colnames(sat_path)=c("Quality", "Value", "Satisfaction", "Loyalty")
##
# measurement model
sat_blocks = list(1:9, 10:11, 12:14, 15:16)
##
# vector of modes (reflective indicators)
##
sat_mod = c("B",rep("A", 3))
##
###################################################################################################################################################################
### plspm run
ECSI.plspm.comp = plspm(ECSIbank_data, sat_path, sat_blocks, modes = sat_mod)
ECSI.plspm.comp
##
##################################
#Partial Least Squares Path Modeling (PLS-PM) 
#---------------------------------------------
#  NAME             DESCRIPTION
#1  $outer_model     outer model
#2  $inner_model     inner model
#3  $path_coefs      path coefficients matrix
#4  $scores          latent variable scores
#5  $crossloadings   cross-loadings
#6  $inner_summary   summary inner model
#7  $effects         total effects
#8  $unidim          unidimensionality
#9  $gof             goodness-of-fit
#10 $boot            bootstrap results
#11 $data            data matrix
#---------------------------------------------
#  You can also use the function 'summary'
##
summary(ECSI.plspm.comp)
##
############################################
ECSI.plspm.summary = list(inner = data.frame(Type=ECSI.plspm.comp$inner_summary$Type,
                                             Mode=ECSI.plspm.comp$unidim$Mode,
                                             MVs=ECSI.plspm.comp$unidim$MVs,
                                             C.alpha=ECSI.plspm.comp$unidim$C.alpha,
                                             DG.rho=ECSI.plspm.comp$unidim$DG.rho,
                                             R2=ECSI.plspm.comp$inner_summary$R2,
                                             Blck_comm=ECSI.plspm.comp$inner_summary$Block_Communality,
                                             row.names=row.names(ECSI.plspm.comp$unidim)),
                          outer = ECSI.plspm.comp$outer_model[,1:5],
                          tot.effs = ECSI.plspm.comp$effects,
                          path.coefs=t(ECSI.plspm.comp$path_coefs))
##
ECSI.plspm.summary$outer$loading[1:9]=NA
ECSI.plspm.summary$outer$communality[1:9]=NA
ECSI.plspm.summary
##
###################################################################################################################################################################
### seminr data
# measurement model
mm.m = constructs(
  composite("Quality", multi_items("QUAL", 1:9),weights= mode_B),
  composite("Value", multi_items("VALU", 1:2),weights= mode_A),
  composite("Satisfaction", multi_items("SATI", 1:3),weights= mode_A),
  composite("Loyalty", multi_items("LOYA", c(1,3)),weights= mode_A)
)
# structural model
sm.m= relationships(
  paths(from="Quality", to = c("Value","Satisfaction")),
  paths(from="Satisfaction", to = c("Loyalty")),
  paths(from= "Value", to=c("Satisfaction"))
)
ECSI.seminr.comp=estimate_pls(data=ECSIbank_data, measurement_model=mm.m, structural_model=sm.m)
##
ECSI.seminr.summary = list(inner = data.frame(Type=c("Exogenous",rep("Endogenous",3)),
                                             Mode=rep("A",4),
                                             MVs=c(9,2,3,2),
                                             C.alpha=summ$reliability[,1],
                                             DG.rho=summ$reliability[,2],
                                             R2=c(0,summ$paths[1,]),
                                             row.names=c("Quality","Value","Satisfaction","Loyalty")),
                          outer = cbind(weight=rowSums(summ$weights),loading=rowSums(summ$loadings)),
                          tot.effs = summ$total_effects,
                          path.coefs=summ$paths[-c(1,2),])
ECSI.seminr.summary$outer[1:9,"loading"]=NA
##################################
### plsExtpm SEM model ###
## measurement model
ECSImm<-read.table("your_path/ECSI example/ECSImm.txt",header=TRUE)
##
## structural model
ECSIsm_all_linear<-as.data.frame(read.csv("your_path/ECSI example/ECSIsm_all_linear.csv",header=TRUE))
##
ECSIsm_all_nlinear<-as.data.frame(read.csv("your_path/ECSI example/ECSIsm_all_nlinear.csv",header=TRUE))
##
ECSIsm_all_nlinear<-as.data.frame(read.csv("your_path/ECSI example/ECSIsm_all_nlinear.csv", header=TRUE))
##
## Reading plsEXTpm package functions
source("your_path/package/plsExtpm.R")
##
ECSI.plsExtpm.comp = plsExtpm(data = ECSIbank_data, strucmod = ECSIsm_all_linear, measuremod = ECSImm)
##
ECSI.plsExtpm.comp.boot = plsExtboot(pls.object = ECSI.plsExtpm.comp, nboot=500, start="old")
##
plsExtplot(pls.object = ECSI.plsExtpm.comp, boot.object = ECSI.plsExtpm.comp.boot, x.LV="Quality", y.LV = "Value")
plsExtplot(pls.object = ECSI.plsExtpm.comp, boot.object = ECSI.plsExtpm.comp.boot, x.LV="Quality", y.LV = "Satisfaction")
plsExtplot(pls.object = ECSI.plsExtpm.comp, boot.object = ECSI.plsExtpm.comp.boot, x.LV="Value", y.LV = "Satisfaction")
plsExtplot(pls.object = ECSI.plsExtpm.comp, boot.object = ECSI.plsExtpm.comp.boot, x.LV="Satisfaction", y.LV = "Loyalty")
##
ECSI.plsExtpm.comp.all_nlinear = plsExtpm(data = ECSIbank_data, strucmod = ECSIsm_all_nlinear, measuremod = ECSImm)
##
ECSI.plsExtpm.comp.all_nlinear.boot = plsExtboot(pls.object = ECSI.plsExtpm.comp.all_nlinear, nboot=500, start="old")
##
plsExtplot(pls.object = ECSI.plsExtpm.comp.all_nlinear, boot.object = ECSI.plsExtpm.comp.all_nlinear.boot, x.LV="Quality", y.LV = "Value")
plsExtplot(pls.object = ECSI.plsExtpm.comp.all_nlinear, boot.object = ECSI.plsExtpm.comp.all_nlinear.boot, x.LV="Quality", y.LV = "Satisfaction")
plsExtplot(pls.object = ECSI.plsExtpm.comp.all_nlinear, boot.object = ECSI.plsExtpm.comp.all_nlinear.boot, x.LV="Value", y.LV = "Satisfaction")
plsExtplot(pls.object = ECSI.plsExtpm.comp.all_nlinear, boot.object = ECSI.plsExtpm.comp.all_nlinear.boot, x.LV="Satisfaction", y.LV = "Loyalty")
##
strucmod_SAT_LOYA_nlinear = ECSIsm_all_nlinear
strucmod_SAT_LOYA_nlinear$type = c("ln","ln","ln","nln")
strucmod_SAT_LOYA_nlinear$K = c(NA,NA,NA,10)
##
ECSI.plsExtpm.comp.SATI_LOYA_nlinear = plsExtpm(data = ECSIbank_data, strucmod = strucmod_SAT_LOYA_nlinear, measuremod = ECSImm)
##
ECSI.plsExtpm.comp.SATI_LOYA_nlinear.boot = plsExtboot(pls.object = ECSI.plsExtpm.comp.SATI_LOYA_nlinear, nboot=5, start="old")
##
plsExtplot(pls.object = ECSI.plsExtpm.comp.SATI_LOYA_nlinear, boot.object = ECSI.plsExtpm.comp.SATI_LOYA_nlinear.boot, x.LV="Quality", y.LV = "Value")
plsExtplot(pls.object = ECSI.plsExtpm.comp.SATI_LOYA_nlinear, boot.object = ECSI.plsExtpm.comp.SATI_LOYA_nlinear.boot, x.LV="Quality", y.LV = "Satisfaction")
plsExtplot(pls.object = ECSI.plsExtpm.comp.SATI_LOYA_nlinear, boot.object = ECSI.plsExtpm.comp.SATI_LOYA_nlinear.boot, x.LV="Value", y.LV = "Satisfaction")
plsExtplot(pls.object = ECSI.plsExtpm.comp.SATI_LOYA_nlinear, boot.object = ECSI.plsExtpm.comp.SATI_LOYA_nlinear.boot, x.LV="Satisfaction", y.LV = "Loyalty")
##
plsExtcheck(pls.object = ECSI.plsExtpm.comp.SATI_LOYA_nlinear)
##
ECSI.plsExtpm.summary = list(inner = data.frame(LVs=c(exogenous(ECSI.plsExtpm.comp),endogenous(ECSI.plsExtpm.comp)),
                                              Type=c(rep("exogenous",length(exogenous(ECSI.plsExtpm.comp))),rep("endogenous",length(endogenous(ECSI.plsExtpm.comp)))),
                                              Mode=str.mode(ECSI.plsExtpm.comp),
                                              MVs=dgrho(object=ECSI.plsExtpm.comp)[,2],
                                              C.alpha=c(cronbach.alpha(ECSIbank_data[,1:9])$alpha,
                                                        cronbach.alpha(ECSIbank_data[,10:11])$alpha,
                                                        cronbach.alpha(ECSIbank_data[,12:14])$alpha,
                                                        cronbach.alpha(ECSIbank_data[,15:16])$alpha),
                                              DG.rho=dgrho(object=ECSI.plsExtpm.comp)[,1],
                                              R2.all_linear=ECSI.plsExtpm.comp$R.Squared[,1],
                                              R2.all_nlinear=ECSI.plsExtpm.comp.all_nlinear$R.Squared[,1],
                                              R2.SATI_LOYA_nlinear=ECSI.plsExtpm.comp.SATI_LOYA_nlinear$R.Squared[,1],
                                              Blck_comm=communality(ECSI.plsExtpm.comp)[,1]),
                           outer = data.frame(name=ECSI.plsExtpm.comp$model$manifest,
                                              block=rep(names(ECSI.plsExtpm.comp$model$blocks),noindicators(ECSI.plsExtpm.comp)),
                                              weight=rowSums(ECSI.plsExtpm.comp$outer_weights),
                                              loading=rowSums(ECSI.plsExtpm.comp$outer_loadings),
                                              communality=rowSums(ECSI.plsExtpm.comp$outer_loadings)^2),
                           tot.effs = ECSI.plsExtpm.comp$total_effects$totEff.summary,
                           path.coefs=ECSI.plsExtpm.comp$path_coefficients$pC)
##
###################################################################################################################################################################
### Here is the code to reproduce de article results
###################################################################################################################################################################
##
ECSI.plsExtpm.all_linear2 <- plsExtpm(data = ECSIbank_data, strucmod = ECSIsm_all_linear, measuremod = ECSImm, scaled=FALSE, sum1=TRUE)
##
ECSI.plsExtpm.all_nlinear2 <- plsExtpm(data = ECSIbank_data, strucmod = ECSIsm_all_nlinear, measuremod = ECSImm, scaled=FALSE, sum1=TRUE)
##
ECSI.plsExtpm.all_nlinear.boot2 = plsExtboot(pls.object = ECSI.plsExtpm.all_nlinear2, nboot=500, start="old")
##
plsExtplot(pls.object = ECSI.plsExtpm.all_nlinear2, boot.object = ECSI.plsExtpm.all_nlinear.boot2, x.LV="Quality", y.LV = "Value", xlim=c(0,10), ylim=c(0,10))
plsExtplot(pls.object = ECSI.plsExtpm.all_nlinear2, boot.object = ECSI.plsExtpm.all_nlinear.boot2, x.LV="Quality", y.LV = "Satisfaction", xlim=c(0,10), ylim=c(0,10))
plsExtplot(pls.object = ECSI.plsExtpm.all_nlinear2, boot.object = ECSI.plsExtpm.all_nlinear.boot2, x.LV="Value", y.LV = "Satisfaction", xlim=c(0,10), ylim=c(0,10))
plsExtplot(pls.object = ECSI.plsExtpm.all_nlinear2, boot.object = ECSI.plsExtpm.all_nlinear.boot2, x.LV="Satisfaction", y.LV = "Loyalty", xlim=c(0,10), ylim=c(0,10))
##
strucmod_SAT_LOYA_nlinear = ECSIsm_all_nlinear
strucmod_SAT_LOYA_nlinear$type = c("ln","ln","ln","nln")
strucmod_SAT_LOYA_nlinear$K = c(NA,NA,NA,10)
ECSI.plsExtpm.SAT_LOYA_nlinear <- plsExtpm(data = ECSIbank_data, strucmod = strucmod_SAT_LOYA_nlinear, measuremod = ECSImm, scaled=FALSE, sum1=TRUE)
##
ECSI.plsExtpm.SAT_LOYA_nlinear.boot = plsExtboot(pls.object = ECSI.plsExtpm.SAT_LOYA_nlinear, nboot=500, start="old")
##
plsExtplot(pls.object = ECSI.plsExtpm.SAT_LOYA_nlinear, boot.object = ECSI.plsExtpm.SAT_LOYA_nlinear.boot, x.LV="Quality", y.LV = "Value", xlim=c(0,10), ylim=c(0,10))
plsExtplot(pls.object = ECSI.plsExtpm.SAT_LOYA_nlinear, boot.object = ECSI.plsExtpm.SAT_LOYA_nlinear.boot, x.LV="Quality", y.LV = "Satisfaction", xlim=c(0,10), ylim=c(0,10))
plsExtplot(pls.object = ECSI.plsExtpm.SAT_LOYA_nlinear, boot.object = ECSI.plsExtpm.SAT_LOYA_nlinear.boot, x.LV="Value", y.LV = "Satisfaction", xlim=c(0,10), ylim=c(0,10))
plsExtplot(pls.object = ECSI.plsExtpm.SAT_LOYA_nlinear, boot.object = ECSI.plsExtpm.SAT_LOYA_nlinear.boot, x.LV="Satisfaction", y.LV = "Loyalty", xlim=c(0,10), ylim=c(0,10))
##
###################################################################################################################################################################
### END OF ECSI EXAMPLE ###########################################################################################################################################
###################################################################################################################################################################

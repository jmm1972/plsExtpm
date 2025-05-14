###########################################################################################################
### This script runs the Monte Carlo process
###########################################################################################################
rm(list=ls())
###
## Set your working directory to store results
##
source("your_path/package/plsExtpm.R")
source("your path/Synthetic example/scriptSimulationFunctions.R")
##
## setting measurement model
mm.m<-as.matrix(cbind(c(rep("xi",5),rep("eta1",5),rep("eta2",5),rep("eta3",5),rep("eta4",5),rep("eta5",5)),
                      c(c("xi.1","xi.2","xi.3","xi.4","xi.5"),
                        c("eta1.1","eta1.2","eta1.3","eta1.4","eta1.5"),
                        c("eta2.1","eta2.2","eta2.3","eta2.4","eta2.5"),
                        c("eta3.1","eta3.2","eta3.3","eta3.4","eta3.5"),
                        c("eta4.1","eta4.2","eta4.3","eta4.4","eta4.5"),
                        c("eta5.1","eta5.2","eta5.3","eta5.4","eta5.5"))))
colnames(mm.m)<-c("source","target")
##
## setting structural model
sm.m<-as.matrix(cbind(c(rep("xi",5)),c("eta1","eta2","eta3","eta4","eta5")))
colnames(sm.m)<-c("source","target")
#
load('your path/Synthetic example/Synthetic data generation/sim.pop.RData')

load('your path/Synthetic example/Synthetic data generation/sample.index.RData')
#
###########################################################################################################
### Results ###############################################################################################
###########################################################################################################
results.h3=run.results(h="h3",pop=sim.pop.h3,sample.idx=sample.index)
save(results.h3,file="your path/Synthetic example/Synthetic data generation/Results/results.h3.RData")
rm(results.h3)
###########################
results.h2=run.results(h="h2",pop=sim.pop.h2,sample.idx=sample.index)
save(results.h2,file="your path/Synthetic example/Synthetic data generation/Results/results.h2.RData")
rm(results.h2)
###########################
results.h1=run.results(h="h1",pop=sim.pop.h1,sample.idx=sample.index)
save(results.h1,file="your path/Synthetic example/Synthetic data generation/Results/results.h1.RData")
rm(results.h1)
##
###########################################################################################################
## END OF Monte Carlo results #############################################################################
###########################################################################################################

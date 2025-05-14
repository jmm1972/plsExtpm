###################################################################################################################################################################
### This scripts was used to generate the synthetic data used in the article
###################################################################################################################################################################
rm(list=ls())
#installpackage("MASS")
library(MASS)
n=10000
var_h1=3 # 1^2/(1^2+3)=0.25 (Communality = 25%)
var_h2=1  # 1^2/(1^2+1)=0.50 (Communality = 50%)
var_h3=1/3  # 1^2/(1^2+1/3)=0.75 (Communality = 75%)
#
xi = sort(rnorm(n))
eta1f = 4.85 + 0.5*xi - 0.35*xi^2
eta2f = 5 - 0.5*xi^2
eta3f = 5 + 0.04*(1-exp(-1.5*xi))
eta4f = 5 + (0.1*xi - 0.3*xi^2)*(xi<0) + (0.1*xi - 0.1*xi^2)*(xi>0)
eta5f = 3.45 + 0.9*xi
#
plot(xi,eta1f, type="l", ylim=c(-5,8))
lines(xi,eta2f)
lines(xi,eta3f)
lines(xi,eta4f)
lines(xi,eta5f)
#
S = matrix(0,5,5); diag(S)=rep(0.5,5)
epsilonLV=mvrnorm(n,mu=rep(0,5),Sigma=S)
#
eta1 = eta1f + epsilonLV[,1]
eta2 = eta2f + epsilonLV[,2]
eta3 = eta3f + epsilonLV[,3]
eta4 = eta4f + epsilonLV[,4]
eta5 = eta5f + epsilonLV[,5]
cor(xi,eta1)^2
cor(xi,eta2)^2
cor(xi,eta3)^2
cor(xi,eta4)^2
cor(xi,eta5)^2
##
sim.pop.latent=data.frame(xi=xi, eta1=eta1,eta2=eta2,eta3=eta3,eta4=eta4,eta5=eta5,
                          eta1f=eta1f,eta2f=eta2f,eta3f=eta3f,eta4f=eta4f,eta5f=eta5f)
grid.xi<-seq(from=min(sim.pop.latent$xi),to=max(sim.pop.latent$xi),length.out=300)
########### replace bt population theoretical values
## this values are documented in the article's suplemantary material
##
center=c(0,4.500000,4.500000,4.916791,4.800000,3.450000)
scale=c(1,0.9974969,0.9974969,0.7929991,0.7810250,1.1445523)
names(center) = c("xi","eta1","eta2","eta3","eta4","eta5")
names(scale) = c("xi","eta1","eta2","eta3","eta4","eta5")
##
##########################
sim.pop.scale=list(grid.xi=grid.xi,center=center,scale=scale)
sim.pop.scale$center=center
sim.pop.scale$scale=scale

#################################################################################################
gen.MV = function(n,LV,var_h,names){
  S = matrix(0,5,5);diag(S)=rep(var_h,5)
  epsilon1 = mvrnorm(n,mu=rep(0,5),Sigma=S)
  epsilon2 = mvrnorm(n,mu=rep(0,5),Sigma=S)
  epsilon3 = mvrnorm(n,mu=rep(0,5),Sigma=S)
  epsilon4 = mvrnorm(n,mu=rep(0,5),Sigma=S)
  epsilon5 = mvrnorm(n,mu=rep(0,5),Sigma=S)
  epsilon6 = mvrnorm(n,mu=rep(0,5),Sigma=S)
  #dim(epsilon)
  LV1 = data.frame(LV[,1],LV[,1],LV[,1],LV[,1],LV[,1],
                   LV[,2],LV[,2],LV[,2],LV[,2],LV[,2],
                   LV[,3],LV[,3],LV[,3],LV[,3],LV[,3],
                   LV[,4],LV[,4],LV[,4],LV[,4],LV[,4],
                   LV[,5],LV[,5],LV[,5],LV[,5],LV[,5],
                   LV[,6],LV[,6],LV[,6],LV[,6],LV[,6])
  head(LV1)
  #dim(LV1)
  colnames(LV1)
  x1 = LV1[,1:5] + epsilon1
  x2 = LV1[,6:10] + epsilon2
  x3 = LV1[,11:15] + epsilon3
  x4 = LV1[,16:20] + epsilon4
  x5 = LV1[,21:25] + epsilon5
  x6 = LV1[,26:30] + epsilon6
  x=data.frame(x1,x2,x3,x4,x5,x6)
  colnames(x)=c(paste0(names[1],1:5),
                paste0(names[2],1:5),
                paste0(names[3],1:5),
                paste0(names[4],1:5),
                paste0(names[5],1:5),
                paste0(names[6],1:5))
  return(x)
}
#################################################################################################
# communality = 25%
#################################################################################################
sim.pop.h1 = gen.MV(n=n,LV=data.frame(xi=xi,eta1=eta1,eta2=eta2,eta3=eta3,eta4=eta4,eta5=eta5),
                    var_h=var_h1,
                    names=c("xi.","eta1.","eta2.","eta3.","eta4.","eta5."))
##
library(ltm)
cronbach.alpha(sim.pop.h1[,1:5])$alpha
cronbach.alpha(sim.pop.h1[,6:10])$alpha
cronbach.alpha(sim.pop.h1[,11:15])$alpha
cronbach.alpha(sim.pop.h1[,16:20])$alpha
cronbach.alpha(sim.pop.h1[,21:25])$alpha
cronbach.alpha(sim.pop.h1[,26:30])$alpha
#################################################################################################
# communality = 50%
#################################################################################################
sim.pop.h2 = gen.MV(n=n,LV=data.frame(xi=xi,eta1=eta1,eta2=eta2,eta3=eta3,eta4=eta4,eta5=eta5),
                    var_h=var_h2,
                    names=c("xi.","eta1.","eta2.","eta3.","eta4.","eta5."))
##
cronbach.alpha(sim.pop.h2[,1:5])$alpha
cronbach.alpha(sim.pop.h2[,6:10])$alpha
cronbach.alpha(sim.pop.h2[,11:15])$alpha
cronbach.alpha(sim.pop.h2[,16:20])$alpha
cronbach.alpha(sim.pop.h2[,21:25])$alpha
cronbach.alpha(sim.pop.h2[,26:30])$alpha
#################################################################################################
# communality = 75%
#################################################################################################
sim.pop.h3 = gen.MV(n=n,LV=data.frame(xi=xi,eta1=eta1,eta2=eta2,eta3=eta3,eta4=eta4,eta5=eta5),
                    var_h=var_h3,
                    names=c("xi.","eta1.","eta2.","eta3.","eta4.","eta5."))
##
cronbach.alpha(sim.pop.h3[,1:5])$alpha
cronbach.alpha(sim.pop.h3[,6:10])$alpha
cronbach.alpha(sim.pop.h3[,11:15])$alpha
cronbach.alpha(sim.pop.h3[,16:20])$alpha
cronbach.alpha(sim.pop.h3[,21:25])$alpha
cronbach.alpha(sim.pop.h3[,26:30])$alpha
#
save(sim.pop.latent,
     sim.pop.h1,
     sim.pop.h2,
     sim.pop.h3,
     sim.pop.scale,
          file="your_path/Synthetic example/Synthetic data generation/sim.pop.RData")

##
#####################################################################################################
plot(xi,eta1f, type="l")
points(xi,sim.pop.h1$eta2.1, col=2)
points(xi,sim.pop.h1$eta2.2, col=3)
points(xi,sim.pop.h1$eta2.3, col=4)
points(xi,sim.pop.h1$eta2.4, col=5)
points(xi,sim.pop.h1$eta2.5, col=6)
#
plot(xi,eta1f, type="l")
points(xi,sim.pop.h2$eta2.1, col=2)
points(xi,sim.pop.h2$eta2.2, col=3)
points(xi,sim.pop.h2$eta2.3, col=4)
points(xi,sim.pop.h2$eta2.4, col=5)
points(xi,sim.pop.h2$eta2.5, col=6)
#
plot(xi,eta1f, type="l")
points(xi,sim.pop.h3$eta2.1, col=2)
points(xi,sim.pop.h3$eta2.2, col=3)
points(xi,sim.pop.h3$eta2.3, col=4)
points(xi,sim.pop.h3$eta2.4, col=5)
points(xi,sim.pop.h3$eta2.5, col=6)
#
plot(xi,eta4f, type="l")
points(xi,sim.pop.h3$eta5.1, col=2)
points(xi,sim.pop.h3$eta5.2, col=3)
points(xi,sim.pop.h3$eta5.3, col=4)
points(xi,sim.pop.h3$eta5.4, col=5)
points(xi,sim.pop.h3$eta5.5, col=6)
##
#####################################################################################################
### SAMPLES
#####################################################################################################
rm(list=ls())
n=10000
load("your_path/Synthetic example/Synthetic data generation/sim_pop.RData")
#
samp.sel = function(size){
  data = list()
  for (ss in 1:length(size))
  {
    for (s in 1:1000)
    {
      index.h1 = sample(1:n,size[ss])
      index.h2 = sample(1:n,size[ss])
      index.h3 = sample(1:n,size[ss])
      data[["h1"]][[paste0("n",as.character(size[ss]))]][[s]] = index.h1
      data[["h2"]][[paste0("n",as.character(size[ss]))]][[s]] = index.h2
      data[["h3"]][[paste0("n",as.character(size[ss]))]][[s]] = index.h3
    }
  }
  return(data)
}
########################
sample.index = samp.sel(size=c(75,100,150,250,300,500,600,750,900))
save(sample.index,file="your_path/Synthetic example/Synthetic data generation/sample.index.RData")

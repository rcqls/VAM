# VAM package instalation
#devtools::install_github("rcqls/VAM")

require(VAM)
#
#####################
# 1 trajectory
# only CM
#####################
#
##########
# Simulation model:
# ARA infinit CM with \rho=0.4 
# Weibull h(t)=0.001*2.5*t^(1.5)
simARAInf<-sim.vam(  ~ (ARAInf(.4) | Weibull(.001,2.5)))
#
n<-30 #number of simulated maintenance
simData<-simulate(simARAInf,n)
head(simData)
tail(simData)
simData<-simulate(simARAInf,n)#another simulation
tail(simData)
#
##########
# Model plots
# in order to plot the corresponding intensity
plot(simARAInf,'i')
# or virtual age
plot(simARAInf,'v')
# or cumulative intensity
plot(simARAInf,'I')
#
# The model associate to simData can be changed
modelAGAN<-model.vam(Time & Type ~ (AGAN() | Weibull(.001,2.5)),data=simData)
# in order to plot
plot(modelAGAN,'i')
#
##########
# MLE:
# Initialisation of maximization algorithm:
rho_0<-0.5
a_0<-1
b_0<-3
# Compute the MLE
mleARAInf <- mle.vam(Time & Type ~ (ARAInf(rho_0) | Weibull(a_0,b_0)),data=simData)
(Est<-coef(mleARAInf))
# The corresponding plug in estimators can be plotted:
plot(mleARAInf,'i-cm')
# or
plot(mleARAInf,'I',col='darkblue',cm.col='red')
#
##########
#likelihood 3D plot
require(rgl)
rhos<-seq(0,1,0.1)
betas<-seq(0.1,6,0.1)
lnL<-c()
for (rho in rhos){
  for (beta in betas)
  {
    lnL<-c(lnL,contrast(mleARAInf,c(1,beta,rho)))
  }
}
lnL<-matrix(data=lnL,nrow=length(rhos),ncol=length(betas),byrow=TRUE)
persp3d(rhos, betas, lnL, col = 'skyblue',zlim=c(-100,-80))

grid3d(c("rho", "b", "lnL"))
spheres3d(Est[3],Est[2],contrast(mleARAInf,c(1,Est[2],Est[3])),r=0.4,alpha=0.5,color="red",add=TRUE)
indMax<-which(lnL==max(lnL), arr.ind = TRUE)
spheres3d(rhos[indMax[1]],betas[indMax[2]],contrast(mleARAInf,c(1,betas[indMax[2]],rhos[indMax[1]])),r=0.3,color="black",add=TRUE)
#
##########
# Monte Carlo study can be done:
nbSim<-1000
aEst<-rep(NaN,nbSim)
bEst<-rep(NaN,nbSim)
rhoEst<-rep(NaN,nbSim)
hasNotConverged<-rep(NaN,nbSim)
simDatas<-list()
for (i in 1:nbSim)
{
  simData<-simulate(simARAInf,n)
  simDatas<-c(simDatas,list(simData))
  mleARAInf <- mle.vam(Time & Type ~ (ARAInf(rho_0) | Weibull(a_0,b_0)),data=simData)
  run.mle.vam(mleARAInf,control=list(maxit=5000),verbose=FALSE)
  Est<-coef(mleARAInf)
  aEst[i]<-Est[1]
  bEst[i]<-Est[2]
  rhoEst[i]<-Est[3]
  hasNotConverged[i]<-mleARAInf$optim$convergence
}
#MLE rho
fivenum(rhoEst)
hist(rhoEst)
abline(v=0.4,lwd=4,col="red")
#
#####################
# multi trajectories
# only CM
#####################
#
# Simulation model:
# ARA 1 CM with \rho=0.4 
# Log linear h(t)=0.01*exp(0.3*t)
simARA1LL_Multi<-sim.vam(  ~ (ARA1(.4) | LogLinear(.01,0.3)))
#
# 1000 trajectories are simultaneously considered
n_Systems<-1000
# Each trajectory is composed of 10 MC
n_Events<-10
simData_Multi<-simulate(simARA1LL_Multi,n_Events,nb.system=n_Systems)
head(simData_Multi)
tail(simData_Multi)
#
# MLE
mleARA1LL_Multi <- mle.vam(System & Time & Type ~ (ARA1(.5) | LogLinear(1,1)),data=simData_Multi)
coef(mleARA1LL_Multi)
#
#####################
# MC et MP
#####################
#
#######
# Simulation model:
# Weibull h(t)=0.001*2.5*t^(1.5)
# CM : ARA1 (rhoMC=0.6)
# PM : ARA Infinit (rhoMP=0.4)
# PM at fixed intensity level
simCMPM<-sim.vam(  ~ (ARA1(.9) | Weibull(.001,2.5)) & (ARAInf(.4) | AtIntensity(0.2)))
#
# 1 trajectory
# with 50 events (CM+PM)
n<-50
#
simData=simulate(simCMPM,n)
head(simData)
tail(simData)
#also available AtVirtualAge and AtFailureProbability
#
# The corresponding intensity plot for the 20 first maintenance times
plot(simCMPM,'i',to=simData$Time[20])
#and for all the maintenance times the cumulative intensity
plot(simCMPM,'I',ylim=c(0,sum(simData$Type==-1)))
#
# MLE
mleCMPM <- mle.vam(Time & Type ~ (ARA1(.5) | Weibull(1,3)) & (ARAInf(0.5)),data=simData)
coef(mleCMPM)
plot(mleCMPM,'I-cm-pm',col='blue',add=TRUE)
#
############
# Simulation model:
# Weibull h(t)=0.001*2.5*t^(1.5)
# CM : ARA infinit (rhoMC=0.3)
# periodic PM randomly choosen upon two types
# PM : ARA Infinit (rhoMP=0.6)
# PM : ARA Infinit (rhoMP=-0.2)
simCMPM_Multi<-sim.vam(  ~ (ARAInf(.3) | Weibull(.001,2.5)) & (ARAInf(.6)+ARAInf(-.2) | Periodic(12,prob=c(0.6,0.4))))
#
# 50 trajectories are simultaneously considered
n_Systems<-50
# Each trajectory is composed of 20 maintenance
n_Events<-20
#
simData_Multi<-simulate(simCMPM_Multi,n_Events,nb.system=n_Systems)
head(simData_Multi)
tail(simData_Multi)
(prop<-table(simData_Multi$Type))
(prop[2][[1]])/(prop[2][[1]]+prop[3][[1]])
#
# The corresponding intensity plot for the second trajectory
plot(simCMPM_Multi,'i',system.index=2)
#
# MLE
mleCMPM_Multi <- mle.vam(System & Time & Type ~  (ARAInf(.5) | Weibull(1,3)) & (ARAInf(.5)+ARAInf(.5)),data=simData_Multi)
coef(mleCMPM_Multi)
plot(mleCMPM_Multi,'i-cm-pm',system.index=2,col='blue',add=TRUE)
#
##########
#Combined PM policiy
######
simCMPM<-sim.vam(  ~ (ABAO() | Weibull(.001,2.6)) & (ARAInf(.6)+ARAInf(-.2)+AGAN() | Periodic(12,prob=c(0.6,0.4))*AtIntensity(0.6)))
simData<-simulate(simCMPM,20)
plot(simCMPM,'i',col='darkblue')
#
###########
# PM policy NOT using the same model as simulation
###########
# The model used for PM plan 
modMPplan <- model.vam(Time & Type ~ (ARA1(.87) | Weibull(.0015,2.6)) & (ARAInf(.44)))
# The model used for failure simulation
simCMPM<-sim.vam( ~ (ARA1(.9) | Weibull(.002,2.5)) & (ARAInf(.25) | AtIntensity(0.35,modMPplan)) )
simData<-simulate(simCMPM,30)
update(modMPplan,simData)
plot(modMPplan,'i',ylim=c(0,0.5))
plot(simCMPM,'i',col='blue',add=TRUE)
abline(h=0.35)
#
########
# Different stoping methods for simulation
#######
simARAInf<-sim.vam(  ~ (ARAInf(.4) | Weibull(.001,2.5)))
(simulate(simARAInf,Time>25 ))
#
(simulate(simARAInf, Time>20 | Size>=5 ))
t_max<-100
s_max<-5
(simulate(simARAInf, Time>t_max | Size>=s_max ))
(simulate(simARAInf, Time>20 & Size>=5 ))
(simulate(simARAInf, T>100 & S>=5 ))
#
(simulate(simARAInf, T>(RightCensorship=40)) )
(simulate(simARAInf, T>(RightCensorship=Unif(20,40)) ,nb.system=3))
#
simCMPM<-sim.vam(  ~ (ARA1(.9) | Weibull(.001,2.5)) & (ARAInf(.4) | AtIntensity(0.2)))
simData<-simulate(simCMPM, Size[-1]>10 & Size[1]>=5 )
table(simData$Type)
simData<-simulate(simCMPM, Size[-1]>10 | Size[1]>=5 )
table(simData$Type)
simData<-simulate(simCMPM, (Size[-1]>=10 & Size[1]>=5) | Size>=15)
table(simData$Type)


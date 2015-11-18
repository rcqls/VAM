rm(list=ls())
require(VAM)

# Partie simulation
formSim <-formula( ~ (ARAInf(.4) | Weibull(.001,2.5)) & (ARA1(.7) + ARA1(.7) | Periodic(10,prob=c(.5,.5))))
simCpp <- sim.vam.cpp(formSim)
n=rep(30,1)
simDf=simulate(simCpp,n)

# Estimation un système
formMle <- update(formSim,Time & Type ~ .)
mleCpp <- mle.vam.cpp(formMle,data=simDf)
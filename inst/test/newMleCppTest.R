require(VAM)

if(!("testExp" %in% ls(1)))  testExp <- 1

## initialize all formulae 
formSim <- switch(testExp,
	~ (ARA1(.4) | Weibull(.001,2.5)),
	~ (ARAInf(.4) | Weibull(.001,2.5)),
	~ (ARA1(.4) | Weibull(.001,2.5)) & (ARA1(.7) + ARA1(.7) | Periodic(1,prob=c(.5,.5)) ),
	~ (ARAInf(.4) | Weibull(.001,2.5)) & (ARA1(.7) + ARA1(.7) | Periodic(10,prob=c(.5,.5)) ),
	~ (ARA1(.8) | LogLinear(exp(-5),0.5))
)


formModel <- update(formSim,Time & Type ~ .)
formMle <- formModel
formModelMulti <- update(formSim, System & Time & Type ~ .)
formMleMulti <- formModelMulti


## and initial parameters values!
par0 <- switch(testExp,
	c(1,2.5,.5),
	c(1,2.5,.5),
	c(1,2.5,.5,.5,.5),
	c(1,2.5,.5,.5,.5),
	c(0,1,.5)
)

## Create simulator
simCpp <- sim.vam.cpp(formSim)

## for single system:
nExp <- 1000
## simulate data for model and mle estimator  
simulate(simCpp,nExp) -> simDf
## create model to plot data with respect to any model
modelCpp <- model.vam.cpp( formModel ,data=simDf)
## create ML Estimator
mleCpp <- mle.vam.cpp( formMle ,data=simDf)

## for multi-systems:
nExp <- rep(1000,2)
## simulate data for model and mle estimator  
simulate(simCpp,nExp) -> simDf
## create model to plot data with respect to any model
modelMultiCpp <- model.vam.cpp( formModelMulti ,data=simDf)
## create ML Estimator
mleCppMulti <- mle.vam.cpp( formMleMulti ,data=simDf)
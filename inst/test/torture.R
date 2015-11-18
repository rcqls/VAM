require(VAM)
testExp <- 1

formSim <- switch(testExp,
	~ (ARA1(.4) | Weibull(.001,2.5)),
	~ (ARAInf(.4) | Weibull(.001,2.5)),
	~ (ARA1(.4) | Weibull(.001,2.5)) & (ARA1(.7) + ARA1(.7) | Periodic(1,prob=c(.5,.5)) ),
	~ (ARAInf(.4) | Weibull(.001,2.5)) & (ARA1(.7) + ARA1(.7) | Periodic(10,prob=c(.5,.5)) ),
	~ (ARA1(.8) | LogLinear(exp(-5),0.5))
)

formModel <- update(formSim,Time & Type ~ .)
formMle <- formModel
formModelMulti <- update(formSim,System & Time & Type ~ .)
formMleMulti <- formModelMulti

simCpp <- sim.vam.cpp(formSim)

nExp <- 1000
simulate(simCpp,nExp) -> simDf

for(i  in 1:2) {
	# cat("simCpp",i,"\n")
	# simCpp <- sim.vam.cpp(formSim)
	# gc()
	# cat("simulate\n")
	# simulate(simCpp,nExp) -> simDf
	# gc()
	cat("mleCpp\n")
	mleCpp <- mle.vam.cpp( formMle ,data=simDf)
	gc()
	cat("coef\n")
	print(coef(mleCpp,
		switch(testExp,
			c(1,2.5,.5),
	 		c(1,2.5,.5),
	 		c(1,2.5,.5,.5,.5),
	 		c(1,2.5,.5,.5,.5),
	 		c(0,1,.5),
	 	)->par0)->res
	)
	gc()

}



# modelCpp <- model.vam.cpp( formModel ,data=simDf)


# mleCpp <- mle.vam.cpp( formMle ,data=simDf)


# nExp <- rep(1000,2)
# simulate(simCpp,nExp) -> simDf
# mleCppMulti <- mle.vam.cpp( formMleMulti ,data=simDf)

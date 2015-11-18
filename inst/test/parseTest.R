require(VAM)
if(!("testExp" %in% ls(1))) testExp <- 1
formSim <- switch(testExp,
	~ ((ARA1(.4) | Weibull(.001,2.5))),
	~ (ARA1(.8) | LogLinear(-5.0,0.5)),
	~ (ARA1(.4) | Weibull(.001,2.5)) & (ARA1(.7) + ARA1(.7) | Periodic(1000,prob=c(.5,.5)) ),
	~ (ARAInf(.4) | Weibull(.001,2.5)) & (ARA1(.7) + ARA1(.7) | Periodic(10,prob=c(.5,.5)))
)

formMle <- update(formSim, Time & Type ~ .)

print(parse.vam.formula(NULL,formSim,Rcpp.mode=TRUE) ->modelSim)
print(parse.vam.formula(NULL,formMle,Rcpp.mode=TRUE)->modelMle)
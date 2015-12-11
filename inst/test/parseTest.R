require(VAM)
if(!("testExp" %in% ls(1))) testExp <- 1
formSim <- switch(testExp,
	~ ((ARA1(.4) | Weibull(.001,2.5))),
	~ (ARA1(.8) | LogLinear(-5.0,0.5)),
	~ (ARA1(.4) | Weibull(.001,2.5)) & (ARA1(.7) + ARA1(.7)),
	~ (ARA1(.4) | Weibull(.001,2.5)) & (ARA1(.7) + ARA1(.7) | Periodic(1000,prob=c(.5,.5)) ),
	~ (ARAInf(.4) | Weibull(.001,2.5)) & (ARA1(.7) + ARA1(.7) | Periodic(10,prob=c(.5,.5))),
	~ (ARAInf(.4) | Weibull(.001,2.5)) & (ARA1(.7) + ARA1(.6) + ARAInf(.6) | Periodic(10,prob=c(.5,.5)) * AtIntensity(1))
)

formMle <- update(formSim, Time & Type ~ .)

print(formSim)
print(parse.vam.formula(NULL,formSim) ->modelSim)
print(parse.vam.formula(NULL,formMle)->modelMle)

#print(parse.vam.formula(NULL,Time & Type ~ (ARAinf(1) | Weibull(.001,2.5))))
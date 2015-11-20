require(VAM)

selectExp <- function(idExp) {
	if(missing(idExp)) {
		cat("Please, execute 'selectExp(idExp)' with idExp chosen among:\n")
		lsExp()
	} else {
		.GlobalEnv$idExp <- idExp
		local({
			## initialize all formulae 
			formSimExps <- list(
				~ (ARA1(.4) | Weibull(.001,2.5)),
				~ (ARAInf(.4) | Weibull(.001,2.5)),
				~ (ARA1(.4) | Weibull(.001,2.5)) & (ARA1(.7) + ARA1(.7) | Periodic(1,prob=c(.5,.5)) ),
				~ (ARA1(.4) | Weibull(.001,2.5)) & (ARAInf(.7)| AtVirtualAge(10)),
				~ (ARA1(.4) | Weibull(.001,2.5)) & (ARAInf(.7)| AtIntensity(1.5)),
				~ (ARAInf(.4) | Weibull(.001,2.5)) & (ARA1(.7) + ARA1(.7) | Periodic(10,prob=c(.5,.5)) ),
				~ (ARA1(.8) | LogLinear(exp(-5),0.5))
			)
			formModelExps <- list(
				~ (ARA1(.4) | Weibull(.001,2.5)),
				~ (ARAInf(.4) | Weibull(.001,2.5)),
				~ (ARA1(.4) | Weibull(.001,2.5)) & (ARA1(.7) + ARA1(.7)),
				~ (ARA1(.4) | Weibull(.001,2.5)) & (ARAInf(.7)),
				~ (ARA1(.4) | Weibull(.001,2.5)) & (ARAInf(.7)),
				~ (ARAInf(.4) | Weibull(.001,2.5)) & (ARA1(.7) + ARA1(.7) ),
				~ (ARA1(.8) | LogLinear(exp(-5),0.5))
			)
			## and initial parameters values!
			par0Exps <- list(
				c(1,2.5,.5),
				c(1,2.5,.5),
				c(1,2.5,.5,.5,.5),
				c(1,2.5,.5,.5),
				c(1,2.5,.5,.5),
				c(1,2.5,.5,.5,.5),
				c(0,1,.5)
			)

			## set the infos related to idExp 
			formSimExp <- formSimExps[[idExp]]
			formModelExp <-  formModelExps[[idExp]]
			formModelMultiExp <- update(formModelExp, System & Time & Type ~ .)
			formModelExp <- update(formModelExp,Time & Type ~ .)
			formMleExp <- formModelExp
			formMleMultiExp <- formModelMultiExp
			par0Exp <- par0Exps[[idExp]]

			## Create simulator
			simCppExp <- sim.vam(formSimExp)

			## for single system:
			nExp <- 100
			## simulate data for model and mle estimator  
			simulate(simCppExp,nExp) -> simDfExp
			## create model to plot data with respect to any model
			modelCppExp <- model.vam( formModelExp ,data=simDfExp)
			## create ML Estimator
			mleCppExp <- mle.vam( formMleExp ,data=simDfExp)

			## for multi-systems:
			nExp <- rep(10,10)
			## simulate data for model and mle estimator  
			simulate(simCppExp,nExp) -> simDfMultiExp
			## create model to plot data with respect to any model
			modelCppMultiExp <- model.vam( formModelMultiExp ,data=simDfMultiExp)
			## create ML Estimator
			mleCppMultiExp <- mle.vam( formMleMultiExp ,data=simDfMultiExp)
		
			## set single system as the default
			nExp <- 100
		},.GlobalEnv)
		showExp()
	}
}

lsExp <-function() {
	for(i in seq_along(.GlobalEnv$formSimExps)) {
		cat(i,") ",deparse(.GlobalEnv$formSimExps[[i]]),"\n",sep="")
	}
}

showExp <- function() {
	cat("Formulae for selected experience (idExp=",.GlobalEnv$idExp,"):\n")
	cat("-> simCppExp: ",deparse(formSimExp),"\n")
	cat("-> modelCppExp: ",deparse(formModelExp),"\n")
	cat("-> modelCppMultiExp: ",deparse(formModelMultiExp),"\n")
	cat("-> mleCppExp: ",deparse(formMleExp),"\n")
	cat("-> mleCppMultiExp: ",deparse(formMleMultiExp),"\n")
}

setNExp <-function(nExp) {
	.GlobalEnv$nExp <- nExp
}

helpExp <- function() {
	cat(
		"The current workflow of this demo is:",
		"-> lsExp() to show all the experiments",
		"-> selectExp(idExp) to select the current experiment",
		"-> showExp() to show the current experiment",
		"-> set nExp (i.e. the number of simulated failures).",
		"   Notice that you are in the case of single or multi system(s) depending on the value of length of nExp.",
		"-> runMleExp() to provide estimates on newly simulated data",
		"-> plotExp(type) to plot model against newly simulated data (with type='v','h' or 'H')",
		sep="\n"
	)
}

runMleExp <- function(multi=length(nExp)>1) {
	if(!multi) {
		cat("Simulating...\n")
		simulate(simCppExp,nExp) -> simDfExp
		cat("Number of events:",nExp,"\n")
		update(mleCppExp,data=simDfExp)
		print(coef(mleCppExp,par0Exp))
	} else {
		cat("Simulating...\n")
		simulate(simCppExp,nExp,as.list=TRUE) -> simDfMultiExp
		cat("Table of number of system:\n")
		print(table(nExp))
		update(mleCppMultiExp,data=simDfMultiExp)
		print(coef(mleCppMultiExp,par0Exp))
	}
}

plotExp <-function(type="v",multi=length(nExp)>1,...) {
	if(!multi) {
			cat("Simulating...\n")
			simulate(simCppExp,nExp) -> simDfExp
			cat("Number of system:",nExp,"\n")
			update(modelCppExp,data=simDfExp)
			plot(modelCppExp,type=type,...)
	}
}

selectExp(1) #default

helpExp()


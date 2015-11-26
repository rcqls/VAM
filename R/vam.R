# Simulation: sim.vam or vam.sim or vam.gen??? 

 
sim.vam <- function(formula,data) {
	 
	self <- newEnv(sim.vam,formula=formula)

	PersistentRcppObject(self,new = {
		model <- parse.vam.formula(NULL,self$formula)
		rcpp <- new(SimVam,model)
		rcpp 
	})

	self
	 
}

simulate.sim.vam <- function(self, n=10, stop.time = Inf,as.list=FALSE) {
	rcpp <- self$rcpp()
	if(length(n)>1) {
		# multisystem
		if(as.list) df<-list()
		for(i in seq_along(n)) {
			df2 <- rcpp$simulate(n[i])[-1,]
			if(as.list) {
				df[[i]] <- df2 #rbind(data.frame(Time=0,Type=1),df2)
			} else {
				df2$System <- i
				df2<-df2[c(3,1:2)]
				df <- if(i==1) df2 else rbind(df,df2)
			}
		}
	} else df <- rcpp$simulate(n)[-1,]
	if(!as.list) rownames(df) <- 1:nrow(df)
	df
}

# Model part

model.vam <- function(formula,data) {
	self <- newEnv(model.vam,formula=formula,data=data)

	PersistentRcppObject(self,new = {
		model <- parse.vam.formula(NULL,self$formula)
		response <- model$response
		data <- data.frame.to.list.multi.vam(self$data,response)
		rcpp <- new(ModelVam,model,data)
		rcpp 
	})

	self
}

data.frame.to.list.multi.vam <- function(data,response) {
	# return data if it is already only a list!
	if(is.list(data) && !is.data.frame(data)) return(lapply(data,function(df) rbind(data.frame(Time=0,Type=1),df)))
	# otherwise
	if(length(response)==2) {
		if(length(intersect(response,names(data))) != 2) stop(paste0("Bad response:",response))
		tmp <- data[[response[1]]]
		data2 <- list(data.frame(Time=c(0,tmp[order(tmp)]),Type=c(1,data[[response[2]]][order(tmp)])))
	} else {
		if(length(intersect(response,names(data))) != 3) stop(paste0("Bad response:",response))
		syst0 <- unique(syst<-data[[response[1]]])
		data2 <- list()
		for(i in seq_along(syst0)) {
			df <- data[syst==syst0[i],response]
			tmp <- df[[response[2]]]
			data2[[i]] <- data.frame(Time=c(0,tmp[order(tmp)]),Type=c(1,df[[response[3]]][order(tmp)]))
		}
	}
	data2
}

# TODO: check data
check.data.vam <-function(data,response) { 
	if(all(data[[response[[-length(response)]]]])) {

	}
}

update.model.vam <- function(self,data) {
	if(!missing(data)) {
		model <- parse.vam.formula(NULL,self$formula)
		response <- model$response
		self$data <- data
		data2 <- data.frame.to.list.multi.vam(self$data,response)
		self$rcpp()$set_data(data2)
	}
}

mle.vam <- function(formula,data) {
	self <- newEnv(mle.vam,formula=formula,data=data)

	PersistentRcppObject(self,new = {
		model <- parse.vam.formula(NULL,self$formula)
		response <- model$response
		data <- data.frame.to.list.multi.vam(self$data,response)
		rcpp <- new(MLEVam,model,data)
		rcpp 
	})

	self
}

# to convert in Rcpp



params.model.vam <- params.sim.vam <- params.mle.vam <- function(self,param) {
	if(missing(param)) {
		 self$rcpp()$get_params()
	} else {
		self$rcpp()$set_params(param)
	}
}

update.mle.vam <- function(self,data) {
	if(!missing(data)) {
		model <- parse.vam.formula(NULL,self$formula)
		response <- model$response
		self$data <- data
		data2 <- data.frame.to.list.multi.vam(self$data,response)
		self$rcpp()$set_data(data2)
	}
}

# alpha is not considered in the estimation!
run.mle.vam <-function(obj,par0,fixed,method=NULL,verbose=TRUE,...) {
	rcpp <- obj$rcpp()
	## save the initial param
	if(is.null(obj$par0)) obj$par0 <- params(obj)
	## parameters stuff!
	if(missing(par0))  {
		if("par" %in% names(obj)) param <- obj$par #not the first run 
		else param<-params(obj)[-1] #first run
	} else if(is.null(par0)) param<-obj$par0[-1] else param<-par0[-1]
	## fixed and functions stuff!
	if(missing(fixed)) fixed<-rep(FALSE,length(param))
	else if(is.numeric(fixed)) {
		fixedInd<-fixed
		fixed<-rep(FALSE,length(param))
		fixed[fixedInd]<-TRUE
	}

	fn<-function(par) {
		##cat("param->");print(par);print(param[!fixed])
		param[!fixed]<-par
		#cat("param->");print(param)
		## All the commented part allows us to save the param when value is NaN
		#res<- 
		   -rcpp$contrast(c(1,param))
		# if(is.nan(res)) {
		# 	mode_param<-"contrast"
		# 	data_param<- rcpp$get_data(0)
		# 	data_param2 <- obj$data
		# 	save(param,mode_param,data_param,data_param2,file="/Users/remy/tmp/VAM/res.RData")
		# }
		# res
	}

 
	gr <- function(par) {
	    param[!fixed]<-par
	    #cat("param2->");print(param)
	    ## All the commented part allows us to save the param when value is NaN
		#res <-
	    -rcpp$gradient(c(1,param))[!fixed]
		# if(any(is.nan(res))) {
		# 	mode_param<-"gradient"
		# 	save(param,mode_param,file="/Users/remy/tmp/VAM/res.RData")
		# }
		# res
	}
  
  ## optim stuff!
  if(is.null(method) || method=="fast") {
    if(length(param[!fixed])>1) param[!fixed]<-(res <- optim(param[!fixed],fn,gr,method="Ne",...))$par
    if(is.null(method)) res<-optim(param[!fixed],fn,gr,method="CG",...)
  } else {
  	res<-optim(param[!fixed],fn,gr,method=method,...)
  }

  #fixed tips
  param[!fixed]<-res$par
  res$par<-param
  
  if(verbose) print(res)

  ## save stuff
  obj$fixed <- fixed
  obj$optim<-res
  obj$par<-res$par
  res$par
}

## Rmk: run.mle.vam is supposed to run many times to get the best estimate!
## Here, par=NULL forces initialisation update but does not ensure that it is the best estimate.
## TODO: try to find a best strategy or many strategies...
coef.mle.vam <- function(obj,par=NULL,method=NULL) {
	res <-run.mle.vam(obj,par,verbose=FALSE,method=method)
	if(obj$optim$convergence>0) cat("convergence=",obj$optim$convergence,"\n",sep="")
	alpha <- obj$rcpp()$alpha_est(c(1,res))
	res <- c(alpha,res)
	params(obj,res)
	res
}

# for both sim and mle

parse.vam.formula <- function(obj,formula) {
	if(formula[[1]] != as.name("~")) stop("Argument has to be a formula")
	if(length(formula) == 2) {
		response <- NULL
		cm <- formula[[2]]
	} else {
		tmp <- formula[[2]]
		## simplify parenthesis
		while(tmp[[1]] == as.name("(")) tmp <- tmp[[2]]
		if(tmp[[1]] != as.name("&") && length(tmp) != 3) stop("Left part of formula of the form 'Time & Type'!")
		if(length(tmp[[2]])==3 && tmp[[2]][[1]]==as.name("&")) {
			response <- c(as.character(tmp[[2]][[2]]),as.character(tmp[[2]][[3]]),as.character(tmp[[3]]))
		} else response <- c(as.character(tmp[[2]]),as.character(tmp[[3]]))
		cm <- formula[[3]]
	}
	## simplify parenthesis
	while(cm[[1]] == as.name("(")) cm <- cm[[2]]
	pms <- list()
	policy <- NULL
	if(there.is.pm <- (cm[[1]] == as.name("&"))) { # there is a PM part
		pm <- cm[[3]]
		cm <- cm[[2]]
		# deal with PM part
		if(pm[[1]] == as.name("(")) {
			pm <- pm[[2]]
			if(pm[[1]] != as.name("|")) {
				#stop("Need a policy to manage Preventive Maintenance")
				policy <- NULL
			} else {
				policy <- pm[[3]]
				if(is.name(policy[[1]])) {
					policy[[1]] <- as.name(paste0(as.character(policy[[1]]),".maintenance.policy"))
				}
				# TODO: add obj as argument of policy when needed
		
				# PMs
				pm <- pm[[2]]
			}
			# parser for pm
			parse.pm <- function(pm) {
				if(is.name(pm[[1]])) {
					pm[[1]] <- as.name(paste0(as.character(pm[[1]]),".va.model"))
				}
				##TO REMOVE (obj deleted): pm[[3]] <- as.name("obj")
				pm
			}
			cpt.pms <- 0
			while(pm[[1]] == as.name("+") ) {
				if(length(pm) == 3) {
					pms[[cpt.pms <- cpt.pms + 1]] <- parse.pm(pm[[3]])
					pm <- pm[[2]]
				}
			}
			pms[[cpt.pms <- cpt.pms + 1]] <- parse.pm(pm)
		} else stop("Need a policy to manage Preventive Maintenance")
	}
	# deal with CM PART
	cms <- list()
	 
	# parser for cm
	parse.cm <- function(cm) {
		# print(there.is.pm)
		# print(cm)
		if(there.is.pm) {
			if(cm[[1]] == as.name("(")) cm <- cm[[2]]
			else stop("CM needs a family!")
		}
		if(cm[[1]] != as.name("|")) stop("CM needs a family!")
		family <- cm[[3]]
		if(is.name(family[[1]])) {
			family[[1]] <- as.name(paste0(as.character(family[[1]]),".family.cm"))
		}
		cm <- cm[[2]]
		if(is.name(cm[[1]])) {
			cm[[1]] <- as.name(paste0(as.character(cm[[1]]),".va.model"))
		}
		##TO REMOVE (obj deleted): cm[[length(cm)+1]] <- as.name("obj")
		list(model=cm,family=family)
	}
	cpt.cms <- 0
	while( cm[[1]] == as.name("+") ) {
		if(length(cm) == 3) {
			cms[[cpt.cms <- cpt.cms + 1]] <- parse.cm(cm[[3]])
			cm <- cm[[2]]
		}
	}
	cms[[cpt.cms <- cpt.cms + 1]] <- parse.cm(cm)

	convert.cm <- function(cm) {

		list(
			model=list(
				name=as.character(cm$model[[1]]),
				##TO REMOVE (obj deleted): params=if(length(cm$model)==2) numeric(0) else sapply(cm$model[2:(length(cm$model)-1)],function(e) as.vector(eval(e)))
				params=as.vector(if(length(cm$model)==1) numeric(0) else sapply(cm$model[2:length(cm$model)],function(e) as.vector(eval(e))))
			),
			family=list(
				name=as.character(cm$family[[1]]),
				params=sapply(cm$family[-1],function(e) as.vector(eval(e)))
				## instead of : params=sapply(cm$family[-1],as.vector)
				## which does not work with negative real since element of tmp[-1] interpreted as call!
			)
		)
	}
	convert.pm <- function(pm) {
		list(
			name=as.character(pm[[1]]),
			##TO REMOVE (obj deleted): params=if(length(pm)==2) numeric(0) else sapply(pm[2:(length(pm)-1)],function(e) as.vector(eval(e)))
			params=as.vector(if(length(pm)==1) numeric(0) else sapply(pm[2:length(pm)],function(e) as.vector(eval(e))))
		)

	}
	convert.mp <- function(mp) {#maintenance policy
		if(is.null(mp)) list(name="None") 
		else {
			
			## The function defining the maintenance policy 
			## (registered in maintenance-policy-register.R or in any other R file)
			mp.fct <- eval(mp[[1]])
			## params used in the call mp
			pars <- as.list(match.call(mp.fct,mp))[-1]

			## Default values are then automatically completed using declaration of maintenance policy
			pars.default <- (as.list(mp.fct)->tmp)[-length(tmp)]
			pars.default <- pars.default[sapply(pars.default,function(e) nchar(as.character(e)))!=0]
			for(e in names(pars.default)) if(is.null(pars[[e]])) pars[[e]] <- pars.default[[e]]
			
			##print(list(pars=pars))

			list(
				name=as.character(mp[[1]]),
				params=lapply(pars,eval)
			)
		}
	}

	cms <- convert.cm(cms[[1]])
	
	list(
		response=response,
		models=c(list(cms$model),lapply(pms[rev(seq(pms))],convert.pm)),
		family=cms$family,
		pm.policy=convert.mp(policy)
	)
	
}




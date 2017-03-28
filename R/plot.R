## Provide cm.type or pm.type  with value "n" to not have cm and pm elements in the plot.
## cm.type or pm.type to NA means default value depending on type value.
plot.model.vam <- function(obj,type=c("v","virtual.age","i","intensity","I","cumulative","F","conditional.cdf","S","conditional.survival","f","conditional.pdf","d","data"),from,to,length.out=101,by,system.index=1,cm.type=NA,pm.type=NA,add=FALSE,preplot,...) {
	## preplot is used for plot.bayesian.vam where preplotting infos are precalculated
	if(missing(preplot)) {
		rcpp <- rcpp(obj)
		## IMPORTANT: sim.vam is now
		# d <- if(inherits(obj,"sim.vam")) rcpp$get_data() else
		d <- rcpp$get_data(system.index-1) #0 since one-system first!
		if(nrow(d)==0) stop("plot failed since data are required!")
		#print("d");print(d)

		mask <- TRUE
		if(missing(from)) from <- min(d$Time) else mask <- mask & d$Time>= from
		if(missing(to)) to <- max(d$Time) else mask <- mask & d$Time <= to

		if(missing(by)) by <- (to-from)/(length.out-1)
		infos <- rcpp$get_virtual_age_infos(by,from,to)
		infos <- infos[sapply(infos,function(e) !is.null(e))]
		#print("infos");print(infos);infos2 <<- infos
	} else {
		d <- preplot$d
		mask <- preplot$mask
		from <- preplot$from
		to <- preplot$to
		by <- preplot$by
		infos <- preplot$infos
	}

	## type
	if(length(grep("-",type))) { # deal with modifier -cm and -pm
		type <- gsub(" ","",type) # clear the white space
		if(length(grep("-cm",type))) {
			cm.type <- "n"
			type <- gsub("-cm","",type)
		}
		if(length(grep("-pm",type))) {
			pm.type <- "n"
			type <- gsub("-pm","",type)
		}
	}
	type <- match.arg(type)

	switch(type,
		virtual.age=,v={
			var <- "v"
			ylab<-"virtual age"
			if(!is.na(cm.type) && cm.type=="s") warning(paste0("cm.type argument could not be 's' for type '",type,"'!"))
			if(is.na(cm.type) || cm.type == "s") cm.type <- "p"
		},
		intensity=,i={
			var <- "i"
			ylab<-"intensity"
			if(!is.na(cm.type) && cm.type=="s") warning(paste0("cm.type argument could not be 's' for type '",type,"'!"))
			if(is.na(cm.type) || cm.type == "s") cm.type <- "p"
		},
		cumulative=,I={
			var <- "I"
			ylab<-"cumulative intensity"
			if(is.na(cm.type)) cm.type <- "s"
		},
		d=,data={
			var <- "d"
			ylab<-"CM counting process"
			if(is.na(cm.type)) cm.type <- "s"
		},
		conditional.cdf=,F={
			var <- "F"
			ylab<-"conditional cdf"
			if(!is.na(cm.type) && cm.type=="s") warning(paste0("cm.type argument could not be 's' for type '",type,"'!"))
			if(is.na(cm.type) || cm.type == "s") cm.type <- "p"
		},
		conditional.survival=,S={
			var <- "S"
			ylab<-"conditional survival function"
			if(!is.na(cm.type) && cm.type=="s") warning(paste0("cm.type argument could not be 's' for type '",type,"'!"))
			if(is.na(cm.type) || cm.type == "s") cm.type <- "p"
		},
		conditional.pdf=,f={
			var <- "f"
			ylab<-"conditional pdf"
			if(!is.na(cm.type) && cm.type=="s") warning(paste0("cm.type argument could not be 's' for type '",type,"'!"))
			if(is.na(cm.type) || cm.type == "s") cm.type <- "p"
		},
	)



	if(is.na(pm.type)) pm.type <- "l"

	## args
	args <- list(...)

	## args for cm
	args.cm <- NULL

	if(cm.type != "n") {
		switch(cm.type,
			p={
				args.cm <- args[names(args) %in% c('cm.pch','cm.col')]
				if(is.null(args.cm[['cm.pch']])) args.cm[['cm.pch']] <- "*"
			},
			l={
				args.cm <- args[names(args) %in% c('cm.lty','cm.lwd','cm.col')]
				if(is.null(args.cm[['cm.lty']])) args.cm[['cm.lty']] <- 2
				if(is.null(args.cm[['cm.lwd']])) args.cm[['cm.lwd']] <- 1
			},
			s={
				args.cm <- args[names(args) %in% c('cm.lty','cm.lwd','cm.col')]
				if(is.null(args.cm[['cm.lty']])) args.cm[['cm.lty']] <- 4
				if(is.null(args.cm[['cm.lwd']])) args.cm[['cm.lwd']] <- 1
			})
			if(is.null(args.cm[['cm.col']])) args.cm[['cm.col']] <- 1
	}
	##DEBUG: print(args.cm)

	## remove cm args
	args <- args[setdiff(names(args),names(args.cm))]

	## args for pm
	args.pm <- NULL
	pm.last.type <- max(d$Type) #Number of Type for PMs
	if(pm.last.type>0 && pm.type!="n") {
		switch(pm.type,
			p={
				args.pm <- args[names(args) %in% c('pm.pch','pm.col')]
				if(is.null(args.pm[['pm.pch']])) args.pm[['pm.pch']] <- 1:pm.last.type
				if(is.null(args.pm[['pm.col']])) args.pm[['pm.col']] <- (1:pm.last.type)+1
			},
			l={
				args.pm <- args[names(args) %in% c('pm.lty','pm.lwd','pm.col')]
				if(is.null(args.pm[['pm.lty']])) args.pm[['pm.lty']] <- rep(2,pm.last.type)
				if(is.null(args.pm[['pm.lwd']])) args.pm[['pm.lwd']] <- rep(1,pm.last.type)
				if(is.null(args.pm[['pm.col']])) args.pm[['pm.col']] <- (1:pm.last.type)+1
			})
			for(n in names(args.pm)) {
				if(length(args.pm[[n]]) != pm.last.type) {
					## recycling
					args.pm[[n]] <- rep(args.pm[[n]],pm.last.type/length(args.pm[[n]])+1)[1:pm.last.type]
				}
			}

	}
	##DEBUG: print(args.pm)

	## remove pm args
	args <- args[setdiff(names(args),names(args.pm))]

	## plot args
	is.args.plot <- names(args) %in% c('main','xlab','ylab','sub','asp','xlim','ylim')
	args.plot<- args[is.args.plot]
	args.lines <- args[!is.args.plot]
	if(is.null(args.plot[['xlab']])) args.plot[['xlab']] <- 'time'
	if(is.null(args.plot[['ylab']])) args.plot[['ylab']] <- ylab
	if(is.null(args.plot[['xlim']])) args.plot[['xlim']] <- c(from,to)

	ymax<-max(unlist(sapply(infos,function(e) e[[var]])))
	if(!add) do.call("plot",c(list(c(from,to),c(0,ymax),type="n"),args.plot))

	## 'v' or 'i' or 'I' plot
	if(var != "d") {
		t <- infos[[1]]$t
		v <- infos[[1]][[var]]
		for(i in seq_along(infos)[-1]) {
			t<-c(t,NA,infos[[i]]$t)
			v <- c(v,NA,infos[[i]][[var]])
		}
		do.call("lines",c(list(t,v),args.lines))
	}
	## cm

	if(cm.type != "n") {
		## IMPORTANT, first remove 'cm.' before calling plot method
		## not removed before because of collision with args.plot and args.lines!
		if(!is.null(args.cm)) names(args.cm) <- substring(names(args.cm),4)

		switch(cm.type,p={
			cm.call<-"points"
			args.cm[["x"]] <- d$Time[d$Type == -1 & mask]
			args.cm[["y"]] <- rep(0,sum(d$Type == -1 & mask))
		},l={
			cm.call<-"abline"
			args.cm[["v"]]<-d$Time[d$Type == -1 & mask]
		},s={
			cm.call<-"lines"
			args.cm[["x"]] <- c(0,d$Time[d$Type == -1 & mask],d$Time[nrow(d)])
			args.cm[["y"]] <- c(0,1:(sum(d$Type == -1 & mask)->tmp),tmp) + sum(d$Type == -1 & d$Time < from)
			args.cm[["type"]] <- "s"
		})
		##DEBUG: print(c(list(call=cm.call),args.cm))
		do.call(cm.call,args.cm)
	}

	## pm

	if(pm.type != "n") {
		## IMPORTANT, first remove 'pm.' before calling plot method
		ind <- d$Type>0 & d$Time>0 & mask
		if(!is.null(args.pm))  {
			names(args.pm) <- substring(names(args.pm),4)
			pm.types <- d$Type[ind]
			for(n in names(args.pm)) args.pm[[n]] <- args.pm[[n]][pm.types]
		}
		switch(pm.type,p={
			pm.call<-"points"

			args.pm[["x"]] <- d$Time[ind]
			args.pm[["y"]] <- rep(0,sum(ind))
		},l={
			pm.call<-"abline"
			args.pm[["v"]]<-d$Time[ind]
		})
		##DEBUG: print(c(list(call=pm.call),args.pm))
		do.call(pm.call,args.pm)
	}

}

plot.mle.vam  <- plot.sim.vam  <- plot.model.vam

preplots.bayesian.vam <- function(obj,from,to,length.out=101,by,system.index=1,type=c("97.5%","mean","2.5%"),filter=c("i","I"),nb.proposal=500) {
	rcpp <- rcpp(obj)
	## IMPORTANT: sim.vam is now
	# d <- if(inherits(obj,"sim.vam")) rcpp$get_data() else
	d <- rcpp$get_data(system.index-1) #0 since one-system first!
	if(nrow(d)==0) stop("plot failed since data are required!")
	#print("d");print(d)

	mask <- TRUE
	if(missing(from)) from <- min(d$Time) else mask <- mask & d$Time>= from
	if(missing(to)) to <- max(d$Time) else mask <- mask & d$Time <= to

	if(missing(by)) by <- (to-from)/(length.out-1)

	param <- if(obj$profile_alpha) obj$par0[-1] else obj$par0
	res <- NULL
	whichIndex<-1:obj$nb_proposal %in% sample(obj$nb_proposal,nb.proposal)
	firstIndex<-TRUE
	for(k in 1:obj$nb_proposal) {
		param[obj$par$ind[k]+!obj$profile_alpha] <- obj$par$estimate[k]
		if(whichIndex[k]){
			rcpp$set_params( if(obj$profile_alpha) c(obj$par$alpha[k],param) else param)
			infos <- rcpp$get_virtual_age_infos(by,from,to)
			infos <- infos[sapply(infos,function(e) !is.null(e))]
			if(firstIndex) {
				res <- lapply(infos,function(e) as.list(e[filter]))
				firstIndex<-FALSE
			}
			else {
				for(i in seq(res)) for(f in filter) {
					res[[i]][[f]] <- cbind(res[[i]][[f]],infos[[i]][[f]])
				}
			}
		}
	}
	preplots <- list()
	for(mode in type) {
		preplots[[mode]] <- list(d=d,mask=mask,from=from,to=to,by=by,infos=list())
		for(i in seq(res)) {
			preplots[[mode]]$infos[[i]] <- infos[[i]][c("t","v")]
			for(f in filter) {
				if(substr(mode,nchar(mode),nchar(mode))=="%") {
					alpha <- as.numeric(substr(mode,1,nchar(mode)-1))/100
					fct <- function(e) quantile(e,alpha)
				} else fct <- mean
				preplots[[mode]]$infos[[i]][[f]] <- apply(res[[i]][[f]],1,fct)
			}
		}
	}
	preplots
}

plot.bayesian.vam <- function(obj,type=c("i","intensity","I","cumulative","F","conditional.cdf","S","conditional.survival","f","conditional.pdf"),from,to,length.out=101,by,system.index=1,cm.type=NA,pm.type=NA,add=FALSE,nb.proposal=500,col=c("blue","black"),lty=c(1,3),lwd=c(1,1),...) {
	type <- match.arg(type)
	if(is.null(obj$par)){
		run(obj,nb=nb.proposal,history=TRUE)
	} else { if((!obj$history) || (obj$nb_proposal<nb.proposal)) run(obj,nb=nb.proposal,history=TRUE,profile.alpha=obj$profile_alpha,fixed=obj$fixed,sigma.proposal=obj$sigma_proposal)}
	if((!is.null(obj$preplots))&&((!missing(from) && from != obj$preplots[[1]]$from) || (!missing(to) && to != obj$preplots[[1]]$to) || (!missing(by) && by != obj$preplots[[1]]$by) || (!missing(length.out) && missing(by) && (obj$preplots[[1]]$by!= (obj$preplots[[1]]$to-obj$preplots[[1]]$from)/(length.out-1))))) obj$preplots <- NULL
	if(is.null(obj$preplots)||(obj$preplot_index!=system.index)||(obj$nb_proposal_preplot!=nb.proposal)) {
		obj$preplots <- preplots.bayesian.vam(obj,from,to,length.out,by,system.index,nb.proposal=nb.proposal)
		obj$preplot_index<-system.index
		obj$nb_proposal_preplot<-nb.proposal
	}
	## first one
	mode <- names(obj$preplots)[1]
	plot.model.vam(obj,type=type,cm.type=cm.type,pm.type=pm.type,preplot=obj$preplots[[mode]],col=(if(mode=="mean") col[1] else col[2]) ,lty=if(mode=="mean") lty[1] else lty[2],lwd=if(mode=="mean") lwd[1] else lwd[2],add=add,...)
	## the other plots
	for(mode in  names(obj$preplots)[-1]) {
		plot.model.vam(obj,type=type,cm.type=cm.type,pm.type=pm.type,preplot=obj$preplots[[mode]],add=TRUE,col=(if(mode=="mean") col[1] else col[2]) ,lty=if(mode=="mean") lty[1] else lty[2],lwd=if(mode=="mean") lwd[1] else lwd[2])
	}
}

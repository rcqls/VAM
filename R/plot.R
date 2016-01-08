## Provide cm.type or pm.type  with value "n" to not have cm and pm elements in the plot.
## cm.type or pm.type to NA means default value depending on type value.
plot.model.vam <- function(obj,type=c("v","virtual.age","i","intensity","I","cumulative"),from,to,by=0.1,system.index=0,cm.type=NA,pm.type=NA,add=FALSE,...) {
	rcpp <- rcpp(obj)
	d <- if(inherits(obj,"sim.vam")) rcpp$get_data() else rcpp$get_data(system.index) #0 since one-system first!
	infos <- rcpp$get_virtual_age_infos(by)
	if(missing(from)) from <- min(d$Time)
	if(missing(to)) to <- max(d$Time)

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
		}
	)



	if(is.na(pm.type)) pm.type <- "l"

	## args
	args <- list(...)

	## args for cm
	args.cm <- NULL
	if(!cm.type != "n") {
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

	}

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
	t <- infos[[1]]$t
	v <- infos[[1]][[var]]
	for(i in seq_along(infos)[-1]) {
		t<-c(t,NA,infos[[i]]$t)
		v <- c(v,NA,infos[[i]][[var]])
	}
	do.call("lines",c(list(t,v),args.lines))

	## cm

	if(cm.type != "n") {
		## IMPORTANT, first remove 'cm.' before calling plot method
		## not removed before because of collision with args.plot and args.lines!
		if(!is.null(args.cm)) names(args.cm) <- substring(names(args.cm),4)

		switch(cm.type,p={
			cm.call<-"points"
			args.cm[["x"]] <- d$Time[d$Type == -1]
			args.cm[["y"]] <- rep(0,sum(d$Type == -1))
		},l={
			cm.call<-"abline"
			args.cm[["v"]]<-d$Time[d$Type == -1]
		},s={
			cm.call<-"lines"
			args.cm[["x"]] <- c(0,d$Time[d$Type == -1],d$Time[nrow(d)])
			args.cm[["y"]] <- c(0,1:(sum(d$Type == -1)->tmp),tmp)
			args.cm[["type"]] <- "s"
		})
		do.call(cm.call,args.cm)
	}

	## pm

	if(pm.type != "n") {
		## IMPORTANT, first remove 'pm.' before calling plot method
		if(!is.null(args.pm))  names(args.pm) <- substring(names(args.pm),4)

		switch(pm.type,p={
			pm.call<-"points"
			args.pm[["x"]] <- d$Time[d$Type == -1]
			args.pm[["y"]] <- rep(0,sum(d$Type == -1))
		},l={
			pm.call<-"abline"
			args.pm[["v"]]<-d$Time[d$Type == -1]
		})
		do.call(pm.call,args.pm)
	}

}

plot.mle.vam  <- plot.sim.vam  <- plot.model.vam

plot.model.vam <- function(obj,type=c("virtual.age","intensity","cumulative","h","H"),from,to,by=0.1,...) {
	rcpp <- rcpp(obj)
	d <- if(inherits(obj,"sim.vam.cpp")) rcpp$get_data() else rcpp$get_data(0) #0 since one-system first!
	infos <- rcpp$get_virtual_age_infos(by)
	if(missing(from)) from <- min(d$Time)
	if(missing(to)) to <- max(d$Time)

	type <- match.arg(type)
	switch(type,
		virtual.age={
			var <- "v"
			ylab<-"virtual age"
		},
		intensity=,h={
			var <- "h"
			ylab<-"intensity"
		},
		cumulative=,H={
			var <- "H"
			ylab<-"cumulative"
		}
	)

	args <- list(...)
	is.args.plot <- names(args) %in% c('main','xlab','ylab','sub','asp','xlim','ylim')	
	args.plot<- args[is.args.plot]
	args.lines <- args[!is.args.plot]
	if(is.null(args.plot[['xlab']])) args.plot[['xlab']] <- 'time'
	if(is.null(args.plot[['ylab']])) args.plot[['ylab']] <- ylab
	if(is.null(args.plot[['xlim']])) args.plot[['xlim']] <- c(from,to)

	ymax<-max(unlist(sapply(infos,function(e) e[[var]])))
	do.call("plot",c(list(c(from,to),c(0,ymax),type="n"),args.plot))
	

	t <- infos[[1]]$t
	v <- infos[[1]][[var]]
	for(i in seq_along(infos)[-1]) {
		t<-c(t,NA,infos[[i]]$t)
		v <- c(v,NA,infos[[i]][[var]])
	}
	do.call("lines",c(list(t,v),args.lines))

}

plot.mle.vam  <- plot.sim.vam  <- plot.model.vam 
EndAt <-function(size,time,type=NULL,...) {# ... can contain cache.size
	if(missing(time)) {
		if(is.null(type)) {
			# Notice that it is the only case where cache.size is fixed by advance!
			obj<-list(name="AtSize.stop.policy",size=size,cache.size=size)
		} else {
			obj <- list(name="AtMSize.stop.policy",size=size,type=type,...)
		}
	} else {
		obj <- list(name="AtTime.stop.policy",time=time,...)
	}
	class(obj)<-c(obj$name,"stop.policy")
	obj
}

"|.stop.policy" <-function(e1,e2) {
	elts<-list()
	if(inherits(e1,"Or.stop.policy")) elts <- c(elts,e1$policies) else elts[[length(elts)+1]] <- e1
	if(inherits(e2,"Or.stop.policy")) elts <- c(elts,e2$policies) else elts[[length(elts)+1]] <- e2
	obj <- list(name="Or.stop.policy",policies=elts)
	class(obj)<-c(obj$name,"stop.policy")
	obj
}

"&.stop.policy" <-function(e1,e2) {
	elts<-list()
	if(inherits(e1,"And.stop.policy")) elts <- c(elts,e1$policies) else elts[[length(elts)+1]] <- e1
	if(inherits(e2,"And.stop.policy")) elts <- c(elts,e2$policies) else elts[[length(elts)+1]] <- e2
	obj <- list(name="And.stop.policy",policies=elts)
	class(obj)<-c(obj$name,"stop.policy")
	obj
}


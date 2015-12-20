EndAt <-function(size,time,...) {# ... can contain cache.size
	if(missing(time)) {
		# Notice that it is the only case where cache.size is fixed by advance!
		obj<-list(name="AtRun.stop.policy",nb=size,cache.size=size)
	} else {
		obj <- list(name="AtTime.stop.policy",time=time,...)
	}
	class(obj)<-c(obj$name,"stop.policy")
	obj
}
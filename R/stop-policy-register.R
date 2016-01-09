## TO REMOVE
EndAt <-function(size,time,type=NULL,...) {# ... can contain cache.size
	if(missing(time)) {
		if(is.null(type)) {
			# Notice that it is the only case where cache.size is fixed by advance!
			obj<-list(name="SizeGreaterThan.stop.policy",size=size,cache.size=size)
		} else {
			obj <- list(name="SizeOfTypeGreaterThan.stop.policy",size=size,type=type,...)
		}
	} else {
		obj <- list(name="TimeGreaterThanCensorship.stop.policy",time=time,...)
	}
	class(obj)<-c(obj$name,"stop.policy")
	obj
}

SizeGreaterThan <- function(size) {
	obj <- list(name="SizeGreaterThan.stop.policy",size=size,cache.size=size)
	class(obj)<-c(obj$name,"stop.policy")
	obj
}

SizeOfTypeGreaterThan <- function(size,type,...) {
	obj <- list(name="SizeOfTypeGreaterThan.stop.policy",size=size,type=type,...)
	class(obj)<-c(obj$name,"stop.policy")
	obj
}

AtCensorship <- TimeGreaterThanCensorship <- function(censorship,...) {# ... can contain cache.size
	if(inherits(censorship,"formula")) {
		obj <- list(name="TimeGreaterThanCensorship.stop.policy",time=-1,time.expr=list(expr=censorship[[2]],env=parent.frame()),...)
	} else {
		# expression("1+1")[[1]] is here a useless  arbitrary call
		obj <- list(name="TimeGreaterThanCensorship.stop.policy",time=censorship,time.expr=list(expr=expression(1+1)[[1]],env=parent.frame()),...)
	}
	class(obj)<-c(obj$name,"stop.policy")
	obj
}


TimeGreaterThan <- function(time,...) {# ... can contain cache.size
	obj <- list(name="TimeGreaterThan.stop.policy",time=time,...)
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

parse.stop.policy <- function(ch) {
  expr <- ch

  ## facility to locally update expr with rules
  Rule <- function(patt,repl) expr<<-gsub(patt,repl,expr,perl=TRUE)

  # Size of Type rules
  Rule("\\s","")
  Rule("S\\[(?!Type=)","S\\[Type=") #'S[' not followed by 'Type='
  Rule("S\\[Type=","Size\\[Type=")
	Rule("Size\\[(?!Type=)","Size\\[Type=")

  #Size rule
  Rule("S(?!(?:\\[Type=|ize))","Size")

  #Time rule
  Rule("T(?=\\s*>)","Time")

  #Right Censorship rule
  Rule("(?:\\()RC=","(RightCensorship=")

  parsed.expr <- parse(text=expr)[[1]]

  simplify.parenthesis <- function(e) {
    if(length(e)==1) return(e)
    if(length(e)==2 && e[[1]] == as.name("(") && length(e[[2]])==2 && e[[2]][[1]] == as.name("(")) e<-e[[2]]
    as.call(c(e[[1]],sapply(e[-1],simplify.parenthesis)))
  }

  parsed.expr <- simplify.parenthesis(parsed.expr)

	expand.distrib <- function(e) {
    if(is.name(e[[1]]) && (substring(e[[1]],1,1) %in% LETTERS)) {
      e[[1]] <- as.name(paste0("r",tolower(substring(e[[1]],1,1)),substring(e[[1]],2)))
      as.call(c(as.name("~"),as.call(c(e[[1]],1,as.list(e)[-1]))))
    } else e
  }

  expand <- function(e) {
    if(length(e)==1) return(e)

    ## Time > ???? or Time > (RightCensorship=???)
    if(e[[1]] == as.name(">") && e[[2]]==as.name("Time")) {
      ## Time > (RightCensorship=???)
      if(length(e[[3]])==2 && e[[3]][[1]]==as.name("(") && length(e[[3]][[2]])==3 && e[[3]][[2]][[1]]==as.name("=") && e[[3]][[2]][[2]]==as.name("RightCensorship")) {
        return(as.call(c(as.name("TimeGreaterThanCensorship"),censorship=expand.distrib(e[[3]][[2]][[3]]))))
      } else {## Time > ????
        return(as.call(c(as.name("TimeGreaterThan"),time=e[[3]])))
      }
    }

    ## Size[Type=???] > ????
    if(e[[1]] == as.name(">") && length(e[[2]])==3 && e[[2]][[1]]==as.name("[") && e[[2]][[2]]==as.name("Size") && names(e[[2]])[[3]]=="Type") {
      return(as.call(c(as.name("SizeOfTypeGreaterThan"),size=e[[3]],type=e[[2]][[3]])))
    }

    ## Size > ???
    if(e[[1]] == as.name(">") && e[[2]]==as.name("Size")) {
      return(as.call(c(as.name("SizeGreaterThan"),size=e[[3]])))
    }

    as.call(c(e[[1]],sapply(e[-1],expand)))
  }

  res <- expand(parsed.expr)
	attr(res,"user.entry") <- ch
	res
}

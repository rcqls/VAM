parse.stop.policy.devel <- function(ch) {
  expr <- ch

  ## facility to locally update e with rules
  Rule <- function(patt,repl) expr<<-gsub(patt,repl,expr,perl=TRUE)

  # Size of Type rules
  Rule("\\s","")
  Rule("S\\[(?!Type=)","S\\[Type=") #'S[' not followed by 'Type='
  Rule("S\\[Type=","Size\\[Type=")

  #Size rule
  Rule("S(?!(?:\\[Type=|ize))","Size")

  #Time rule
  Rule("T(?!(?:ime|ype))","Time")

  #Right Censorship rule
  Rule("(?:\\()RC=","(RightCensorship=")

  pe <- parse(text=expr)[[1]]

  simplify.parenthesis <- function(e) {
    if(length(e)==1) return(e)
    if(length(e)==2 && e[[1]] == as.name("(") && length(e[[2]])==2 && e[[2]][[1]] == as.name("(")) e<-e[[2]]
    as.call(c(e[[1]],sapply(e[-1],simplify.parenthesis)))
  }

  pe <- simplify.parenthesis(pe)

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
        print(e[[3]][[2]][[3]])
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



  expand(pe)


}

test.parse.stop.policy <- function(ch,devel=FALSE) {
  o <- deparse(substitute(ch))
  if(devel) {
    c(ch=o,parse.stop.policy.devel(o))
  } else {
    c(ch=o,parse.stop.policy(o))
  }
}
require(VAM)
print(test.parse.stop.policy((S[3] > 10 & ((S > 10))) & Size > 12 & Time > 100 & T > (RC=30))) #,devel=T))
print(test.parse.stop.policy((S[3] > 10 & ((S > 10))) & Size > 12 & Time > 100 & T > (RC=Unif(50,100)))) #,devel=T))

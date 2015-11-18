## This file is copied here but will be deleted as soon as the package CqlsRcppPersistentObject will be available at the CRAN.
newEnv <- function (...,class.as.character) {
    args.call <- as.list(match.call())[-1]
    names.call <- names(args.call)
    if (is.null(names.call)) 
        names.call <- rep("", length(args.call))
    class <- if(missing(class.as.character)) as.character(args.call[nchar(names.call) == 0]) else class.as.character
    names.call <- names.call[nchar(names.call) > 0]
    obj <- new.env()
    for (nm in names.call) assign(nm,eval.parent(args.call[[nm]]),envir=obj)
    class(obj) <- class
    obj
}
params <- function(obj,...) UseMethod("params")

run <- function(obj,...) UseMethod("run")

## Use logLik which is standad method
##logLikelihood <- function(obj,...) UseMethod("logLikelihood")

contrast <- function(obj,...) UseMethod("contrast")

## mean exits as a method but not sd. So introducing sigma
sigma <- function(obj,...) UseMethod("sigma")

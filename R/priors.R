## register priors but also

#this update methods help for substituting mle
mean.Beta.prior <- function(obj) {
  obj$params[1]/(obj$params[1]+obj$params[2])
}

mean.Unif.prior <- function(obj) {
  (obj$params[1]+obj$params[2])/2
}

mean.Gamma.prior <- function(obj) {
  obj$params[1]*obj$params[2]
}

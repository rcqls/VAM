## register priors but also

#this update methods help for substituting mle
update.Beta.prior <- function(obj) {
  obj$params[1]/(obj$params[1]+obj$params[2])
}

update.Unif.prior <- function(obj) {
  (obj$params[1]+obj$params[2])/2
}

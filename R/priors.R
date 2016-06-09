## register priors but also

#this update methods help for substituting mle
mean.Beta.prior <- function(obj) {
  obj$params[1]/(obj$params[1]+obj$params[2])
}

sigma.Beta.prior <- function(obj) {
  a<-obj$params[1];b<-obj$params[2]
  sqrt(a*b/(a+b)^2*(a+b+1))
}

mean.Unif.prior <- function(obj) {
  (obj$params[1]+obj$params[2])/2
}

sigma.Unif.prior <- function(obj) {
  a<-obj$params[1];b<-obj$params[2]
  sqrt((b-a)^2/12)
}

mean.Gamma.prior <- function(obj) {
  obj$params[1]*obj$params[2]
}

sigma.Gamma.prior <- function(obj) {
  a<-obj$params[1];b<-obj$params[2]
  sqrt(a*b^2)
}

mean.NonInform.prior <- function(obj) {
  obj$params[1]/(obj$params[1]+obj$params[2])
}

sigma.NonInform.prior <- function(obj) {
  a<-obj$params[1];b<-obj$params[2]
  sqrt(a*b/(a+b)^2*(a+b+1))
}

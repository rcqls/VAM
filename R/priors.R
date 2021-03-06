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

mean.Norm.prior <- function(obj) {
  obj$params[1]
}

sigma.Norm.prior <- function(obj) {
  obj$params[2]
}

mean.LNorm.prior <- function(obj) {
  exp(obj$params[1]+(obj$params[2])^2/2)
}

sigma.LNorm.prior <- function(obj) {
  sqrt((exp((obj$params[2])^2)-1)*exp(2*obj$params[1]+(obj$params[2])^2))
}

mean.NonInform.prior <- function(obj) {
  obj$params[1]
}

sigma.NonInform.prior <- function(obj) {
  obj$params[2]
}

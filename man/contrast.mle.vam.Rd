\name{contrast.mle.vam}

\alias{contrast.mle.vam}

\title{Compute the contrast corresponding to the maximum likelihood estimation method for a virtual age model}

\description{\code{contrast.mle.vam} computes the contrast corresponding to the maximum likelihood estimation method for a virtual age model with Corrective Maintenance (CM) and planned Preventive Maintenance (PM). 
The difference between the log-likelihood \code{\link{logLikelihood.mle.vam}} and the contrast is due to the parameter
\eqn{\alpha} which represents the scale parameter of the time to failure distribution of the new unmaintained system, \eqn{h(t)} (see \code{\link{sim.vam}} for more details). In fact, the value of \eqn{\alpha} that maximizes the log-likelihood has a closed-form solution, function of the others parameters values. The contrast corresponds to the value of the log-likelihood evaluated at this particular value of \eqn{\alpha}. Consequently, the contrast does not depend on \eqn{\alpha}, that is to say of \code{par0[1]}.}

\usage{
contrast.mle.vam(obj,par0,with_value=TRUE,with_gradient=FALSE,with_hessian=FALSE)
}

\arguments{
\item{\code{obj}}{an object of class \code{\link{mle.vam}} specifying the virtual age model considered and the data set of observations.}
\item{\code{par0}}{an optional argument specifying the parameter values at which the contrast is computed.
 \code{par0} is a vector whose length is equal to the total number of parameters in the formula used for defining the object \code{obj}. The successive values of \code{par0} respectively refer to the parameters of time to failure distribution of the new unmaintained system, of the CM effect model and finally of the PM effect models (if defined and in the same order as they appear in the \code{obj} formula).

If \code{par0} is not specified and the optimization algorithm has never been executed on \code{obj}(via \code{\link{coef.mle.vam}} or \code{\link{run.mle.vam}}), the contrast is computed on the parameter values specified in the \code{obj} formula.
Otherwise, the contrast is computed on the last parameters values obtained by the optimization algorithm for \code{obj}. }
\item{\code{with_value}}{a logical which indicates if the value of the contrast has to be computed.}
\item{\code{with_gradient}}{a logical which indicates if the gradient of the contrast has to be computed.}
\item{\code{with_hessian}}{a logical which indicates if the hessian of the contrast has to be computed.}
}

\value{\itemize{
\item If only \code{with_value} is \code{TRUE}, the method produces the contrast value.
\item If only \code{with_gradient} is \code{TRUE}, the method produces a vector corresponding to the gradient of the contrast,
\item If only \code{with_hessian} is \code{TRUE}, the method produces a matrix corresponding to the hessian of the contrast.
 }
 Otherwise, the method produces a list of the contrast characteristics for which the corresponding argument is \code{TRUE}.}

 \author{L. Doyen and R. Drouilhet}

\seealso{
\code{\link{run.mle.vam}} to compute the MLE.
\code{\link{coef.mle.vam}} to extract the parameters value of the MLE.
\code{\link{logLikelihood.mle.vam}} to compute the log-likelihood.
}

\examples{
simARAInf<-sim.vam(  ~ (ARAInf(.4) | Weibull(.001,2.5)))
simData<-simulate(simARAInf,30)
mleARAInf <- mle.vam(Time & Type ~ (ARAInf(0.5) | Weibull(1,3)),data=simData)
Est<-coef(mleARAInf)
contrast(mleARAInf)
contrast(mleARAInf,Est,c(TRUE,TRUE,TRUE))
logLikelihood(mleARAInf,Est,c(TRUE,TRUE,TRUE))

require(rgl)
rhos<-seq(0,1,0.1)
betas<-seq(0.1,6,0.1)
lnL<-c()
for (rho in rhos){
  for (beta in betas)
  {
    lnL<-c(lnL,contrast(mleARAInf,c(1,beta,rho)))
  }
}
lnL<-matrix(data=lnL,nrow=length(rhos),ncol=length(betas),byrow=TRUE)
persp3d(rhos, betas, lnL, col = 'skyblue',zlim=c(-100,max(lnL)))

grid3d(c("rho", "b", "lnL"))
spheres3d(Est[3],Est[2],contrast(mleARAInf,c(Est[1],Est[2],Est[3])),r=0.4,alpha=0.5,color="red",add=TRUE)
indMax<-which(lnL==max(lnL), arr.ind = TRUE)
spheres3d(rhos[indMax[1]],betas[indMax[2]],contrast(mleARAInf,c(1,betas[indMax[2]],rhos[indMax[1]])),r=0.3,color="black",add=TRUE)

}

\keyword{virtual age}
\keyword{imperfect maintenance}
\keyword{corrective maintenance}
\keyword{preventive maintenance}
\keyword{ARA}
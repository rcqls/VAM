\name{logLikelihood.mle.vam}

\alias{logLikelihood.mle.vam}

\title{Compute the log-likelihood for a virtual age model}

\description{\code{logLikelihood.mle.vam} computes the log-likelihood of the parameters of a virtual age model for Corrective Maintenance (CM) and planned Preventive Maintenance (PM).}

\usage{
logLikelihood.mle.vam(obj,par0,type=c(TRUE,FALSE,FALSE))
}

\arguments{
\item{\code{obj}}{an object of class \code{\link{mle.vam}} specifying the virtual age model considered and the data set of observations.}
\item{\code{par0}}{an optional argument specifying the parameter values at which the log-likelihood is computed.
 \code{par0} is a vector whose length is equal to the total number of parameters in the formula used for defining the object \code{obj}. The successive values of \code{par0} respectively refer to the parameters of time to failure distribution of the new unmaintained system, of the CM effect model and finally of the PM effect models (if defined and in the same order as they appear in the \code{obj} formula).

If \code{par0} is not specified and the optimization algorithm has never been executed on \code{obj}(via \code{\link{coef.mle.vam}} or \code{\link{run.mle.vam}}), the log-likelihood is computed on the parameter values specified in the \code{obj} formula.
Otherwise, the log-likelihood is computed on the last parameters values obtained by the optimization algorithm for \code{obj}. }
\item{\code{type}}{a vector of logical of length 3, which indicates respectively what to compute: the log-likelihood value, its gradient and/or  its hessian.}
}

\value{If \code{sum(type)==1} and, \itemize{
\item \code{type[1]==TRUE}, the method produces the log-likelihood value,
\item \code{type[2]==TRUE}, the method produces a vector corresponding to the gradient of the log-likelihood,
\item \code{type[3]==TRUE}, the method produces a matrix corresponding to the hessian of the log-likelihood.
 }
 Otherwise, the method produces a list of the log-likelihood characteristics according to the \code{type} argument, as described above.}

 \author{L. Doyen and R. Drouilhet}

\seealso{
\code{\link{run.mle.vam}} to compute the MLE.
\code{\link{coef.mle.vam}} to extract the parameters value of the MLE.
\code{\link{contrast.mle.vam}} to compute the contrast associated to the MLE.
}

\examples{
simARAInf<-sim.vam(  ~ (ARAInf(.4) | Weibull(.001,2.5)))
simData<-simulate(simARAInf,30)
mleARAInf <- mle.vam(Time & Type ~ (ARAInf(0.5) | Weibull(1,3)),data=simData)
logLikelihood(mleARAInf,c(0.02,2.4,0.4))
Est<-coef(mleARAInf)
contrast(mleARAInf)
logLikelihood(mleARAInf)
logLikelihood(mleARAInf,Est,c(TRUE,TRUE,TRUE))
}

\keyword{virtual age}
\keyword{imperfect maintenance}
\keyword{corrective maintenance}
\keyword{preventive maintenance}
\keyword{ARA}
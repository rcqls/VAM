\name{coef.mle.vam}

\alias{coef.mle.vam}
\alias{coef}

\title{Extract the maximum likelihood estimator of a virtual age model}

\description{\code{coef.mle.vam} extract the Maximum Likelihood Estimator (MLE) of the parameters of a virtual age model for Corrective Maintenance (CM) and planned Preventive Maintenance (PM). If the optimization method has never been applied it previously call the method \code{\link{run.mle.vam}}.}

\usage{
\method{coef}{mle.vam}(obj,par=NULL,method=NULL,verbose=FALSE)
}

\arguments{
\item{\code{obj}}{an object of class \code{\link{mle.vam}} specifying the virtual age model considered and the data set of observations.}
\item{\code{par}}{an optional argument specifying the initial parameter values for the optimization algorithm.
 \code{par} is a vector whose length is equal to the total number of parameters in the formula used for defining the object \code{obj}. The successive values of \code{par} respectively refer to the parameters of time to failure distribution of the new unmaintained system, of the CM effect model and finally of the PM effect models (if defined and in the same order as they appear in the \code{obj} formula).

If \code{par} is not specified and the optimization algorithm has never been executed on \code{obj} (via \code{\link{coef.mle.vam}} or \code{\link{run.mle.vam}}), the optimization algorithm used the parameter values specified in the \code{obj} formula as initialization point.
Otherwise, the last parameters values obtained by the optimization algorithm for \code{obj} are used as initialization point. Then, if one execution of \code{run.mle.vam} is not sufficient to make the optimization algorithm converged (\code{$convergence!=1}), the method can be re-executed until convergence. }
\item{\code{method}}{an optimization method of \code{\link{optim}} function used for optimization. \code{method} can also be equal to the string \code{"fast"} in order to let the function choose the method known as the faster one.}
\item{\code{verbose}}{if \code{TRUE}, the value returned by \code{\link{optim}} function is printed.}
}

\value{The function extract the vector of parameter values obtained by the maximization algorithm applied to the likelihood. The successive values respectively refer to the parameters of time to failure distribution of the new unmaintained system (except the scale parameter), of the CM effect model and finally of the PM effect models (if defined and in the same order as they appear in the \code{obj} formula).}

 \author{L. Doyen and R. Drouilhet}

\seealso{
\code{\link{run.mle.vam}} to compute the MLE.
\code{\link{contrast.mle.vam}} to compute the contrast associated to the MLE.
\code{\link{logLikelihood.mle.vam}} to compute the log-likelihood.
}

\examples{
simARAInf<-sim.vam(  ~ (ARAInf(.4) | Weibull(.001,2.5)))
simData<-simulate(simARAInf,30)
mleARAInf <- mle.vam(Time & Type ~ (ARAInf(0.5) | Weibull(1,3)),data=simData)
coef(mleARAInf)

mleARAInf2 <- mle.vam(Time & Type ~ (ARAInf(0.5) | Weibull(1,3)),data=simData)
run.mle.vam(mleARAInf2,fixed=c(TRUE,FALSE))
coef(mleARAInf2)
}

\keyword{virtual age}
\keyword{imperfect maintenance}
\keyword{corrective maintenance}
\keyword{preventive maintenance}
\keyword{ARA}
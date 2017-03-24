\name{run.mle.vam}

\alias{run.mle.vam}
\alias{run}

\title{Compute the maximum likelihood estimator of a virtual age model}

\description{\code{run.mle.vam} computes the Maximum Likelihood Estimator (MLE) of the parameters of a virtual age model for Corrective Maintenance (CM) and planned Preventive Maintenance (PM). It also updates the parameters values with those obtained by the optimization algorithm.
The optimization function used is \code{\link{optim}}. }

\usage{
run.mle.vam(obj, par0, fixed, method=NULL, verbose=TRUE,...)
}

\arguments{
\item{\code{obj}}{an object of class \code{\link{mle.vam}} specifying the virtual age model considered and the data set of observations.}
\item{\code{par0}}{an optional argument specifying the initial parameter values for the optimization algorithm.
 \code{par0} is a vector whose length is equal to the total number of parameters in the formula used for defining the object \code{obj}. The successive values of \code{par0} respectively refer to the parameters of time to failure distribution of the new unmaintained system, of the CM effect model and finally of the PM effect models (if defined and in the same order as they appear in the \code{obj} formula).

If \code{par0} is not specified and the optimization algorithm has never been executed on \code{obj} (via \code{\link{coef.mle.vam}} or \code{\link{run.mle.vam}}), the optimization algorithm used the parameter values specified in the \code{obj} formula as initialization point.
Otherwise, the last parameters values obtained by the optimization algorithm for \code{obj} are used as initialization point. Then, if one execution of \code{run.mle.vam} is not sufficient to make the optimization algorithm converged (\code{$convergence!=1}), the method can be re-executed until convergence. }
\item{\code{fixed}}{an optional argument specifying the parameters for which the value is fixed to initialization. \code{fixed} can be a vector of logical whose length is equal the total number of parameters in the formula used for defining the object \code{obj}. The successive values of \code{fixed} respectively refer to the parameters of time to failure distribution of the new unmaintained system, of the CM effect model and finally of the PM effect models (if defined and in the same order as they appear in the \code{obj} formula). Otherwise \code{fixed} can be equivalently a vector of integers correspoding to the indexes of the true values in the previous \code{fixed} version.}
\item{\code{method}}{an optimization method of \code{\link{optim}} function used for optimization. \code{method} can also be equal to the string \code{"fast"} in order to let the function choose the method known as the faster one.}
\item{\code{verbose}}{if \code{TRUE}, the value returned by \code{\link{optim}} function is printed.}
\item{\code{...}}{further arguments passed to the \code{\link{optim}} function.}
}

\value{The method returns the results of the \code{optim} function called inside the \code{run} function.}

 \author{L. Doyen and R. Drouilhet}

\seealso{
\code{\link{coef.mle.vam}} to extract the parameters value of the MLE.
\code{\link{contrast.mle.vam}} to compute the contrast associated to the MLE.
\code{\link{logLik.mle.vam}} to compute the log-likelihood.
}

\examples{
simARAInf<-sim.vam(  ~ (ARAInf(.4) | Weibull(.001,2.5)))
simData<-simulate(simARAInf,30)
mleARAInf <- mle.vam(Time & Type ~ (ARAInf(0.5) | Weibull(1,3)),data=simData)
run(mleARAInf)

mleARAInf2 <- mle.vam(Time & Type ~ (ARAInf(0.5) | Weibull(1,3)),data=simData)
run(mleARAInf2,fixed=c(FALSE,TRUE,FALSE))
}

\keyword{virtual age}
\keyword{imperfect maintenance}
\keyword{corrective maintenance}
\keyword{preventive maintenance}
\keyword{ARA}

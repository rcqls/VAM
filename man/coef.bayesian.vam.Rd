\name{coef.bayesian.vam}
\alias{coef.bayesian.vam}
\alias{coef}

\title{Compute a Bayesian point estimate of the parameters of a virtual age model}

\description{For \code{bayesian.vam} object, compute a Bayesian estimate of the parameters. This estimate corresponds to the mean of the posterior distribution sampled with \code{\link{run.bayesian.vam}}.}

\usage{
coef.bayesian.vam(obj,new.run=FALSE,...)
}

\arguments{
	\item{\code{obj}}{an object of class \code{\link{bayesian.vam}} specifying the virtual age model considered, the prior distribution of the parameters and the data set of observations.}
  	\item{\code{new.run}}{an optional argument specifying if \code{\link{run.bayesian.vam}} methods has to be relaunched even if it has already been executed.}
  	\item{\code{...}}{some supplementary arguments for the \code{\link{run.bayesian.vam}} method.}
}

\value{
 A numeric vector of the estimates corresponding respectively to the parameters of time to failure distribution of the new unmaintained system, of the CM effect model and finally of the PM effect models (if defined and in the same order as they appear in the \code{obj} formula).
 }

 \author{R. Drouilhet et L. Doyen}

\seealso{
\code{\link{bayesian.vam}} to define the Bayesian object.
\code{\link{run.bayesian.vam}} to compute the Bayesian method.
\code{\link{summary.bayesian.vam}} to produce a result summary of the Bayesian method.
\code{\link{hist.bayesian.vam}} for plotting the histogram of the posterior distribution of the parameters.
%\code{\link{formula.bayesian.vam}} to extract the original and estimated model.
\code{\link{plot.bayesian.vam}} for plotting estimating characteristics of the model.
}

\examples{
simARAInf<-sim.vam(  ~ (ARAInf(.4) | Weibull(.001,2.5)))
simData<-simulate(simARAInf,30)
bayesARAInf <- bayesian.vam(Time & Type ~ (ARAInf(~Unif(0,1)) | Weibull(~Unif(0,1),~Unif(2,4))),data=simData)
coef(bayesARAInf)
run(bayesARAInf,profile.alpha=TRUE)
coef(bayesARAInf)
coef(bayesARAInf,new.run=TRUE,par0=c(1e-2,2.5,0.5),fixed=2)
}

\keyword{virtual age}
\keyword{imperfect maintenance}
\keyword{corrective maintenance}
\keyword{preventive maintenance}
\keyword{ARA}
\keyword{Bayesian}
\name{summary.bayesian.vam}
\alias{summary.bayesian.vam}
\alias{summary}

\title{Summarizing Bayesian method for a virtual age model}

\description{For \code{bayesian.vam} object, produce a summary of the Bayesian estimation method corresponding to the posterior distribution sampled obtained with \code{\link{run.bayesian.vam}}.}

\usage{
summary.bayesian.vam(obj,alpha=0.05,new.run=FALSE,...)
}

\arguments{
	\item{\code{obj}}{an object of class \code{\link{bayesian.vam}} specifying the virtual age model considered, the prior distribution of the parameters and the data set of observations.}
	\item{\code{alpha}}{an optional argument specifying the level (\code{1-alpha}) of the credibility intervals for the parameters.}
  	\item{\code{new.run}}{an optional argument specifying if \code{\link{run.bayesian.vam}} methods has to be relaunched even if it has already been executed.}
  	\item{\code{...}}{some supplementary arguments for the \code{\link{run.bayesian.vam}} method.}
}

\value{
Print a summary of the Bayesian method \code{\link{run.bayesian.vam}} applied to the object \code{obj}:
\itemize{
	\item{The initialization values of the Gibbs algorithm.}
	\item{The Bayesian point estimates of the parameters, corresponding to the means of the sampled marginal posterior distributions.}
	\item{A (\code{1-alpha}) credibility interval for the parameters, corresponding to the \code{alpha/2} and \code{1-alpha/2} quantiles of the marginal posterior distribution.}
	\item{The number of accepted marginal simulated parameters values of the posterior distribution in the Metropolis Hasting step.}
	\item{The acceptation rate of the Metropolis Hasting step.}
}
 }

\details{The acceptation rate of the Metropolis Hasting step can be used to calibrate the standard deviation of the instrumental distribution (see the argument \code{sigma.proposal} of the \code{\link{run.bayesian.vam}} method).}

 \author{R. Drouilhet et L. Doyen}

\seealso{
\code{\link{bayesian.vam}} to define the Bayesian object.
\code{\link{run.bayesian.vam}} to compute the Bayesian method.
\code{\link{coef.bayesian.vam}} to extract the parameters estimation values of the Bayesian method.
\code{\link{hist.bayesian.vam}} for plotting the histogram of the posterior distribution of the parameters.
%\code{\link{formula.bayesian.vam}} to extract the original and estimated model.
\code{\link{plot.bayesian.vam}} for plotting estimating characteristics of the model.
}

\examples{
simARAInf<-sim.vam(  ~ (ARAInf(.4) | Weibull(.001,2.5)))
simData<-simulate(simARAInf,30)
bayesARAInf <- bayesian.vam(Time & Type ~ (ARAInf(~Unif(0,1)) | Weibull(~Unif(0,1),~Unif(2,4))),data=simData)
summary(bayesARAInf)
run(bayesARAInf,profile.alpha=TRUE)
summary(bayesARAInf)
summary(bayesARAInf,par0=c(1e-2,2.5,0.5),fixed=2)
}

\keyword{virtual age}
\keyword{imperfect maintenance}
\keyword{corrective maintenance}
\keyword{preventive maintenance}
\keyword{ARA}
\keyword{Bayesian}
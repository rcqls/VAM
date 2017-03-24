\name{hist.bayesian.vam}
\alias{hist.bayesian.vam}
\alias{hist}

\title{Plot the posterior distribution of the parameters for a virtual age model}

\description{For \code{bayesian.vam} object, plot the posterior distribution of the parameters. This posterior distribution corresponds to the one sampled with \code{\link{run.bayesian.vam}}.}

\usage{hist.bayesian.vam(obj,i=1,...)}

\arguments{
	\item{\code{obj}}{an object of class \code{\link{bayesian.vam}} specifying the virtual age model considered, the prior distribution of the parameters and the data set of observations.}
  	\item{i}{an optional argument specifying the index of the parameter for which the posterior distribution is plotted.}
  	\item{\code{...}}{some supplementary arguments for the generic \code{\link{hist}} function.}
}

 \author{R. Drouilhet et L. Doyen}

\seealso{
\code{\link{bayesian.vam}} to define the Bayesian object.
\code{\link{run.bayesian.vam}} to compute the Bayesian method.
\code{\link{summary.bayesian.vam}} to produce a result summary of the Bayesian method.
\code{\link{coef.bayesian.vam}} to extract the parameters estimation values of the Bayesian method.
%\code{\link{formula.bayesian.vam}} to extract the original and estimated model.
\code{\link{plot.bayesian.vam}} for plotting estimating characteristics of the model.
}

\examples{
simARAInf<-sim.vam(  ~ (ARAInf(.4) | Weibull(.001,2.5)))
simData<-simulate(simARAInf,30)
bayesARAInf <- bayesian.vam(Time & Type ~ (ARAInf(~Unif(0,1)) | Weibull(~Unif(0,1),~Unif(2,4))),data=simData)
run(bayesARAInf,profile.alpha=TRUE)
hist(bayesARAInf,3,main="Posterior distribution of rho")
}

\keyword{virtual age}
\keyword{imperfect maintenance}
\keyword{corrective maintenance}
\keyword{preventive maintenance}
\keyword{ARA}
\keyword{Bayesian}
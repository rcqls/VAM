\name{bayesian.vam}
\alias{bayesian.vam}

\title{Virtual age model for Bayesian estimation}

\description{\code{bayesian.vam} is used to define a virtual age model for Corrective Maintenance (CM) and planned Preventive Maintenance (PM). The object define with \code{bayesian.vam} can be used to compute Bayesian estimators of the parameters thanks to the \code{\link{run.bayesian.vam}} method.}

\usage{
bayesian.vam(formula, data)
}

\arguments{
  \item{\code{formula}}{a symbolic description of the virtual age model and observations. The details of formula specifications are given under `Details'.}
  \item{\code{data}}{a data frame or possibly a list (when several system are considered together) containing the observations.}
}

\details{
The symbolic description of the model done in \code{formula} has the form \code{response ~ model}.
\code{response} is a symbolic description of the data considered. The specifications are the same as those of \code{\link{model.vam}} function.
\code{model} is a symbolic description of the virtual age model considered. The specifications are globally similar those of \code{\link{sim.vam}} function. The difference with \code{\link{model.vam}} is that each parameter value of \code{model} has to be replaced by a symbolic description of the prior distribution of the parameter. This symbolic description has the form \code{~ Prior}. The available prior distribution for \code{Prior} are:
\itemize{
	\item{\code{Beta(a,b)}}{or \code{Be(a,b)} or \code{B(a,b)} for a beta distribution with parameters \code{a} and \code{b} (see \code{\link{dbeta}} for supplementary informations about this distribution, \code{a} and \code{b} represent respectively \code{shape} and \code{scale} arguments of \code{\link{dbeta}}),}
	\item{\code{Gamma(a,s)}}{or \code{G(a,s)} for a gamma distribution with parameters \code{a} and \code{s} (see \code{\link{dgamma}} for supplementary informations about this distribution),}
	\item{\code{Unif(a=0,b=1)}}{or \code{U(a=0,b=1)} for a uniform distribution on the interval [\code{a},\code{b}] (see \code{\link{dunif}} for supplementary informations about this distribution),}
	\item{\code{Norm(m=0,s=1)}}{or \code{N(m=0,s=1)} for a normal distribution with mean \code{m} and standard deviation \code{s} (see \code{\link{dnorm}} for supplementary informations about this distribution),}
		\item{\code{LNorm(m=0,s=1)}}{or \code{LN(m=0,s=1)} or \code{LogNorm(m=0,s=1)} for a log-normal distribution, that is to say Y=exp(X) follows a log-normal distribution if X follows a normal distribution with mean \code{m} and standard deviation \code{s},}
	\item{\code{NonInform(init=1,init_sigma=1)}}{or \code{NInf(init=1,init_sigma=1)} or \code{NI(init=1,init_sigma=1)} for a a standard Jeffrey's non informative prior distribution for Weibull parameters (with a prior density in x proportional to 1/x). \code{Init} represents the initialization of the maximum likelihood method used to compute the initialization of the Gibbs sampling algorithm in the \code{\link{run.bayesian.vam}} method (in order to evaluate the posterior distribution). \code{init_sigma} represents the standard deviation of instrumental distribution of the Metropolis Hasting step of Gibbs sampling algorithm used in \code{\link{run.bayesian.vam}} to evaluate the posterior distribution. \code{init_sigma} is used only if this standard deviation is not specified in \code{\link{run.bayesian.vam}}. Those characteristics are defined since this non informative prior distribution is not a proper distribution, then it does not have a mean and standard deviation.}
}
In addition, in this case the PM policy is useless, so it has not to be necessarily defined.
}

\value{
 The function produces an object of class \code{bayesian.vam} which contains the virtual age model considered, the prior distribution of the parameters and the corresponding observations.
 }

 \author{R. Drouilhet}

\seealso{
\code{\link{run.bayesian.vam}} to compute the Bayesian method.
\code{\link{coef.bayesian.vam}} to extract the parameters estimation values of the Bayesian method.
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

}

\keyword{virtual age}
\keyword{imperfect maintenance}
\keyword{corrective maintenance}
\keyword{preventive maintenance}
\keyword{ARA}
\keyword{Bayesian}

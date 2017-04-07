\name{run.bayesian.vam}
\alias{run.bayesian.vam}
\alias{run}

\title{Compute the Bayesian method for a virtual age model}

\description{For \code{bayesian.vam} object, the method \code{run.bayesian.vam} produces a sample from the posterior distribution of the parameters using a Gibbs sampling algorithm with a Metropolis Hasting step. The generated sample is memorized in the \code{bayesian.vam} object and can be used in particular by the \code{\link{coef.bayesian.vam}} and \code{\link{summary.bayesian.vam}} to produce Bayesian estimators of the parameters.}

\usage{
run.bayesian.vam(obj,par0,fixed,sigma.proposal,nb=100000,burn=10000,profile.alpha=FALSE,method=NULL,verbose=FALSE,history=FALSE,proposal='norm',...)

}

\arguments{
	\item{\code{obj}}{an object of class \code{\link{bayesian.vam}} specifying the virtual age model considered, the prior distribution of the parameters and the data set of observations.}
  \item{\code{par0}}{an optional argument specifying the initial parameter values for the Gibbs sampling method.
 	\code{par0} is a vector whose length is equal to the total number of parameters in the formula used for defining the object \code{obj}. The successive values of \code{par0} respectively refer to the parameters of time to failure distribution of the new unmaintained system, of the CM effect model and finally of the PM effect models (if defined and in the same order as they appear in the \code{obj} formula).

	If \code{par0} is not specified, the Gibbs sampling method is initialized at the maximum likelihood estimators parameters values (see \code{\link{mle.vam}}). In order to compute this maximum likelihood estimators, the optimization method is initialized at the means of the prior distributions.}
	\item{\code{fixed}}{an optional argument specifying the parameters for which the value is fixed to initialization. \code{fixed} can be a vector of logical whose length is equal the total number of parameters in the formula used for defining the object \code{obj}. The successive values of \code{fixed} respectively refer to the parameters of time to failure distribution of the new unmaintained system, of the CM effect model and finally of the PM effect models (if defined and in the same order as they appear in the \code{obj} formula). Otherwise \code{fixed} can be equivalently a vector of integers corresponding to the indexes of the true values in the previous \code{fixed} version.}
	\item{\code{sigma.proposal}}{an optional argument specifying the standard deviation of the instrumental distribution of the Metropolis Hasting step of the Gibbs sampling algorithm. That can be a seldom value and then this standard distribution will be the same for all the parameters. Otherwise, it can be vector whose length is equal to the total number of parameters.}
  	\item{\code{nb}}{an optional argument specifying when to stop the Gibbs sampling algorithm. If \code{history=FALSE}, \code{nb} represents the number of iterations of the Gibbs sampling algorithm. Otherwise, it represents the number of accepted simulated values for the different marginal posterior distributions of the parameters in the Metropolis Hasting step.}
  	\item{\code{burn}}{an optional argument specifying the number of burn in iterations of the Gibbs sampling algorithm.}
  	\item{\code{profile.alpha}}{an optional argument specifying if the likelihood is profiled in \eqn{\alpha} (the log likelihood is replaced by the contrast function, see \code{\link{contrast.mle.vam}}). If the likelihood is profiled the prior distribution of parameter \eqn{\alpha} is not used. However it must have been specified in the \code{\link{bayesian.vam}} object even if it is not used.}
  	\item{\code{method}}{an optional argument to be used by the maximum likelihood method for computing the initialization of the Gibbs algorithm (see \code{run.mle.vam}).}
  	\item{\code{verbose}}{an optional argument to be used by the maximum likelihood method for computing the initialization of the Gibbs algorithm (see \code{run.mle.vam}).}
  	\item{\code{history}}{an optional argument defining how to stop the Gibbs sampling method and how to memorized its result. The plotting method \code{\link{plot.bayesian.vam}} needs to run the Bayesian method with \code{history=TRUE}, but that takes more computing time.}
  	\item{\code{proposal}}{an optional argument specifying the instrumental distribution of the Metropolis Hasting step of the Gibbs sampling algorithm. Possible distributions are \code{'norm'} for normal distribution and \code{'lnorm'} for log-normal distribution. That can be a seldom value and then the used distribution will be the same for all the parameters. Otherwise, it can be vector whose length is equal to the total number of parameters.}
  	\item{\code{...}}{some supplementary arguments for the maximum likelihood method used to compute the initialization of the Gibbs algorithm (see \code{run.mle.vam}).}
}

\details{
The standard deviation of the instrumental distribution of the Metropolis Hasting play an important part. It can be specified with the \code{sigma.proposal} argument. If it is too big or too small regarding to the parameter value space, the sampling method can be very long to converge. In fact the standard deviations, \code{sigma.proposal}, specifies how the marginal parameters spaces are explored. If the standard deviation is too small, only a little part of the parameter space is explored. If it is too big, the simulated parameters values would be ever rejected in the Metropolis Hasting step.
}

\value{
 The function produces an sample of the posterior parameter value which is memorized in \code{mle.vam} object. If \code{history=FALSE}, the  method produces a list whose length is equal to the total number of parameters. Each element of this list corresponds to the sampled values for a parameter (in the same order as they appears in \code{par0}. If \code{history=TRUE}, the  method produces a data frame. Each accepted marginal sampled value of a parameter is alliteratively memorizes in this data frame. The column \code{ind} specifies the index of the corresponding parameter (starting from 0). The sampled values is in the column \code{estimate}. If \code{profile.alpha=TRUE}, a third column called \code{alpha} provides the corresponding estimate of \eqn{\alpha}.
 }

 \author{R. Drouilhet et L. Doyen}

\seealso{
\code{\link{bayesian.vam}} to define the Bayesian object.
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
run(bayesARAInf,profile.alpha=TRUE)
coef(bayesARAInf)
coef(bayesARAInf,par0=c(1e-2,2.5,0.5),fixed=2)
}

\keyword{virtual age}
\keyword{imperfect maintenance}
\keyword{corrective maintenance}
\keyword{preventive maintenance}
\keyword{ARA}
\keyword{Bayesian}

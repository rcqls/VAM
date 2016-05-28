\name{bayesian.vam}
\alias{bayesian.vam}

\title{Virtual age model for Bayesian estimation}

\description{\code{bayesian.vam} is used to define a virtual age model for Corrective Maintenance (CM) and planned Preventive Maintenance (PM). The object define with \code{bayesian.vam} can be used to compute the Maximum Likelihood Estimator (MLE) of the parameters thanks to the \code{\link{run.bayesian.vam}} method.}

\usage{
bayesian.vam(formula, data)
}

\arguments{
  \item{\code{formula}}{a symbolic description of the virtual age model and observations, or a \code{\link{bayesian.vam}} class object for which the estimation method has been launched at least one time. When \code{formula} is \code{bayesian.vam} object, the model considered corresponds to the plug in estimator, that is to say the output of the \code{\link{formula.bayesian.vam}} function. Otherwise, the details of formula specifications are given under `Details'.}
  \item{\code{data}}{a data frame or possibly a list (when several system are considered together) containing the observations.}
}

\details{
The symbolic description of the model done in \code{formula} has the form \code{response ~ model}.
\code{response} is a symbolic description of the data considered. The specifications are the same as those of \code{\link{model.vam}} function.
\code{model} is a symbolic description of the virtual age model considered. The specifications are the same as those of \code{\link{sim.vam}} function. In this case the PM policy is useless, so it has not to be necessarily defined. The parameter values specify in \code{model} for the maintenance effect models, and the time to failure distribution of the new unmaintained system, are used as initialization values for the first run of the likelihood maximization method.
}

\value{
 The function produces an object of class \code{bayesian.vam} which contains the virtual age model considered and the corresponding observations.
 }

 \author{R. Drouilhet}

\seealso{
\code{\link{run.bayesian.vam}} to compute the Bayesian method.
\code{\link{coef.bayesian.vam}} to extract the parameters value of the Bayesian method.
%\code{\link{formula.bayesian.vam}} to extract the original and estimated model.
%\code{\link{plot.bayesian.vam}} for plotting characteristics of the model.
}

\examples{
simARAInf<-sim.vam(  ~ (ARAInf(.4) | Weibull(.001,2.5)))
simData<-simulate(simARAInf,30)
bayesARAInf <- bayesian.vam(Time & Type ~ (ARAInf(~Unif(0,1)) | Weibull(~Unif(1,1.5),~Unif(2,4))),data=simData)
coef(bayesARAInf)

}

\keyword{virtual age}
\keyword{imperfect maintenance}
\keyword{corrective maintenance}
\keyword{preventive maintenance}
\keyword{ARA}

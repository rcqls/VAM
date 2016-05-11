\name{mle.vam}
\alias{mle.vam}

\title{Virtual age model for maximum likelihood estimation}

\description{\code{mle.vam} is used to define a virtual age model for Corrective Maintenance (CM) and planned Preventive Maintenance (PM). The object define with \code{mle.vam} can be used to compute the Maximum Likelihood Estimator (MLE) of the parameters thanks to the \code{\link{run.mle.vam}} method.}

\usage{
mle.vam(formula, data)
}

\arguments{
  \item{\code{formula}}{a symbolic description of the virtual age model and observations, or a \code{\link{mle.vam}} class object for which the estimation method has been launched at least one time. When \code{formula} is \code{mle.vam} object, the model considered corresponds to the plug in estimator, that is to say the output of the \code{\link{formula.mle.vam}} function. Otherwise, the details of formula specifications are given under `Details'.}
  \item{\code{data}}{a data frame or possibly a list (when several system are considered together) containing the observations.}
}

\details{
The symbolic description of the model done in \code{formula} has the form \code{response ~ model}.
\code{response} is a symbolic description of the data considered. The specifications are the same as those of \code{\link{model.vam}} function.
\code{model} is a symbolic description of the virtual age model considered. The specifications are the same as those of \code{\link{sim.vam}} function. In this case the PM policy is useless, so it has not to be necessarily defined. The parameter values specify in \code{model} for the maintenance effect models, and the time to failure distribution of the new unmaintained system, are used as initialization values for the first run of the likelihood maximization method.
}

\value{
 The function produces an object of class \code{mle.vam} which contains the virtual age model considered and the corresponding observations.
 }

 \author{L. Doyen and R. Drouilhet}

\seealso{
\code{\link{run.mle.vam}} to compute the MLE.
\code{\link{coef.mle.vam}} to extract the parameters value of the MLE.
\code{\link{formula.mle.vam}} to extract the original and estimated model.
\code{\link{plot.mle.vam}} for plotting characteristics of the model.
\code{\link{update.mle.vam}} to change the associated data set.
\code{\link{contrast.mle.vam}} to compute the contrast associated to the MLE.
\code{\link{logLikelihood.mle.vam}} to compute the log-likelihood.}

\examples{
simARAInf<-sim.vam(  ~ (ARAInf(.4) | Weibull(.001,2.5)))
simData<-simulate(simARAInf,30)
mleARAInf <- mle.vam(Time & Type ~ (ARAInf(.5) | Weibull(1,3)),data=simData)
coef(mleARAInf)

simCMPM_Multi<-sim.vam(  ~ (ARAInf(.3) | Weibull(.001,2.5)) & (ARAInf(.6)+ARAInf(-.2) | Periodic(12,prob=c(0.6,0.4))))
simData_Multi<-simulate(simCMPM_Multi,5000,nb.system=5)
mleCMPM_Multi <- mle.vam(System & Time & Type ~  (ARAInf(.5) | Weibull(1,3)) & (ARAInf(.5)+ARAInf(.5)),data=simData_Multi)
coef(mleCMPM_Multi)

}

\keyword{virtual age}
\keyword{imperfect maintenance}
\keyword{corrective maintenance}
\keyword{preventive maintenance}
\keyword{ARA}

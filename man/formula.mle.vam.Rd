\name{formula.model.vam}

\alias{formula.mle.vam}
\alias{formula}

\title{Return the plug-in and original model of \code{mle.vam} object}

\description{\code{formula} returns, for a \code{mle.vam} object, the original formula and the formula in which the estimated values of the parameters have been plug-in.}

\usage{
\method{formula}{mle.vam}(self,origin=FALSE)
}

\arguments{
\item{\code{self}}{an object of class \code{\link{mle.vam}}.}
\item{\code{origin}}{a logical which indicates if the original formula has also to be returned.}
}

\value{
If \code{origin==FALSE}, \code{formula} function returns a formula similar to the one used to construct the \code{self} object, but with parameter values which have been replaced by the ones previously estimated in the \code{self} object. Otherwise, it retruns a two dimentional list. The first argument \code{$formula} is the same plug-in estimated formula. Ans the second argument \code{$origin} is the original formula.
}


\author{L. Doyen and R. Drouilhet}

\seealso{
\code{\link{run.mle.vam}} to compute the MLE.
\code{\link{contrast.mle.vam}} to compute the contrast associated to the MLE.
\code{\link{coef.mle.vam}} to extract the parameters value of the MLE.
\code{\link{logLikelihood.mle.vam}} to compute the log-likelihood.
\code{\link{formula.mle.vam}} to extract the original and estimated model.
\code{\link{plot.mle.vam}} for plotting characteristics of the model.
\code{\link{update.mle.vam}} to change the associated data set.
}

\examples{
simCMPM<-sim.vam(  ~ (ARAm(.7|3) | Weibull(.001,2.5)) & (ARAInf(.4) | AtIntensity(0.2)))
simData=simulate(simCMPM,50)
mleCMPM <- mle.vam(Time & Type ~ (ARAm(.5|3) | Weibull(1,3)) & (ARAInf(0.5)),data=simData)
formula(mleCMPM)
coef(mleCMPM)
}


\keyword{virtual age}
\keyword{imperfect maintenance}
\keyword{corrective maintenance}
\keyword{preventive maintenance}
\keyword{ARA}

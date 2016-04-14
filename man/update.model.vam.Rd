\name{update.model.vam}

\alias{update.model.vam}
\alias{update.mle.vam}
\alias{update}

\title{Updates the data set in \code{model.vam} or \code{mle.vam} object}

\description{\code{update} changes the data set in \code{model.vam} or \code{mle.vam} object.}

\usage{
\method{update}{model.vam}(model,data)
\method{update}{mle.vam}(model,data)
}

\arguments{
\item{\code{model}}{an object of class \code{\link{model.vam}} or \code{\link{mle.vam}}.}
\item{\code{data}}{a data frame or possibly a list (when several system are considered together) containing the observations.}
}

\details{
The specifications of \code{data} must be compatible with the formula used to create the object \code{model}. See \code{\link{model.vam}} function for more details.
}

\author{L. Doyen and R. Drouilhet}

\examples{
# The model used to plan PM times
modMPplan <- model.vam(Time & Type ~ (ARA1(.87) | Weibull(.0015,2.6)) & (ARAInf(.44)))
# The model used for failure times simulation
simCMPM<-sim.vam( ~ (ARA1(.9) | Weibull(.002,2.5)) & (ARAInf(.25) | AtIntensity(0.35,modMPplan)) )
simData<-simulate(simCMPM,30)
# The simulated data set is updated in the model used for the PM policy
update(modMPplan,simData)
# The failureintensity of the failure times simulation model
plot(modMPplan,'i',ylim=c(0,0.5))
# And the intensity of the model used for the PM policy
plot(simCMPM,'i',col='blue',add=TRUE)
# The thershold at which PM are done
abline(h=0.35)
}


\keyword{virtual age}
\keyword{imperfect maintenance}
\keyword{corrective maintenance}
\keyword{preventive maintenance}
\keyword{ARA}

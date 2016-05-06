\name{model.vam}
\alias{model.vam}

\title{Virtual age model}

\description{\code{model.vam} is used to define a virtual age model for Corrective Maintenance (CM) and planned Preventive Maintenance (PM). The object define with \code{model.vam} can be used to plot characteristics of the CM-PM process.}

\usage{
model.vam(formula, data)
}

\arguments{
  \item{\code{formula}}{a symbolic description of the virtual age model and observations. The details of formula specifications are given under `Details'.}
  \item{\code{data}}{a data frame or possibly a list (when several system are considered together) containing the observations.}
}

\details{
The symbolic description of the model done in \code{formula} has the form \code{response ~ model}.

If only one system is considered, \code{response} has the form \code{Time & Type}. \code{Time} and \code{Type} are the column names of \code{data} containing respectively the successive maintenance times and the corresponding maintenance types. CM type is denoted -1. The successive different types of PM effects defined in \code{model} are denoted 1, 2, ... Censorship time type is denoted 0.

If several systems are considered together and \code{data} is a data frame, \code{response} has the form \code{System & Time & Type}. \code{Time} and \code{Type} have the same specifications as previously. \code{System} is the column names of \code{data} containing for each maintenance the system to which it refers. The different system are denoted 1, 2, ...

If several systems are considered together and \code{data} is a list, \code{response}  has the form \code{Time & Type}. And each term of the list \code{data} is a data frame corresponding to the observations of one different system. Each of these data frame has a column \code{Time} and \code{Type} with the same specification has previously.

\code{model} is a symbolic description of the virtual age model considered. The specifications are the same as those of \code{\link{sim.vam}} function. In this case the PM policy is useless, so it has not to be necessarily defined.
}

\value{
 The function produces an object of class \code{model.vam} which contains the virtual age model considered and the corresponding observations.
 }

 \author{L. Doyen and R. Drouilhet}

\seealso{\code{\link{plot.model.vam}} for plotting characteristics of the model.
\code{\link{update.model.vam}} to change the associated data set.}

\examples{
simARAInf<-sim.vam(  ~ (ARAInf(.4) | Weibull(.001,2.5)))
simData<-simulate(simARAInf,30)
plot(simARAInf,'i')#for plotting the failure intensity of the model used for simulating the maintenances
# But, the model associated to simData can be changed
modelAGAN<-model.vam(Time & Type ~ (AGAN() | Weibull(.001,2.5)),data=simData)
plot(modelAGAN,'i',col='darkblue',add=TRUE)
}
\keyword{virtual age}
\keyword{imperfect maintenance}
\keyword{corrective maintenance}
\keyword{preventive maintenance}
\keyword{ARA}

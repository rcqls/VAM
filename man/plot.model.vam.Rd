\name{plot.model.vam}

\alias{plot.model.vam}
\alias{plot.sim.vam}
\alias{plot.mle.vam}
\alias{plot}

\title{Plot characteristics of virtual age model}

\description{\code{plot} plots characteristics of a virtual age model for Corrective Maintenance (CM) and planned Preventive Maintenance (PM).}

\usage{
\method{plot}{model.vam}(obj, type, from, to, length.out=101, by, system.index=1 ,cm.type=NA, pm.type=NA, add=FALSE, ....)
\method{plot}{sim.vam}(obj, type, from, to, length.out=101, by, system.index=1 ,cm.type=NA, pm.type=NA, add=FALSE, ....)
\method{plot}{mle.vam}(obj, type, from, to, length.out=101, by, system.index=1 ,cm.type=NA, pm.type=NA, add=FALSE, ....)
}

\arguments{
\item{\code{obj}}{an object of class \code{\link{model.vam}} or \code{\link{sim.vam}} or \code{\link{mle.vam}} specifying the virtual age model considered. This object must contain a data set. Then, in the case of a \code{\link{sim.vam}} object, the \code{\link[VAM]{simulate}} method must have been previously applied, at least one time, on the object.}
\item{\code{type}}{specifies which characteristic to plot:
\itemize{
	\item "i" or "intensity" for failure intensity,
	\item "I" or "cumulative" for cumulative failure intensity also called compensator of the CM-PM process,
	\item "v" or "virtual.age" for virtual age.
}
By default, PM times are also represented, but "-pm" can be added at the end of the string "\code{type}" in order to avoid this representation. By default, CM times are also represented (or the cumulative number of CM in the "I" case), but "-cm" can be added at the end of the string "\code{type}" in order to avoid this representation. For example, the types \code{"I-pm"} and \code{"I-pm-cm"} both correspond to the cumulative failure intensity plot, but in the first case the cumulative number of failure will also be plotted.
}
\item{\code{from, to}}{optional arguments specifying the range time over which the characteristic will be plotted.}
\item{\code{by}}{time increment between two successive  maintenance time at which the characteristic will be plotted.}
\item{\code{length.out}}{When provided, \code{by} is defined as \code{(to - from)/(length.out - 1)}.}
\item{\code{system.index}}{the index of the system for which to plot the characteristic.}
\item{\code{cm.type}}{how to additionally represent the CM times or the cumulative number of CM in the "I" case. Possible types are "p" for points and "l" for lines.}
\item{\code{pm.type}}{how to additionally represent the PM times. Possible types are "p" for points and "l" for lines.}
\item{\code{add}}{if "\code{add==FALSE}" the characteristic is plotted on a new graphical, otherwise it is added to the current one.}
\item{\code{...}}{Further classical graphical parameters specifying how the plot of the characteristic. Others non usual arguments can also be added in order to specify the additional representation of the CM and PM times: \code{cm.col, cm.pch, cm.lty, cm.lwd, pm.col, pm.pch, pm.lty, pm.lwd}. For PM, those arguments can be vectors in order to specify a specific value of the argument for each PM type.}
}

\author{L. Doyen and R. Drouilhet}

\examples{
simARAInf<-sim.vam(  ~ (ARAInf(.4) | Weibull(.05,2.5)))
simData<-simulate(simARAInf,10,nb.system=50,as.list=TRUE)
plot(simARAInf,'i',system.index=2)
mleARAInf <- mle.vam(Time & Type ~ (ARAInf(0.5) | Weibull(1,3)),data=simData)
run.mle.vam(mleARAInf)
plot(mleARAInf,'i-cm',col='blue',add=TRUE,system.index=2)

simCMPM<-sim.vam(  ~ (ARAInf(.3) | Weibull(.001,2.5)) & (ARAInf(.6)+ARAInf(-.2) | Periodic(12,prob=c(0.6,0.4))))
(simData<-simulate(simCMPM,50))
plot(simCMPM,'I',from=simData$Time[10],to=simData$Time[40],col='darkblue',xlab='t',ylim=c(0,sum(simData$Type[1:40]==-1)),cm.col='red',pm.type='p',pm.pch=c(1,2),pm.col=c("green","black"),main="Fantastic!")
}


\keyword{virtual age}
\keyword{imperfect maintenance}
\keyword{corrective maintenance}
\keyword{preventive maintenance}
\keyword{ARA}

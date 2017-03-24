\name{plot.bayesian.vam}
\alias{plot.bayesian.vam}
\alias{plot}

\title{Plot the Bayesian estimtions of characteristics of a virtual age model}

\description{For \code{bayesian.vam} object, plots the Bayesian estimations (mean and 0.95 credibility interval) of the characteristics of a virtual age model for Corrective Maintenance (CM) and planned Preventive Maintenance (PM).}

\usage{plot.bayesian.vam(obj,type,from,to,length.out=101,by,system.index=1,cm.type=NA,pm.type=NA,add=FALSE,nb.proposal=500,...)
}

\arguments{
	\item{\code{obj}}{an object of class \code{\link{bayesian.vam}} specifying the virtual age model considered, the prior distribution of the parameters and the data set of observations.}
	\item{\code{type}}{specifies which characteristic to plot:
\itemize{
	\item "i" or "intensity" for failure intensity,
	\item "I" or "cumulative" for cumulative failure intensity also called compensator of the CM-PM process,
	\item "v" or "virtual.age" for virtual age,
	\item "F" or "conditional.cdf" for conditional distribution function of the next failure time,
	\item "S" or "conditional.survival" for conditional survival function of the next failure time,
	\item "f" or "conditional.pdf" for conditional probability density function of the next failure time.
	\item "d" or "data" for data representation as the CM counting process.
}
By default, PM times are also represented, but "-pm" can be added at the end of the string \code{type} in order to avoid this representation. By default, CM times are also represented (or the cumulative number of CM in the "I" case), but "-cm" can be added at the end of the string "\code{type}" in order to avoid this representation. For example, the types \code{"I-pm"} and \code{"I-pm-cm"} both correspond to the cumulative failure intensity plot, but in the first case the cumulative number of failure will also be plotted.
}
\item{\code{from, to}}{optional arguments specifying the range time over which the characteristic will be plotted.}
\item{\code{by}}{time increment between two successive  maintenance time at which the characteristic will be plotted.}
\item{\code{length.out}}{When provided, \code{by} is defined as \code{(to - from)/(length.out - 1)}.}
\item{\code{system.index}}{the index of the system for which to plot the characteristic.}
\item{\code{cm.type}}{how to additionally represent the CM times or the cumulative number of CM in the "I" case. Possible types are "p" for points and "l" for lines.}
\item{\code{pm.type}}{how to additionally represent the PM times. Possible types are "p" for points and "l" for lines.}
\item{\code{add}}{if \code{add==FALSE} the characteristic is plotted on a new graphical, otherwise it is added to the current one.}
  	\item{nb.proposal}{an optional argument specifying the size of the simulated sample of the posterior distribution (used to compute the estimation and credibility interval of the plotted characteristic).}
\item{\code{...}}{Further classical graphical parameters specifying how the plot of the characteristic. Others non usual arguments can also be added in order to specify the additional representation of the CM and PM times: \code{cm.col, cm.pch, cm.lty, cm.lwd, pm.col, pm.pch, pm.lty, pm.lwd}. For PM, those arguments can be vectors in order to specify a specific value of the argument for each PM type.}
}

\details{
	If the \code{\link{run.bayesian.vam}} method has already been executed with arguments \code{history=TRUE} and \code{nb>=nb.proposal}, it is not executed again and the simulated sample of the posterior distribution memorized in \code{obj} is used (only the \code{nb.proposal} first simulated values are used). Otherwise the \code{\link{run.bayesian.vam}} method is executed.
}

 \author{R. Drouilhet et L. Doyen}

\seealso{
\code{\link{bayesian.vam}} to define the Bayesian object.
\code{\link{run.bayesian.vam}} to compute the Bayesian method.
\code{\link{summary.bayesian.vam}} to produce a result summary of the Bayesian method.
\code{\link{coef.bayesian.vam}} to extract the parameters estimation values of the Bayesian method.
%\code{\link{formula.bayesian.vam}} to extract the original and estimated model.
\code{\link{hist.bayesian.vam}} for plotting the histogram of the posterior distribution of the parameters.
}

\examples{
simARAInf<-sim.vam(  ~ (ARAInf(.4) | Weibull(.001,2.5)))
simData<-simulate(simARAInf,30)
bayesARAInf <- bayesian.vam(Time & Type ~ (ARAInf(~Unif(0,1)) | Weibull(~Unif(0,1),~Unif(2,4))),data=simData)
run(bayesARAInf,profile.alpha=TRUE,history=TRUE)
plot(bayesARAInf,'i')
}

\keyword{virtual age}
\keyword{imperfect maintenance}
\keyword{corrective maintenance}
\keyword{preventive maintenance}
\keyword{ARA}
\keyword{Bayesian}
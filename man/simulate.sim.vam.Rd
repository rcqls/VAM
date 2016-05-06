\name{simulate.sim.vam}
\alias{simulate.sim.vam}
\alias{simulate}

\title{Simulation of virtual age model}

\description{\code{simulate} generates a realization of Corrective Maintenance (CM) and planned Preventive Maintenance (PM) times of a virtual age model. The last simulated data set is also memorized inside the \code{\link{sim.vam}} object used for specifying the simulation model.}

\usage{
\method{simulate}{sim.vam}(sim, stop.policy = 10, nb.system=1, cache.size=500, as.list=FALSE, data=NULL)
}

\arguments{
  \item{\code{sim}}{an object of class \code{\link{sim.vam}} specifying the virtual age model used for simulation.}
  \item{\code{stop.policy}}{indicates when to stop the simulation. \code{stop.policy} can be a scalar, then the simulation stop as soon as the number of maintenance (CM or PM) times per system (per CM-PM trajectory) reaches \code{stop.policy}. Otherwise \code{stop.policy} can be a symbolic formula, details of formula specifications are given under `Details'.}
  \item{\code{nb.system}}{the number of parallel systems (trajectories of the CM-PM process) to simulate.}
  \item{\code{cache.size}}{if the number of events to generate can not be directly determined thanks to the \code{stop.policy} and need to be dynamically allocated, \code{cache.size} specifies the length of each size reallocation of the table of results.}
  \item{\code{as.list}}{by default the function returns a data frame. But, if \code{nb.system>1} and \code{as.list==TRUE}, the result returned is a list. The inside code of the function use a list, then if \code{as.list==TRUE} a data set class conversion is avoid and the executing time of the function can be faster when \code{nb.system} is large.}
  \item{\code{data}}{NOT YET IMPLEMENTED ! but the idea is to simulate the future of the data set \code{data} if this one is specified.}
}

\details{
Successive maintenance time simulation is done iteratively on each system.
A symbolic formula can be written in \code{stop.policy} in order to describe when to stop the simulation. This criteria is applied to each simulated system. The symbolic formula of \code{stop.policy} can combined  as many different stop conditions as required thanks to the operators \code{|} and \code{&}:
\describe{
	\item{\code{cond1 | cond2}}{stops as soon as \code{cond1} or \code{cond2} is verified.}
	\item{\code{cond1 & cond2}}{stops as soon as \code{cond1} and \code{cond2} are verified.}
}
The stopping conditions available for \code{cond1}, \code{cond2} and so on are of the form:
\describe{
	\item{\code{Time>a}}{stops as soon as the last simulated maintenance time is greater than \code{a}.}
	\item{\code{T>a}}{shortcut for \code{Time>a}.}
	\item{\code{Size==a}}{stops as soon as the total number of simulated maintenance times is equal to \code{a}.}
	\item{\code{S==a} or \code{Size>=a} or \code{S>=a}}{shortcut for \code{Size==a}.}
	\item{\code{Size[b]==a}}{stops as soon as the total number of simulated maintenance times of type \code{b}
	is equal to \code{a}. CM type is denoted -1. And the successive different types of PM effects defined in the formula of \code{SimModel} are denoted 1, 2, ...}
	\item{\code{S[b]==a} or \code{Size[b]>=a} or \code{S[b]>=a} }{shortcut for \code{Size[b]==a}.}
	\item{\code{Time>(RightCensorship=a)}}{stops as soon as the last simulated maintenance time is greater than \code{a}. In addition, all the possible simulated maintenance time greater than \code{a} (at least the last maintenance time simulated) are censored at time \code{a}: their times are equal to \code{a} and their types to 0. In this case, \code{a} can also be replaced by a random time distribution name, with a capital first letter and followed by its arguments: for example \code{Unif(20,30)}. Then, the censorship times of each system will be simulated with the corresponding random distribution, in the previous example \code{runif(nb.system,20,30)}.}
	\item{\code{T>(RC=a)}}{shortcut for \code{Time>(RightCensorship=a)}.}
}
}

 \value{
 By default, the function produces one data frame containing all the simulated maintenances. Each line
corresponds to a maintenance. If \code{nb.system==1}, the generated data frame has two columns. The first one, called \code{Time}, indicates the successive maintenance times. The second one, called \code{Type}, indicates the corresponding maintenance types. CM type is denoted -1. The successive different types of PM effects defined in the formula of \code{SimModel} are denoted 1, 2, ... Censorship time type is denoted 0. If \code{nb.system>1}, a column is added at the beginning. Its name is \code{System}, and it precises for each event to which system it refers to. The different system are denoted 1, 2, ...

If \code{as.list=FALSE} and  \code{nb.system>1}, the function produces a list. Each element of the list corresponds to a different system and is a data frame with the two columns \code{Time} and \code{Type}.
 }

\author{L. Doyen and R. Drouilhet}

\seealso{\code{\link{sim.vam}} for model definition.}

\examples{
simARAInf<-sim.vam(  ~ (ARAInf(.4) | Weibull(.001,2.5)))
(simulate(simARAInf,Time>25))
(simulate(simARAInf,Time>(RightCensorship=25)))
(simData_Multi<-simulate(simARAInf,T>(RC=Unif(20,30)),nb.system=5))

(simulate(simARAInf,Time>25 | Size>5))
time_max<-100
size_max<-5.5
(simulate(simARAInf,Time>time_max|Size>size_max))
(simulate(simARAInf,Time>25 & Size>5))
(simulate(simARAInf,T>100 & S>5))
(simulate(simARAInf,T>(RC=25) & S>5))

simCMPM<-sim.vam(  ~ (ARA1(.9) | Weibull(.001,2.5)) & (ARAInf(.4) | AtIntensity(0.2)))
simData<-simulate(simCMPM,Size[-1]>10&Size[1]>5)
table(simData$Type)
simData<-simulate(simCMPM,Size[-1]>10|Size[1]>5)
table(simData$Type)
simData<-simulate(simCMPM,(Size[-1]>10&Size[1]>5)|Size>15)
table(simData$Type)
}


\keyword{virtual age}
\keyword{imperfect maintenance}
\keyword{corrective maintenance}
\keyword{preventive maintenance}
\keyword{ARA}

\name{sim.vam}
\alias{sim.vam}

\title{Virtual age simulation model}

\description{\code{sim.vam} is used to define a virtual age model for Corrective Maintenance (CM) and planned Preventive Maintenance (PM). The object define with \code{sim.vam} can be used to simulate realizations of the CM-PM process thanks to the \code{\link{simulate}} function. The object define with \code{sim.vam} can also be used to define a PM policy.}

\usage{
sim.vam(formula)
}

\arguments{
  \item{\code{formula}}{ specifies the virtual age model considered. For simulation, with \code{\link{simulate}} function, "\code{formula}" has necessarily one of following form:
\describe{
	\item{Only CM:}{\code{~(CMeff|InitDist)}}
	\item{CM and one type of PM effect:}{\code{~(CMeff|InitDist)&(PMeff|PMpolicy)}}
	\item{CM and several types of PM effect:}{\code{~(CMeff|InitDist)&(PMeff_1+...+PMeff_n|PMpol)} 
	
	\code{~(CMeff|InitDist)&(PMeff_1+...+PMeff_n|PMpol_1*...*PMpol_k)}}
}
In others cases than simulation (where the PM policy is useless), "\code{formula}" can possibly be simplify to:
 \itemize{
	\item \code{~(CMeff|InitDist)&(PMeff)}
	\item \code{~(CMeff|InitDist)&(PMeff_1+...+PMeff_n)} 
	\item \code{~(CMeff|InitDist)&(PMeff_1+...+PMeff_n)}
}
 }}
 
\details{
From a mathematical point of view the CM-PM process can be specified by its failure intensity. The failure intensity of a virtual age model verifies for \eqn{t>=0}
\deqn{\lambda(t)=h(V(t))}
where \eqn{h} is a deterministic function characterizing the time to failure hazard rate of the new unmaintained system (specified by "\code{InitDist}" in "\code{formula}"). And \eqn{V} is an adapted process, called virtual age, that takes into account the effects of the different maintenances. In order to be able to simulate the CM-PM process, one has also to specify how PM are planned. This is done thanks to the PM policy denoted "\code{PMpol}" or "\code{PMpol_1}", ...., "\code{PMpol_k}" in "\code{formula}".

\subsection{Maintenance effects}{
"\code{CMeff}", "\code{PMeff}", "\code{PMeff_1}", ....,, "\code{PMeff_n}" specified respectively the effect of CM, of PM (where only one type of PM effect is considered) and of the different types of PM effect. Notice that only one type of CM effect can be considered. The available maintenance effects are :
\describe{
	\item{\code{AGAN()}}{for As Good As New maintenance. An AGAN maintenance renews the system. That is to say,
	after an AGAN maintenance done at time \eqn{C[k]}, the virtual age is equal to \eqn{V(t)=t-C[k]} up to the next maintenance time \eqn{C[k+1]} (for \eqn{C[k]<t<= C[k+1]}).}
	\item{\code{ABAO()}}{for As Bad As Old maintenance. An ABAO maintenance has no effect on the propensity of the system to break done, it continues to evolve similarly to what it was before maintenance. That is to say,
	after an ABAO maintenance  done at time \eqn{C[k]}, the virtual age is equal to \eqn{V(t)=t-C[k]+V(C[k])} up to the newt maintenance time \eqn{C[k+1]} (for \eqn{C[k]<t<= C[k+1]}).}
	\item{\code{ARAInf(rho)}}{for Arithmetic Reduction of Age model with infinite memory. This is also equivalent to a Kijima type II model with deterministic constant effect. The effect of such maintenance is to reduce the virtual age from a quantity proportional to its value just before maintenance. That is to say, after an ARAInf maintenance  done at time \eqn{C[k]}, the virtual age is equal to
	 \eqn{V(t)=t-C[k]+(1-\rho)V(C[k])} up to the newt maintenance time \eqn{C[k+1]} (for \eqn{C[k]<t<= C[k+1]}). "\code{rho}" is a scalar smaller or equal to 1, representing the maintenance effect \eqn{\rho}.}
	\item{\code{ARA1(rho)}}{for Arithmetic Reduction of Age model with memory one. This also equivalent to a Kijima type I model with deterministic constant effect. The effect of such maintenance is to reduce the virtual age from a quantity proportional to the supplement of age accumulated since the last maintenance. That is to say, after an ARA1 maintenance  done at time \eqn{C[k]} (\eqn{C[k-1]} denotes the previous maintenance time), the virtual age is equal to
	\eqn{V(t)=t-C[k]+(1-\rho)(V(C[k])-V(C[k-1]))} up to the newt maintenance time \eqn{C[k+1]} (for \eqn{C[k]<t<= C[k+1]}). "\code{rho}" is a scalar smaller or equal to 1, representing the maintenance effect \eqn{\rho}.}
	}}

\subsection{Time to failure distribution of the new unmaintained system}{
	"\code{InitDistr}" characterizes the hazard rate of the new unmaintained system \eqn{h(.)}. The available distributions are:
	\describe{
	\item{\code{Weibull(alpha,beta)}}{for a Weibull initial hazard rate with parametrization \eqn{h(t)=\alpha\beta t^(\beta-1)} with \eqn{\alpha>0} and \eqn{\beta>0}.}
	\item{\code{LogLinear(alpha,beta)}}{for a log-linear hazard rate with parametrization \eqn{h(t)=\alpha exp(\beta t)} with \eqn{\alpha>0}.}
	}
}

\subsection{PM policies}{
	"\code{PMpol}" or "\code{PMpol_1}", ...., "\code{PMpol_k}" characterize the PM policy. Some of the PM policies available can manage simultaneously several PM effect types. That is why the number of PM effect arguments ("\code{PMeff_1}", ...., "\code{PMeff_n}") in "\code{formula}" is not necessary the same as the number of PM policy arguments ("\code{PMpol_1}", ...., "\code{PMpol_k}"). Whatever the first PM policy argument refers to the corresponding first PM effect arguments and so on. The different PM policies available are:
	\describe{
	\item{\code{periodic(by,from=0,prob=1)}}{a PM is done every "\code{by}" units of times from time "\code{from}". This PM policy can manage \eqn{k>1} different PM effect types if "\code{prob}" is a vector of length \eqn{k}. In this case, each PM time as a probability \code{prob[j]} to be of \eqn{j} th type (for the different PM effect considered) for every \eqn{j} between 1 and \eqn{k}.}
	\item{\code{AtIntensity(level,model=NULL)}}{a PM is done as soon as the failure intensity reaches "\code{level}". By default "\code{level}" refers to the failure intensity computed with the model defined by "\code{formula}". But, it can also refer to another model denoted "\code{model}". "\code{model}" is a \code{sim.vam} object depending of another formula denoted "\code{formula_bis}". The number of different PM effects has to be the same in "\code{formula}" and "\code{formula_bis}". The grammar used to write "\code{formula_bis}" is the same as this of "\code{formula}". But in "\code{formula_bis}" the PM policy is useless, so it has not to be necessarily defined. For example in the case where only one type of PM effect is considered in "\code{formula}", then "\code{formula_bis}" can equivalently follows the form : \code{~(CMeff_bis|InitDist_bis)&(PMeff_bis|PMpolicy_bis)} or \code{~(CMeff_bis|InitDist_bis)&(PMeff_bis)}.}
	\item{\code{AtVirtualAge(level,model=NULL)}}{a PM is done as soon as the virtual age reaches "\code{level}". By default "\code{level}" refers to the virtual age computed with the model defined by "\code{formula}". But, it can also refers to another model similarly to "\code{AtIntensity}" PM policy.}
	\item{\code{AtFailureProbability(level,model=NULL)}}{a PM is done as soon as the conditional probability of failure reaches "\code{level}". By default "\code{level}" refers to the conditional failure probability computed with the model defined by "\code{formula}". But, it can also refers to another model similarly to "\code{AtIntensity}" PM policy.}
	}
}
 }

 \value{
 The function produces an object of class "\code{sim.vam}" which contain the virtual age model considered.
 }

 \references{Classes of imperfect repair models based on reduction of failure intensity or virtual age:
L. Doyen, O. Gaudoin,
Reliability Engineering and System Safety, Elsevier, 2004, 84 (1), pp.45-56.

Modelling and assessment of aging and efficiency of corrective and planned preventive maintenance:
L. Doyen, O. Gaudoin, IEEE Transactions on Reliability, 2011, 60 (4), pp.759-769.
}

\author{L. Doyen and R. Drouilhet}

\seealso{\code{\link{simulate}} for simulation.}


\examples{
##########
# Simulation model:
# ARA Infinite CM with \rho=0.4 
# Weibull initial intensity h(t)=0.001*2.5*t^(1.5)
simARAInf<-sim.vam(  ~ (ARAInf(.4) | Weibull(.001,2.5)))
simData<-simulate(simARAInf,30)

##########
# Simulation model:
# Weibull h(t)=0.001*2.5*t^(1.5)
# CM ARA1 (rhoMC=0.6)
# PM ARA Infinite (rhoMP=0.4)
# PM at fixed failure intensity level
simCMPM<-sim.vam(  ~ (ARA1(.9) | Weibull(.001,2.5)) & (ARAInf(.4) | AtIntensity(0.2)))
simData=simulate(simCMPM,50)

###########
# Simulation model:
# Weibull h(t)=0.001*2.5*t^(1.5)
# CM ARA infinite (rhoMC=0.3)
# periodic PM randomly chosen upon two types
# PM ARA Infinite (rhoMP=0.6) with probability 0.6
# PM ARA Infinite (rhoMP=-0.2) with probability 0.4
simCMPM_Multi<-sim.vam(  ~ (ARAInf(.3) | Weibull(.001,2.5)) & (ARAInf(.6)+ARAInf(-.2) | Periodic(12,prob=c(0.6,0.4))))
simData<-simulate(simCMPM_Multi,20,nb.system=50)

##########
# Combined PM policies:
# Weibull h(t)=0.001*2.6*t^(1.6)
# CM ABAO
# periodic PM randomly chosen upon two types
# PM ARA Infinite (rhoMP=0.6) with probability 0.6
# PM ARA Infinite (rhoMP=-0.2) with probability 0.4
# and AGAN PM renew the system as soon as the failure intensity reaches the level 0.3
simCMPM<-sim.vam(  ~ (ABAO() | Weibull(.001,2.6)) & (ARAInf(.6)+ARAInf(-.2)+AGAN() | Periodic(12,prob=c(0.6,0.4))*AtIntensity(0.3)))
simData<-simulate(simCMPM,15)

###########
# PM policy not using the same model as for simulation:
# The model used for PM planning
modMPplan <- model.vam( ~ (ARA1(.87) | Weibull(.0015,2.7)) & (ARAInf(.44)))
# The model used for failure time simulation
# a PM is done as soon as the failure intensity of the previous model reaches level 0.3 
simCMPM<-sim.vam( ~ (ARA1(.9) | Weibull(.001,2.5)) & (ARAInf(.35) | AtIntensity(0.3,modMPplan)) )
simData<-simulate(simCMPM,20)
modSim_plot<-model.vam(Time & Type ~ (ARA1(.9) | Weibull(.001,2.5)) & (ARAInf(.35)),data=simData)
modMPplan_plot<-model.vam(Time & Type ~ (ARA1(.87) | Weibull(.0015,2.7)) & (ARAInf(.44)),data=simData)
# The failure intensity of the model used for PM planning
# a PM is done when this failure intensity reaches the level 0.3
plot(modMPplan_plot,'h',col='blue')
# The failure intensity of the model used for failure simulation
plot(modSim_plot,'h',col='darkblue',add=TRUE)
}


\keyword{virtual age}
\keyword{imperfect maintenance}
\keyword{corrective maintenance}
\keyword{preventive maintenance}
\keyword{ARA}
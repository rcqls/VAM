\name{sim.vam}
\alias{sim.vam}

\title{Virtual age simulation model}

\description{\code{sim.vam} is used to define a virtual age model for Corrective Maintenance (CM) and planned Preventive Maintenance (PM). The object define with \code{sim.vam} can be used to simulate realizations of the CM-PM process thanks to \code{\link[VAM]{simulate}} method. The last simulated data set, produced with the \code{\link[VAM]{simulate}} method, is also memorized inside the \code{sim.vam} object, in such way that the characteristics of the model can be plotted thanks to \code{\link[VAM]{plot}} method. The object define with \code{sim.vam} can also be used to define a PM policy.}

\usage{
sim.vam(formula)
}

\arguments{
  \item{\code{formula}}{a symbolic description of the virtual age model. 
For simulation, with \code{\link[VAM]{simulate}} method, "\code{formula}" has necessarily one of following form:
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
More details of formula specifications are given under `Details'.
 }}
 
\details{
From a mathematical point of view the CM-PM process can be specified by its failure intensity. The general failure intensity of a virtual age model verifies for \eqn{t >= 0}
\deqn{\lambda(t) = V'(t) h(V(t))}
where \eqn{h} is a deterministic function characterizing the time to failure hazard rate of the new unmaintained system (specified by "\code{InitDist}" in "\code{formula}"). \eqn{V(.)} is an adapted process, called virtual age, that takes into account the effects of the different maintenances, and \eqn{V'(.)} denotes its derivative. Let us denote \eqn{C[k]} the successive maintenance times and \eqn{V(t | k)} the virtual age at time \eqn{t} given that \eqn{k} maintenance times have already been observed. \eqn{V(t | 0)} denotes the virtual age of a new never maintained system, it is assumed to be equal to \eqn{t}. Consequently, the cumulative intensity or compensator of the CM-PM process, \eqn{\Lambda} verifies for \eqn{C[k] <= t <= C[k+1]},
\deqn{\Lambda(t)-\Lambda(C[k]) = H(V(t | k)) - H(V(C[k] | k))}
with \eqn{H(t) = \int_{0}^{t} h(s) ds}.

In order to be able to simulate the CM-PM process, one has also to specify how PM are planned. This is done thanks to the PM policy denoted "\code{PMpol}" or "\code{PMpol_1}", ...., "\code{PMpol_k}" in "\code{formula}". 

\subsection{Maintenance effects}{
"\code{CMeff}", "\code{PMeff}", "\code{PMeff_1}", ....,, "\code{PMeff_n}" specified respectively the effect of CM, of PM (where only one type of PM effect is considered) and of the different types of PM effect. Notice that only one type of CM effect can be considered. The available maintenance effects are :
\describe{
	\item{\code{AGAN()}}{for As Good As New maintenance. An AGAN maintenance renews the system. That is to say,
	after an AGAN maintenance done at time \eqn{C[k]}, the virtual age is equal to \eqn{V(t | k) = V(t-C[k] | 0) = t-C[k]} up to the next maintenance time \eqn{C[k+1]} (for \eqn{C[k] < t <= C[k+1]}).}
	\item{\code{ABAO()}}{for As Bad As Old maintenance. An ABAO maintenance has no effect on the propensity of the system to break done, it continues to evolve similarly to what it was before maintenance. That is to say,
	after an ABAO maintenance  done at time \eqn{C[k]}, the virtual age is equal to \eqn{V(t | k) = V(t | k-1)} up to the next maintenance time \eqn{C[k+1]} (for \eqn{C[k] < t <= C[k+1]}).}
	\item{\code{AGAP()}}{for As Good As Previous maintenance. An AGAP maintenance puts the system exactly in the same state as it was just after the last previous maintenance. That is to say, after an AGAP maintenance  done at time \eqn{C[k]}, the virtual age is equal to \eqn{V(t | k) = V(t-(C[k]-C[k-1]) | k-1)} up to the next maintenance time \eqn{C[k+1]} (for \eqn{C[k] < t <= C[k+1]}).}
	\item{\code{QAGAN()}}{for Quasi As Good As New maintenance. A QAGAN maintenance puts back the virtual age to zero without necessary renewing the system. That is to say, after a QAGAN maintenance  done at time \eqn{C[k]}, the virtual age is equal to \eqn{V(t | k) = V(t | k-1)-V(C[k] | k-1)} up to the next maintenance time \eqn{C[k+1]} (for \eqn{C[k] < t <= C[k+1]}).}
	\item{\code{ARAInf(rho)}}{for Arithmetic Reduction of Age model with infinite memory. This is also equivalent to a Kijima type II model with deterministic constant effect. The effect of such maintenance is to reduce the virtual age from a quantity proportional to its value just before maintenance. That is to say, after an ARAInf maintenance  done at time \eqn{C[k]}, the virtual age is equal to
	 \eqn{V(t | k) = V(t | k-1) - \rho V(C[k] | k-1)} up to the next maintenance time \eqn{C[k+1]} (for \eqn{C[k] < t <= C[k+1]}). "\code{rho}" is the corresponding maintenance effect parameter vector, it is of size 1 and  represents \eqn{\rho}.}
	\item{\code{ARA1(rho)}}{for Arithmetic Reduction of Age model with memory one. This also equivalent to a Kijima type I model with deterministic constant effect. The effect of such maintenance is to reduce the virtual age from a quantity proportional to the supplement of age accumulated since the last maintenance. That is to say, after an ARA1 maintenance  done at time \eqn{C[k]}, the virtual age is equal to
	\eqn{V(t | k) = V(t | k-1) - \rho ( V(C[k] | k-1) - V(C[k-1] | k-1) )} up to the next maintenance time \eqn{C[k+1]} (for \eqn{C[k] < t <= C[k+1]}). "\code{rho}" is the corresponding maintenance effect parameter vector, it is of size 1 and  represents \eqn{\rho}.}
	\item{\code{ARAm(rho) or ARAm(rho | m)}}{for Arithmetic Reduction of Age model with memory "\code{m}". In the ARAInf model, maintenance is supposed to reduce the whole virtual age, whereas in the ARA1 model, it is supposed to reduce only the supplement of age accumulated since the last maintenance. The idea of the ARAm model is to proposed an intermediate effect on the "\code{m}" last times between maintenances. Consequently, an ARAm model with \eqn{m=1} corresponds to an ARA1 model, and an ARAm model with \eqn{m} bigger than the number of maintenance time observations corresponds to an ARAInf model. But, the bigger \code{m} is, in the ARAm model, the less efficient (from an algorithmic point of view) are the different algorithms of the VAM package.
	The default value for \code{m} is 1. "\code{rho}" is the corresponding maintenance effect parameter vector, it is of size 1 and has the same meaning as for ARAInf and ARA1 models.}
	\item{\code{QR(rho)}}{for Quasi Renewal maintenance model. This also equivalent to the geometric model. Let \eqn{Y_1, Y_2, Y_3, ...} be iid random variables. The classical geometric process assumes that the successive times between failures follows the same distribution as \eqn{Y_1, a Y_2, a^2 Y_3, ...}. Consequently, after a QR maintenance  done at time \eqn{C[k]}, the virtual age is equal to
	\eqn{V(t | k) = \rho ( V(t | k-1) -  V(C[k] | k-1) )} up to the next maintenance time \eqn{C[k+1]} (for \eqn{C[k] < t <= C[k+1]}). A QR maintenance modifies the virtual age time speed. "\code{rho}" is the corresponding maintenance effect parameter vector, it is of size 1 and  represents \eqn{\rho = 1/a}.}
	\item{\code{GQR(rho)} or \code{GQR(rho | fun)}}{for Generalized Quasi Renewal maintenance model. This model is similar to the extended geometric process. The geometric evolution of the times between maintenances in the QR model often seems to much explosive to be realistic from a practical point of view. The idea of the GQR model is simply to replace \eqn{a^n} by \eqn{a^(f(n))} where \eqn{f()} stands for a non decreasing function that can be specified by the "\code{fun}" argument. The possible values for "\code{fun}" are 
	\itemize{
		\item{\code{identity}} : (default value) \eqn{f(x)=x}, in this case the GQR model corresponds to a QR model,
		\item{\code{sqrt}} : \eqn{f(x)=x^(0.5)},
		\item{\code{log}} : \eqn{f(x)=log(x+1)}.
	}
	A GQR maintenance modifies the virtual age time speed.
	"\code{rho}" is the corresponding maintenance effect parameter vector, it is of size 1.}
	\item{\code{GQR_ARAInf(rho)} or \code{GQR_ARAInf(rho | fun)}}{
	A GQR maintenance puts back the virtual age to zero and modifies the virtual age time speed.
	An ARAInf maintenance does not modify the virtual age time speed but restores the virtual age into an intermediate condition between its value before maintenance and zero. The GQR_ARAInf model mixes both effects.
	"\code{rho}" is the corresponding maintenance effect parameter vector, it is of size 2 : "\code{rho[1]}" characterizes the GQR effect and "\code{rho[2]}" characterizes the ARAInf effect.}
	\item{\code{GQR_ARA1(rho)} or \code{GQR_ARA1(rho | fun)}}{
	A GQR maintenance puts back the virtual age to zero and modifies the virtual age time speed.
	An ARA1 maintenance does not modify the virtual age time speed but restores the virtual age into an intermediate condition between its value before maintenance and zero. The GQR_ARA1 model mixes both effects.
	"\code{rho}" is the corresponding maintenance effect parameter vector, it is of size 2 : "\code{rho[1]}" characterizes the GQR effect and "\code{rho[2]}" characterizes the ARA1 effect.}
	\item{\code{GQR_ARAm(rho)} or \code{GQR_ARAm(rho | fun)} or \code{GQR_ARAm(rho | m)} or \code{GQR_ARAm(rho | m, fun)} or \code{GQR_ARAm(rho | fun, m)}}{
	The GQR_ARAm model mixes simultaneously GQR effect and ARAm effect.
	"\code{rho}" is the corresponding maintenance effect parameter vector, it is of size 2 : "\code{rho[1]}" characterizes the GQR effect and "\code{rho[2]}" characterizes the ARAm effect.}
	}
}

\subsection{Time to failure distribution of the new unmaintained system}{
	"\code{InitDistr}" characterizes the hazard rate of the new unmaintained system \eqn{h(.)}. The available distributions are:
	\describe{
	\item{\code{Weibull(alpha,beta)}}{for a Weibull initial hazard rate with parametrization \eqn{h(t) = \alpha \beta t^(\beta-1)} with \eqn{\alpha > 0} and \eqn{\beta > 0}.}
	\item{\code{LogLinear(alpha,beta)}}{for a log-linear hazard rate with parametrization \eqn{h(t) = \alpha exp(\beta t)} with \eqn{\alpha > 0}.}
	\item{\code{Weibull3(alpha,beta,c)}}{for a three parameter Weibull initial hazard rate with parametrization \eqn{h(t) = \alpha \beta (t+c)^(\beta-1)} with \eqn{\alpha > 0}, \eqn{\beta > 0} and \eqn{c > 0}.}
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
 The function produces an object of class "\code{sim.vam}" which contains the virtual age model considered.
 }

 \references{Classes of imperfect repair models based on reduction of failure intensity or virtual age:
L. Doyen, O. Gaudoin,
Reliability Engineering and System Safety, Elsevier, 2004, 84 (1), pp.45-56.

Modelling and assessment of aging and efficiency of corrective and planned preventive maintenance:
L. Doyen, O. Gaudoin, IEEE Transactions on Reliability, 2011, 60 (4), pp.759-769.
}

\author{L. Doyen and R. Drouilhet}

\seealso{\code{\link{simulate.sim.vam}} for simulation.

\code{\link{plot.sim.vam}} for plotting characteristics of the model.
}


\examples{
##########
# Simulation model:
# ARA Infinite CM with \rho=0.4 
# Weibull initial intensity h(t)=0.001*2.5*t^(1.5)
simCM<-sim.vam(  ~ (ARAInf(.4) | Weibull(.001,2.5)))
(simulate(simCM,30))
plot(simCM,'i')

##########
# Simulation model:
# Weibull h(t)=0.001*2.5*t^(1.5)
# CM ARA1 (rhoMC=0.6)
# PM ARA Infinite (rhoMP=0.4)
# PM at fixed failure intensity level
simCMPM<-sim.vam(  ~ (ARA1(.9) | Weibull(.001,2.5)) & (ARAInf(.4) | AtIntensity(0.2)))
(simulate(simCMPM,50))
plot(simCMPM,'i')

###########
# Simulation model:
# Weibull h(t)=0.001*2.5*t^(1.5)
# CM ARA infinite (rhoMC=0.3)
# periodic PM randomly chosen upon two types
# PM ARA Infinite (rhoMP=0.6) with probability 0.6
# PM ARA Infinite (rhoMP=-0.2) with probability 0.4
simCMPM_Multi<-sim.vam(  ~ (ARAInf(.3) | Weibull(.001,2.5)) & (ARAInf(.6)+ARAInf(-.2) | Periodic(12,prob=c(0.6,0.4))))
(simulate(simCMPM_Multi,20,nb.system=3))
plot(simCMPM_Multi,'i',system.index=3)

##########
# Combined PM policies:
# Weibull h(t)=0.001*2.6*t^(1.6)
# CM ABAO
# periodic PM randomly chosen upon two types
# PM ARA Infinite (rhoMP=0.6) with probability 0.6
# PM ARA Infinite (rhoMP=-0.2) with probability 0.4
# and AGAN PM renew the system as soon as the failure intensity reaches the level 0.3
simCMPM<-sim.vam(  ~ (ABAO() | Weibull(.001,2.6)) & (ARAInf(.6)+ARAInf(-.2)+AGAN() | Periodic(12,prob=c(0.6,0.4))*AtIntensity(0.3)))
(simulate(simCMPM,15))
plot(simCMPM,'i')

###########
# PM policy not using the same model as for simulation:
# The model used for PM planning
modMPplan <- model.vam( ~ (ARA1(.87) | Weibull(.0015,2.7)) & (ARAInf(.44)))
# The model used for failure time simulation
# a PM is done as soon as the failure intensity of the previous model reaches level 0.3 
simCMPM<-sim.vam( ~ (ARA1(.9) | Weibull(.001,2.5)) & (ARAInf(.35) | AtIntensity(0.3,modMPplan)) )
simData<-simulate(simCMPM,20)
modMPplan_plot<-model.vam(Time & Type ~ (ARA1(.87) | Weibull(.0015,2.7)) & (ARAInf(.44)),data=simData)
# The failure intensity of the model used for PM planning
# a PM is done when this failure intensity reaches the level 0.3
plot(modMPplan_plot,'i',col='blue')
# The failure intensity of the model used for failure simulation
plot(simCMPM,'i-cm-pm',col='darkblue',add=TRUE)
}


\keyword{virtual age}
\keyword{imperfect maintenance}
\keyword{corrective maintenance}
\keyword{preventive maintenance}
\keyword{ARA}
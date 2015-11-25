require(VAM)

#####################
# 1 trajectory
# only CM
#####################

# Simulation model:
# ARA 1 CM with \rho=0.4
# Weibull h(t)=0.001*2.5*t^(1.5)
simARA1<-sim.vam(  ~ (ARA1(.4) | Weibull(.001,2.5)))

n<-30 #number of simulated maintenance
(simData<-simulate(simARA1,n))

# MLE:
# Initialisation of maximization algorithm:
rho_0<-0.5
theta_0<-c(1,1.5)
# Compute the MLE
mleARA1 <- mle.vam(Time & Type ~ (ARA1(rho_0) | Weibull(theta_0[1],theta_0[2])),data=simData)
coef(mleARA1)
# The corresponding plug in estimators can be plotted:
plot(mleARA1,'h',col='darkblue')
points(simData$Time,rep(0,n),type='p',pch='*')
# PROBLEME !!!
# ca ne correspond pas à l'intensité du modèle avec les paramètres estimés
testmodel<-model.vam(Time & Type ~ (ARA1((coef(mleARA1))[3]) | Weibull((coef(mleARA1))[1:2])),data=simData)
plot(testmodel,'h',col='darkblue')
# pour autant ce n'est pas non plus avec les valeurs d'initialisation
testmodel2<-model.vam(Time & Type ~ (ARA1(rho_0) | Weibull(theta_0)),data=simData)
plot(testmodel,'h',col='darkblue')
# or
plot(mleARA1,'H',col='darkblue')
lines(stepfun(simData$Time,0:n))
# MEME PROBLEME POUR H !!!!
plot(testmodel,'H',col='darkblue')
lines(stepfun(simData$Time,0:n))
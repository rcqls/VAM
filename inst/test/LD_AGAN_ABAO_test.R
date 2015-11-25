require(VAM)

#####################
# 1 trajectory
# only CM
#####################

#Simulation model:
# ARA 1 CM with \rho=0.4
# Weibull h(t)=0.001*2.5*t^(1.5)
simARA1<-sim.vam(  ~ (ARA1(.4) | Weibull(.001,2.5)))

n<-10 #number of simulated maintenance
(simData<-simulate(simARA1,n))

#Associate a model to _Time and _Type of simData
modelARA1<-model.vam(Time & Type ~ (ARA1(.4) | Weibull(.001,2.5)),data=simData)
#in order to plot the corresponding intensity
plot(modelARA1,'h',col='darkblue')
points(simData$Time,rep(0,n),type='p',pch='*')
#or cumulative intensity (or virtual age with 'v')
plot(modelARA1,'H',col='darkblue') #problème dans le plot à corriger:saut vers le bas
lines(stepfun(simData$Time,0:n))

#The model associate to simData can be changed
modelAGAN<-model.vam(Time & Type ~ (AGAN() | Weibull(.001,2.5)),data=simData)
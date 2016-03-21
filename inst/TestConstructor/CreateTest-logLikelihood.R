# R script used to develop the testthat test for log-likelihood.R 

# The number of the test
nbtest<-"T1"

switch(nbtest,
  
  T01={
    #Weibull + CM AGAN
    simData<-data.frame(Time=c(3.36,4.04,4.97,5.16),Type=c(-1,-1,-1,-1),row.names=1:4)
    mle <- mle.vam(Time & Type ~ (AGAN() | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.3,0.8,1)
         
    rho<-1
    h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
    H<-function(t) theta[1]*t^(theta[2])
    T<-simData$Time
    Lcalc<-log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rho*T[1]))-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc+log(h(T[3]-rho*T[2]-rho*(1-rho)*T[1]))-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    Lcalc<-Lcalc+log(h(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T02={
    #LogLinear + CM ABAO
    simData<-data.frame(Time=c(3.36,4.04,4.97,5.16),Type=c(-1,-1,-1,-1),row.names=1:4)
    mle <- mle.vam(Time & Type ~ (ABAO() | LogLinear(0.001,2.5)),data=simData)
    theta<-c(0.3,0.8,1)
    
    rho<-0
    h<-function(t) theta[1]*exp(theta[2]*t)
    H<-function(t) theta[1]/theta[2]*(exp(theta[2]*t)-1)
    T<-simData$Time
    Lcalc<-log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rho*T[1]))-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc+log(h(T[3]-rho*T[2]-rho*(1-rho)*T[1]))-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    Lcalc<-Lcalc+log(h(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T03={
    #Weibull3
    simData<-data.frame(Time=c(3.36),Type=c(-1),row.names=1:1)
    mle <- mle.vam(Time & Type ~ (ARAInf(0.4) | Weibull3(0.001,2.5,5)),data=simData)
    theta<-c(0.3,1.8,4,0.6)
    
    c<-theta[3]
    rho<-theta[4]
    h<-function(t) theta[1]*theta[2]*(t+c)^(theta[2]-1)
    H<-function(t) theta[1]*((t+c)^(theta[2])-c^(theta[2]))
    T<-simData$Time
    Lcalc<-log(h(T[1]))-H(T[1])
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T1={
    #Weibull
    simData<-data.frame(Time=c(3.36),Type=c(-1),row.names=1:1)
    mle <- mle.vam(Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.3,1.8,0.6)
    
    rho<-theta[3]
    h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
    H<-function(t) theta[1]*t^(theta[2])
    T<-simData$Time
    Lcalc<-log(h(T[1]))-H(T[1])
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T2={
    #LogLinear
    simData<-data.frame(Time=c(3.36),Type=c(-1),row.names=1:1)
    mle <- mle.vam(Time & Type ~ (ARAInf(0.4) | LogLinear(0.001,2.5)),data=simData)
    theta<-c(0.3,0.8,0.6)
    
    rho<-theta[3]
    h<-function(t) theta[1]*exp(theta[2]*t)
    H<-function(t) theta[1]/theta[2]*(exp(theta[2]*t)-1)
    T<-simData$Time
    Lcalc<-log(h(T[1]))-H(T[1])
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T3={
    #Weibull + CM ARAInf
    simData<-data.frame(Time=c(3.36,4.04,4.97,5.16),Type=c(-1,-1,-1,-1),row.names=1:4)
    mle <- mle.vam(Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.3,0.8,0.6)
    
    rho<-theta[3]
    h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
    H<-function(t) theta[1]*t^(theta[2])
    T<-simData$Time
    Lcalc<-log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rho*T[1]))-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc+log(h(T[3]-rho*T[2]-rho*(1-rho)*T[1]))-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    Lcalc<-Lcalc+log(h(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T4={
    #Weibull + CM ARA1
    simData<-data.frame(Time=c(3.36,4.04,4.97,5.16),Type=c(-1,-1,-1,-1),row.names=1:4)
    mle <- mle.vam(Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.3,0.8,0.6)
    
    rho<-theta[3]
    h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
    H<-function(t) theta[1]*t^(theta[2])
    T<-simData$Time
    Lcalc<-log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rho*T[1]))-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc+log(h(T[3]-rho*T[2]))-(H(T[3]-rho*T[2])-H(T[2]-rho*T[2]))
    Lcalc<-Lcalc+log(h(T[4]-rho*T[3]))-(H(T[4]-rho*T[3])-H(T[3]-rho*T[3]))
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T5={
    #LogLinear + CM ARAInf
    simData<-data.frame(Time=c(3.36,4.04,4.97,5.16),Type=c(-1,-1,-1,-1),row.names=1:4)
    mle <- mle.vam(Time & Type ~ (ARAInf(0.4) | LogLinear(0.001,2.5)),data=simData)
    theta<-c(0.3,0.8,0.6)
    
    rho<-theta[3]
    h<-function(t) theta[1]*exp(theta[2]*t)
    H<-function(t) theta[1]/theta[2]*(exp(theta[2]*t)-1)
    T<-simData$Time
    Lcalc<-log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rho*T[1]))-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc+log(h(T[3]-rho*T[2]-rho*(1-rho)*T[1]))-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    Lcalc<-Lcalc+log(h(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T6={
    #Weibull + CM ARA1 + Censorship
    simData<-data.frame(Time=c(3.36,4.04,4.97,5.16),Type=c(-1,-1,-1,0),row.names=1:4)
    mle <- mle.vam(Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.3,0.8,0.6)
    
    rho<-theta[3]
    h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
    H<-function(t) theta[1]*t^(theta[2])
    T<-simData$Time
    Lcalc<-log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rho*T[1]))-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc+log(h(T[3]-rho*T[2]))-(H(T[3]-rho*T[2])-H(T[2]-rho*T[2]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3])-H(T[3]-rho*T[3]))
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T7={
    #Weibull + multiSystems
    simData<-data.frame(System=c(1,2),Time=c(3.36,2.34),Type=c(-1,-1),row.names=1:2)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.3,1.8,0.6)
    
    rho<-theta[3]
    h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
    H<-function(t) theta[1]*t^(theta[2])
    T<-simData$Time[simData$System==1]
    Lcalc<-log(h(T[1]))-H(T[1])
    T<-simData$Time[simData$System==2]
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T8={
    #Weibull + CM ARAInf + Censorship + multiSystems
    simData<-data.frame(System=c(1,1,1,1,2,2,2,3),Time=c(3.36,4.04,4.97,5.16,2.34,3.46,5.02,4),Type=c(-1,-1,-1,0,-1,-1,-1,0),row.names=1:8)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.3,1.8,0.6)
    
    rho<-theta[3]
    h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
    H<-function(t) theta[1]*t^(theta[2])
    T<-simData$Time[simData$System==1]
    Lcalc<-log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rho*T[1]))-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc+log(h(T[3]-rho*T[2]-rho*(1-rho)*T[1]))-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
    T<-simData$Time[simData$System==2]
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rho*T[1]))-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc+log(h(T[3]-rho*T[2]-rho*(1-rho)*T[1]))-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    T<-simData$Time[simData$System==3]
    Lcalc<-Lcalc-H(T[1])
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T9={
    #Weibull + CM ARAInf + PM ARAInf + mutlisystems
    simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (ARAInf(0.5)),data=simData)
    theta<-c(0.3,1.8,0.3,0.8)
    
    rhoMC<-theta[3]
    rhoMP<-theta[4]
    h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
    H<-function(t) theta[1]*t^(theta[2])
    T<-simData$Time[simData$System==1]
    rho<-rhoMP
    Lcalc<--H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
    T<-simData$Time[simData$System==2]
    rho<-rhoMC
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rho*T[1]))-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc+log(h(T[3]-rho*T[2]-rho*(1-rho)*T[1]))-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
    T<-simData$Time[simData$System==3]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP*T[1]))-(H(T[2]-rhoMP*T[1])-H(T[1]-rhoMP*T[1]))
    V<-(1-rhoMC)*(T[2]-rhoMP*T[1])
    Lcalc<-Lcalc+log(h(T[3]-T[2]+V))-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMC)*(T[3]-T[2]+V)
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==4]
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rhoMC*T[1])-H(T[1]-rhoMC*T[1]))
    V<-(1-rhoMP)*(T[2]-rhoMC*T[1])
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP)*(T[3]-T[2]+V)
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==5]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP*T[1]))-(H(T[2]-rhoMP*T[1])-H(T[1]-rhoMP*T[1]))
    V<-(1-rhoMC)*(T[2]-rhoMP*T[1])
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP)*(T[3]-T[2]+V)
    Lcalc<-Lcalc+log(h(T[4]-T[3]+V))-(H(T[4]-T[3]+V)-H(V))
    V<-(1-rhoMC)*(T[4]-T[3]+V)
    Lcalc<-Lcalc+log(h(T[5]-T[4]+V))-(H(T[5]-T[4]+V)-H(V))
    V<-(1-rhoMC)*(T[5]-T[4]+V)
    Lcalc<-Lcalc-(H(T[6]-T[5]+V)-H(V))
    V<-(1-rhoMP)*(T[6]-T[5]+V)
    Lcalc<-Lcalc-(H(T[7]-T[6]+V)-H(V))
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T91={
    #Weibull + CM ARAInf + PM AGAN + mutlisystems
    simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (AGAN()),data=simData)
    theta<-c(0.3,1.8,0.3)
    
    rhoMC<-theta[3]
    rhoMP<-1
    h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
    H<-function(t) theta[1]*t^(theta[2])
    T<-simData$Time[simData$System==1]
    rho<-rhoMP
    Lcalc<--H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
    T<-simData$Time[simData$System==2]
    rho<-rhoMC
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rho*T[1]))-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc+log(h(T[3]-rho*T[2]-rho*(1-rho)*T[1]))-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
    T<-simData$Time[simData$System==3]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP*T[1]))-(H(T[2]-rhoMP*T[1])-H(T[1]-rhoMP*T[1]))
    V<-(1-rhoMC)*(T[2]-rhoMP*T[1])
    Lcalc<-Lcalc+log(h(T[3]-T[2]+V))-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMC)*(T[3]-T[2]+V)
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==4]
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rhoMC*T[1])-H(T[1]-rhoMC*T[1]))
    V<-(1-rhoMP)*(T[2]-rhoMC*T[1])
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP)*(T[3]-T[2]+V)
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==5]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP*T[1]))-(H(T[2]-rhoMP*T[1])-H(T[1]-rhoMP*T[1]))
    V<-(1-rhoMC)*(T[2]-rhoMP*T[1])
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP)*(T[3]-T[2]+V)
    Lcalc<-Lcalc+log(h(T[4]-T[3]+V))-(H(T[4]-T[3]+V)-H(V))
    V<-(1-rhoMC)*(T[4]-T[3]+V)
    Lcalc<-Lcalc+log(h(T[5]-T[4]+V))-(H(T[5]-T[4]+V)-H(V))
    V<-(1-rhoMC)*(T[5]-T[4]+V)
    Lcalc<-Lcalc-(H(T[6]-T[5]+V)-H(V))
    V<-(1-rhoMP)*(T[6]-T[5]+V)
    Lcalc<-Lcalc-(H(T[7]-T[6]+V)-H(V))
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T92={
    #Weibull + CM ARAInf + PM ABAO + mutlisystems
    simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (ABAO()),data=simData)
    theta<-c(0.3,1.8,0.3)
    
    rhoMC<-theta[3]
    rhoMP<-0
    h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
    H<-function(t) theta[1]*t^(theta[2])
    T<-simData$Time[simData$System==1]
    rho<-rhoMP
    Lcalc<--H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
    T<-simData$Time[simData$System==2]
    rho<-rhoMC
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rho*T[1]))-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc+log(h(T[3]-rho*T[2]-rho*(1-rho)*T[1]))-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
    T<-simData$Time[simData$System==3]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP*T[1]))-(H(T[2]-rhoMP*T[1])-H(T[1]-rhoMP*T[1]))
    V<-(1-rhoMC)*(T[2]-rhoMP*T[1])
    Lcalc<-Lcalc+log(h(T[3]-T[2]+V))-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMC)*(T[3]-T[2]+V)
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==4]
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rhoMC*T[1])-H(T[1]-rhoMC*T[1]))
    V<-(1-rhoMP)*(T[2]-rhoMC*T[1])
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP)*(T[3]-T[2]+V)
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==5]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP*T[1]))-(H(T[2]-rhoMP*T[1])-H(T[1]-rhoMP*T[1]))
    V<-(1-rhoMC)*(T[2]-rhoMP*T[1])
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP)*(T[3]-T[2]+V)
    Lcalc<-Lcalc+log(h(T[4]-T[3]+V))-(H(T[4]-T[3]+V)-H(V))
    V<-(1-rhoMC)*(T[4]-T[3]+V)
    Lcalc<-Lcalc+log(h(T[5]-T[4]+V))-(H(T[5]-T[4]+V)-H(V))
    V<-(1-rhoMC)*(T[5]-T[4]+V)
    Lcalc<-Lcalc-(H(T[6]-T[5]+V)-H(V))
    V<-(1-rhoMP)*(T[6]-T[5]+V)
    Lcalc<-Lcalc-(H(T[7]-T[6]+V)-H(V))
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T93={
    #Weibull + CM AGAN + PM ARAInf + mutlisystems
    simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (AGAN() | Weibull(0.001,2.5)) & (ARAInf(0.5)),data=simData)
    theta<-c(0.3,1.8,0.6)
    
    rhoMC<-1
    rhoMP<-theta[3]
    h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
    H<-function(t) theta[1]*t^(theta[2])
    T<-simData$Time[simData$System==1]
    rho<-rhoMP
    Lcalc<--H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
    T<-simData$Time[simData$System==2]
    rho<-rhoMC
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rho*T[1]))-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc+log(h(T[3]-rho*T[2]-rho*(1-rho)*T[1]))-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
    T<-simData$Time[simData$System==3]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP*T[1]))-(H(T[2]-rhoMP*T[1])-H(T[1]-rhoMP*T[1]))
    V<-(1-rhoMC)*(T[2]-rhoMP*T[1])
    Lcalc<-Lcalc+log(h(T[3]-T[2]+V))-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMC)*(T[3]-T[2]+V)
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==4]
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rhoMC*T[1])-H(T[1]-rhoMC*T[1]))
    V<-(1-rhoMP)*(T[2]-rhoMC*T[1])
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP)*(T[3]-T[2]+V)
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==5]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP*T[1]))-(H(T[2]-rhoMP*T[1])-H(T[1]-rhoMP*T[1]))
    V<-(1-rhoMC)*(T[2]-rhoMP*T[1])
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP)*(T[3]-T[2]+V)
    Lcalc<-Lcalc+log(h(T[4]-T[3]+V))-(H(T[4]-T[3]+V)-H(V))
    V<-(1-rhoMC)*(T[4]-T[3]+V)
    Lcalc<-Lcalc+log(h(T[5]-T[4]+V))-(H(T[5]-T[4]+V)-H(V))
    V<-(1-rhoMC)*(T[5]-T[4]+V)
    Lcalc<-Lcalc-(H(T[6]-T[5]+V)-H(V))
    V<-(1-rhoMP)*(T[6]-T[5]+V)
    Lcalc<-Lcalc-(H(T[7]-T[6]+V)-H(V))
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T94={
    #Weibull + CM ABAO + PM ARAInf + mutlisystems
    simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ABAO() | Weibull(0.001,2.5)) & (ARAInf(0.5)),data=simData)
    theta<-c(0.3,1.8,0.6)
    
    rhoMC<-0
    rhoMP<-theta[3]
    h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
    H<-function(t) theta[1]*t^(theta[2])
    T<-simData$Time[simData$System==1]
    rho<-rhoMP
    Lcalc<--H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
    T<-simData$Time[simData$System==2]
    rho<-rhoMC
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rho*T[1]))-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc+log(h(T[3]-rho*T[2]-rho*(1-rho)*T[1]))-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
    T<-simData$Time[simData$System==3]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP*T[1]))-(H(T[2]-rhoMP*T[1])-H(T[1]-rhoMP*T[1]))
    V<-(1-rhoMC)*(T[2]-rhoMP*T[1])
    Lcalc<-Lcalc+log(h(T[3]-T[2]+V))-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMC)*(T[3]-T[2]+V)
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==4]
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rhoMC*T[1])-H(T[1]-rhoMC*T[1]))
    V<-(1-rhoMP)*(T[2]-rhoMC*T[1])
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP)*(T[3]-T[2]+V)
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==5]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP*T[1]))-(H(T[2]-rhoMP*T[1])-H(T[1]-rhoMP*T[1]))
    V<-(1-rhoMC)*(T[2]-rhoMP*T[1])
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP)*(T[3]-T[2]+V)
    Lcalc<-Lcalc+log(h(T[4]-T[3]+V))-(H(T[4]-T[3]+V)-H(V))
    V<-(1-rhoMC)*(T[4]-T[3]+V)
    Lcalc<-Lcalc+log(h(T[5]-T[4]+V))-(H(T[5]-T[4]+V)-H(V))
    V<-(1-rhoMC)*(T[5]-T[4]+V)
    Lcalc<-Lcalc-(H(T[6]-T[5]+V)-H(V))
    V<-(1-rhoMP)*(T[6]-T[5]+V)
    Lcalc<-Lcalc-(H(T[7]-T[6]+V)-H(V))
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T10={
    #LogLinear + CM ARA1 + PM AR1 + mutlisystems
    simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5)) & (ARA1(0.5)),data=simData)
    theta<-c(0.3,1.8,0.3,0.8)
    
    rhoMC<-theta[3]
    rhoMP<-theta[4]
    h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
    H<-function(t) theta[1]*t^(theta[2])
    T<-simData$Time[simData$System==1]
    rho<-rhoMP
    Lcalc<--H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc-(H(T[3]-rho*T[2])-H(T[2]-rho*T[2]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3])-H(T[3]-rho*T[3]))
    T<-simData$Time[simData$System==2]
    rho<-rhoMC
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rho*T[1]))-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc+log(h(T[3]-rho*T[2]))-(H(T[3]-rho*T[2])-H(T[2]-rho*T[2]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3])-H(T[3]-rho*T[3]))
    T<-simData$Time[simData$System==3]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP*T[1]))-(H(T[2]-rhoMP*T[1])-H(T[1]-rhoMP*T[1]))
    V<-(1-rhoMC)*(T[2]-T[1])+T[1]-rhoMP*T[1]
    Lcalc<-Lcalc+log(h(T[3]-T[2]+V))-(H(T[3]-T[2]+V)-H(V))
    V<-V+(1-rhoMC)*(T[3]-T[2])
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==4]
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rhoMC*T[1])-H(T[1]-rhoMC*T[1]))
    V<-(1-rhoMP)*(T[2]-T[1])+T[1]-rhoMC*T[1]
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP)*(T[3]-T[2])+V
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==5]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP*T[1]))-(H(T[2]-rhoMP*T[1])-H(T[1]-rhoMP*T[1]))
    V<-(1-rhoMC)*(T[2]-T[1])+T[1]-rhoMP*T[1]
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP)*(T[3]-T[2])+V
    Lcalc<-Lcalc+log(h(T[4]-T[3]+V))-(H(T[4]-T[3]+V)-H(V))
    V<-(1-rhoMC)*(T[4]-T[3])+V
    Lcalc<-Lcalc+log(h(T[5]-T[4]+V))-(H(T[5]-T[4]+V)-H(V))
    V<-(1-rhoMC)*(T[5]-T[4])+V
    Lcalc<-Lcalc-(H(T[6]-T[5]+V)-H(V))
    V<-(1-rhoMP)*(T[6]-T[5])+V
    Lcalc<-Lcalc-(H(T[7]-T[6]+V)-H(V))
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T11={
    #Weibull + CM ARA1 + PM ARAInf + mutlisystems
    simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5)) & (ARAInf(0.5)),data=simData)
    theta<-c(0.3,1.8,0.3,0.8)
    #n<-11
    #Data<-data.frame(System=simData$System[1:n],Time=simData$Time[1:n],Type=simData$Type[1:n])
    #mle2 <- mle.vam(System & Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5)) & (ARAInf(0.5)),data=Data)
    #L2<-logLikelihood(mle2,theta,c(TRUE,FALSE,FALSE))
    rhoMC<-theta[3]
    rhoMP<-theta[4]
    h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
    H<-function(t) theta[1]*t^(theta[2])
    T<-simData$Time[simData$System==1]
    rho<-rhoMP
    Lcalc<--H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
    T<-simData$Time[simData$System==2]
    rho<-rhoMC
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rho*T[1]))-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc+log(h(T[3]-rho*T[2]))-(H(T[3]-rho*T[2])-H(T[2]-rho*T[2]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3])-H(T[3]-rho*T[3]))
    T<-simData$Time[simData$System==3]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP*T[1]))-(H(T[2]-rhoMP*T[1])-H(T[1]-rhoMP*T[1]))
    V<-(1-rhoMC)*(T[2]-T[1])+T[1]-rhoMP*T[1]
    Lcalc<-Lcalc+log(h(T[3]-T[2]+V))-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMC)*(T[3]-T[2])+V
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==4]
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rhoMC*T[1])-H(T[1]-rhoMC*T[1]))
    V<-(1-rhoMP)*(T[2]-rhoMC*T[1])
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP)*(T[3]-T[2]+V)
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==5]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP*T[1]))-(H(T[2]-rhoMP*T[1])-H(T[1]-rhoMP*T[1]))
    V<-(1-rhoMC)*(T[2]-T[1])+T[1]-rhoMP*T[1]
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP)*(T[3]-T[2]+V)
    Lcalc<-Lcalc+log(h(T[4]-T[3]+V))-(H(T[4]-T[3]+V)-H(V))
    V<-(1-rhoMC)*(T[4]-T[3])+V
    Lcalc<-Lcalc+log(h(T[5]-T[4]+V))-(H(T[5]-T[4]+V)-H(V))
    V<-(1-rhoMC)*(T[5]-T[4])+V
    Lcalc<-Lcalc-(H(T[6]-T[5]+V)-H(V))
    V<-(1-rhoMP)*(T[6]-T[5]+V)
    Lcalc<-Lcalc-(H(T[7]-T[6]+V)-H(V))
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T12={
    #Weibull + CM ARAInf + PM ARA1 + mutlisystems
    simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (ARA1(0.5)),data=simData)
    theta<-c(0.3,1.8,0.3,0.8)
    
    rhoMC<-theta[3]
    rhoMP<-theta[4]
    h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
    H<-function(t) theta[1]*t^(theta[2])
    T<-simData$Time[simData$System==1]
    rho<-rhoMP
    Lcalc<--H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc-(H(T[3]-rho*T[2])-H(T[2]-rho*T[2]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3])-H(T[3]-rho*T[3]))
    T<-simData$Time[simData$System==2]
    rho<-rhoMC
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rho*T[1]))-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc+log(h(T[3]-rho*T[2]-rho*(1-rho)*T[1]))-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
    T<-simData$Time[simData$System==3]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP*T[1]))-(H(T[2]-rhoMP*T[1])-H(T[1]-rhoMP*T[1]))
    V<-(1-rhoMC)*(T[2]-rhoMP*T[1])
    Lcalc<-Lcalc+log(h(T[3]-T[2]+V))-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMC)*(T[3]-T[2]+V)
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==4]
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rhoMC*T[1])-H(T[1]-rhoMC*T[1]))
    V<-(1-rhoMP)*(T[2]-T[1])+T[1]-rhoMC*T[1]
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP)*(T[3]-T[2])+V
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==5]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP*T[1]))-(H(T[2]-rhoMP*T[1])-H(T[1]-rhoMP*T[1]))
    V<-(1-rhoMC)*(T[2]-rhoMP*T[1])
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP)*(T[3]-T[2])+V
    Lcalc<-Lcalc+log(h(T[4]-T[3]+V))-(H(T[4]-T[3]+V)-H(V))
    V<-(1-rhoMC)*(T[4]-T[3]+V)
    Lcalc<-Lcalc+log(h(T[5]-T[4]+V))-(H(T[5]-T[4]+V)-H(V))
    V<-(1-rhoMC)*(T[5]-T[4]+V)
    Lcalc<-Lcalc-(H(T[6]-T[5]+V)-H(V))
    V<-(1-rhoMP)*(T[6]-T[5])+V
    Lcalc<-Lcalc-(H(T[7]-T[6]+V)-H(V))
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T13={
    #Weibull + CM ARAInf + PM ARAInf + PM ARA1 + mutlisystems
    simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,10)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43),Type=c(1,2,2,1, -1,-1,-1,0, 1,-1,-1,2, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,0),row.names=1:26)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (ARAInf(0.5)+ARA1(-1)),data=simData)
    theta<-c(0.3,1.8,0.3,0.8,-1)
    
    rhoMC<-theta[3]
    rhoMP1<-theta[4]
    rhoMP2<-theta[5]
    h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
    H<-function(t) theta[1]*t^(theta[2])
    T<-simData$Time[simData$System==1]
    Lcalc<--H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rhoMP1*T[1])-H(T[1]-rhoMP1*T[1]))
    V<-(1-rhoMP2)*(T[2]-T[1])+T[1]-rhoMP1*T[1]
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP2)*(T[3]-T[2])+V
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==2]
    rho<-rhoMC
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rho*T[1]))-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc+log(h(T[3]-rho*T[2]-rho*(1-rho)*T[1]))-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
    T<-simData$Time[simData$System==3]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP1*T[1]))-(H(T[2]-rhoMP1*T[1])-H(T[1]-rhoMP1*T[1]))
    V<-(1-rhoMC)*(T[2]-rhoMP1*T[1])
    Lcalc<-Lcalc+log(h(T[3]-T[2]+V))-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMC)*(T[3]-T[2]+V)
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==4]
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rhoMC*T[1])-H(T[1]-rhoMC*T[1]))
    V<-(1-rhoMP2)*(T[2]-T[1])+T[1]-rhoMC*T[1]
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP2)*(T[3]-T[2])+V
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==5]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP1*T[1]))-(H(T[2]-rhoMP1*T[1])-H(T[1]-rhoMP1*T[1]))
    V<-(1-rhoMC)*(T[2]-rhoMP1*T[1])
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP1)*(T[3]-T[2]+V)
    Lcalc<-Lcalc+log(h(T[4]-T[3]+V))-(H(T[4]-T[3]+V)-H(V))
    V<-(1-rhoMC)*(T[4]-T[3]+V)
    Lcalc<-Lcalc+log(h(T[5]-T[4]+V))-(H(T[5]-T[4]+V)-H(V))
    V<-(1-rhoMC)*(T[5]-T[4]+V)
    Lcalc<-Lcalc-(H(T[6]-T[5]+V)-H(V))
    V<-(1-rhoMP2)*(T[6]-T[5])+V
    Lcalc<-Lcalc-(H(T[7]-T[6]+V)-H(V))
    V<-(1-rhoMP1)*(T[7]-T[6]+V)
    Lcalc<-Lcalc-(H(T[8]-T[7]+V)-H(V))
    V<-(1-rhoMP2)*(T[8]-T[7])+V
    Lcalc<-Lcalc+log(h(T[9]-T[8]+V))-(H(T[9]-T[8]+V)-H(V))
    V<-(1-rhoMC)*(T[9]-T[8]+V)
    Lcalc<-Lcalc-(H(T[10]-T[9]+V)-H(V))
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T14={
    #Weibull + CM ARAInf + PM AGAN + PM ARA1 + mutlisystems
    simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,10)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43),Type=c(1,2,2,1, -1,-1,-1,0, 1,-1,-1,2, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,0),row.names=1:26)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (AGAN()+ARA1(-1)),data=simData)
    theta<-c(0.3,1.8,0.3,-1)
    
    rhoMC<-theta[3]
    rhoMP1<-1
    rhoMP2<-theta[4]
    h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
    H<-function(t) theta[1]*t^(theta[2])
    T<-simData$Time[simData$System==1]
    Lcalc<--H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rhoMP1*T[1])-H(T[1]-rhoMP1*T[1]))
    V<-(1-rhoMP2)*(T[2]-T[1])+T[1]-rhoMP1*T[1]
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP2)*(T[3]-T[2])+V
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==2]
    rho<-rhoMC
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rho*T[1]))-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc+log(h(T[3]-rho*T[2]-rho*(1-rho)*T[1]))-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
    T<-simData$Time[simData$System==3]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP1*T[1]))-(H(T[2]-rhoMP1*T[1])-H(T[1]-rhoMP1*T[1]))
    V<-(1-rhoMC)*(T[2]-rhoMP1*T[1])
    Lcalc<-Lcalc+log(h(T[3]-T[2]+V))-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMC)*(T[3]-T[2]+V)
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==4]
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rhoMC*T[1])-H(T[1]-rhoMC*T[1]))
    V<-(1-rhoMP2)*(T[2]-T[1])+T[1]-rhoMC*T[1]
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP2)*(T[3]-T[2])+V
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==5]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP1*T[1]))-(H(T[2]-rhoMP1*T[1])-H(T[1]-rhoMP1*T[1]))
    V<-(1-rhoMC)*(T[2]-rhoMP1*T[1])
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP1)*(T[3]-T[2]+V)
    Lcalc<-Lcalc+log(h(T[4]-T[3]+V))-(H(T[4]-T[3]+V)-H(V))
    V<-(1-rhoMC)*(T[4]-T[3]+V)
    Lcalc<-Lcalc+log(h(T[5]-T[4]+V))-(H(T[5]-T[4]+V)-H(V))
    V<-(1-rhoMC)*(T[5]-T[4]+V)
    Lcalc<-Lcalc-(H(T[6]-T[5]+V)-H(V))
    V<-(1-rhoMP2)*(T[6]-T[5])+V
    Lcalc<-Lcalc-(H(T[7]-T[6]+V)-H(V))
    V<-(1-rhoMP1)*(T[7]-T[6]+V)
    Lcalc<-Lcalc-(H(T[8]-T[7]+V)-H(V))
    V<-(1-rhoMP2)*(T[8]-T[7])+V
    Lcalc<-Lcalc+log(h(T[9]-T[8]+V))-(H(T[9]-T[8]+V)-H(V))
    V<-(1-rhoMC)*(T[9]-T[8]+V)
    Lcalc<-Lcalc-(H(T[10]-T[9]+V)-H(V))
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  T15={
    #Weibull3 + CM ARAInf + PM AGAN + PM ARA1 + mutlisystems
    simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,10)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43),Type=c(1,2,2,1, -1,-1,-1,0, 1,-1,-1,2, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,0),row.names=1:26)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull3(0.001,2.5,3)) & (AGAN()+ARA1(-1)),data=simData)
    theta<-c(0.3,1.8,4,0.3,-1)
    
    c<-theta[3]
    rhoMC<-theta[4]
    rhoMP1<-1
    rhoMP2<-theta[5]
    h<-function(t) theta[1]*theta[2]*(t+c)^(theta[2]-1)
    H<-function(t) theta[1]*((t+c)^(theta[2])-(c)^(theta[2]))
    T<-simData$Time[simData$System==1]
    Lcalc<--H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rhoMP1*T[1])-H(T[1]-rhoMP1*T[1]))
    V<-(1-rhoMP2)*(T[2]-T[1])+T[1]-rhoMP1*T[1]
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP2)*(T[3]-T[2])+V
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==2]
    rho<-rhoMC
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rho*T[1]))-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
    Lcalc<-Lcalc+log(h(T[3]-rho*T[2]-rho*(1-rho)*T[1]))-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
    Lcalc<-Lcalc-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
    T<-simData$Time[simData$System==3]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP1*T[1]))-(H(T[2]-rhoMP1*T[1])-H(T[1]-rhoMP1*T[1]))
    V<-(1-rhoMC)*(T[2]-rhoMP1*T[1])
    Lcalc<-Lcalc+log(h(T[3]-T[2]+V))-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMC)*(T[3]-T[2]+V)
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==4]
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc-(H(T[2]-rhoMC*T[1])-H(T[1]-rhoMC*T[1]))
    V<-(1-rhoMP2)*(T[2]-T[1])+T[1]-rhoMC*T[1]
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP2)*(T[3]-T[2])+V
    Lcalc<-Lcalc-(H(T[4]-T[3]+V)-H(V))
    T<-simData$Time[simData$System==5]
    Lcalc<-Lcalc-H(T[1])
    Lcalc<-Lcalc+log(h(T[2]-rhoMP1*T[1]))-(H(T[2]-rhoMP1*T[1])-H(T[1]-rhoMP1*T[1]))
    V<-(1-rhoMC)*(T[2]-rhoMP1*T[1])
    Lcalc<-Lcalc-(H(T[3]-T[2]+V)-H(V))
    V<-(1-rhoMP1)*(T[3]-T[2]+V)
    Lcalc<-Lcalc+log(h(T[4]-T[3]+V))-(H(T[4]-T[3]+V)-H(V))
    V<-(1-rhoMC)*(T[4]-T[3]+V)
    Lcalc<-Lcalc+log(h(T[5]-T[4]+V))-(H(T[5]-T[4]+V)-H(V))
    V<-(1-rhoMC)*(T[5]-T[4]+V)
    Lcalc<-Lcalc-(H(T[6]-T[5]+V)-H(V))
    V<-(1-rhoMP2)*(T[6]-T[5])+V
    Lcalc<-Lcalc-(H(T[7]-T[6]+V)-H(V))
    V<-(1-rhoMP1)*(T[7]-T[6]+V)
    Lcalc<-Lcalc-(H(T[8]-T[7]+V)-H(V))
    V<-(1-rhoMP2)*(T[8]-T[7])+V
    Lcalc<-Lcalc+log(h(T[9]-T[8]+V))-(H(T[9]-T[8]+V)-H(V))
    V<-(1-rhoMC)*(T[9]-T[8]+V)
    Lcalc<-Lcalc-(H(T[10]-T[9]+V)-H(V))
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
)

L<-logLikelihood(mle,theta,c(TRUE,FALSE,FALSE))
dL<-logLikelihood(mle,theta,c(FALSE,TRUE,FALSE))
d2L<-logLikelihood(mle,theta,c(FALSE,FALSE,TRUE))
C<-contrast(mle,theta,c(TRUE,FALSE,FALSE))
dC<-contrast(mle,theta,c(FALSE,TRUE,FALSE))
d2C<-contrast(mle,theta,c(FALSE,FALSE,TRUE))

epsilon<-0.000001
EstdL<-dL
EstdC<-dC
for(i in (1:length(dL))){
  theta1<-theta
  theta1[i]<-theta1[i]+epsilon
  EstdL[i]<-(logLikelihood(mle,theta1,c(TRUE,FALSE,FALSE))-logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)))/epsilon
}
for(i in (1:length(dC))){
  theta1<-theta
  theta1[i+1]<-theta1[i+1]+epsilon
  EstdC[i]<-(contrast(mle,theta1,c(TRUE,FALSE,FALSE))-contrast(mle,theta,c(TRUE,FALSE,FALSE)))/epsilon
}
Estd2L<-d2L
Estd2C<-d2C
for(i in (1:(dim(d2L)[1]))){
  theta1<-theta
  theta1[i]<-theta1[i]+epsilon
  Estd2L[i,]<-(logLikelihood(mle,theta1,c(FALSE,TRUE,FALSE))-logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)))/epsilon
}
if(!(is.null(dim(d2C)))){
  for(i in (1:(dim(d2C)[1]))){
    theta1<-theta
    theta1[i+1]<-theta1[i+1]+epsilon
    Estd2C[i,]<-(contrast(mle,theta1,c(FALSE,TRUE,FALSE))-contrast(mle,theta,c(FALSE,TRUE,FALSE)))/epsilon
  }
} else {
  theta1<-theta
  theta1[2]<-theta1[2]+epsilon
  Estd2C<-(contrast(mle,theta1,c(FALSE,TRUE,FALSE))-contrast(mle,theta,c(FALSE,TRUE,FALSE)))/epsilon
}
print(L)
print(Lcalc)
print(C)
print(Ccalc)
print(dL)
print(EstdL)
print(dC)
print(EstdC)
print(d2L)
print(Estd2L)
print(d2C)
print(Estd2C)

ecrit<-function(aprinter,nl,nc){
  res<-""
  if((nl>0)&&(nc>0)){
  for(i in 1:nc){
    for(j in 1:nl){
      res<-paste(res,aprinter[j,i],sep=',')
    }
  }}
  else if((nl>0)||(nc>0)){
    for(j in 1:max(nl,nc)){
      res<-paste(res,aprinter[j],sep=',')
    }
  }
  else res<-paste(res,aprinter,sep=',')
  print(res)
  
}

ecrit(L,0,0)
ecrit(dL,0,length(dL))
ecrit(d2L,dim(d2L)[1],dim(d2L)[2])

ecrit(C,0,0)
ecrit(dC,0,length(dC))
if(!(is.null(dim(d2C)))){
  ecrit(d2C,dim(d2C)[1],dim(d2C)[2])
} else {
  ecrit(d2C,0,0)
}
# R script used to develop the testthat test for log-likelihood.R 

# The number of the test
nbtest<-"TQR4"

switch(nbtest,
       TGQRARA_5={
         #Weibull + MC ARA1 + MP GQ_RARAInf-log + MultiSystems + Censorship
         simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
         #simData<-simData[1:12,]
         mle <- mle.vam(System & Time & Type ~ (ARA1(0.2) | Weibull(0.001,2.5)) & (GQR_ARAInf(0.7,-1.3|log)),data=simData)
         theta<-c(0.3,2.2,0.7,0.4,-0.9)
         #L<-logLikelihood(mle,theta,c(TRUE,FALSE,FALSE))
         
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         rhoMC<-theta[3]
         rhoMP<-theta[4]
         rhoMP2<-theta[5]
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         Lcalc<-0
         T<-c(0,simData$Time[simData$System==1])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i])+V); A<-rhoMP^(log(2)); C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i])+V); A<-rhoMP^(log(3)); C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i])+V); A<-rhoMP^(log(4)); C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==2])
         V<-0; C<-0; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==3])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i])+V); A<-rhoMP^(log(2)); C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==4])
         V<-0; C<-0; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i])+V); A<-rhoMP^(log(2));  C<-1; i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i])+V); A<-rhoMP^(log(3)); C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==5])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i])+V); A<-rhoMP^(log(2)); C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i])+V); A<-rhoMP^(log(3)); C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i])+V); A<-rhoMP^(log(4)); C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         fix<-rep(TRUE,length(theta))
         fix[1]=FALSE
         alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
         Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
       },
       TGQRARA_4={
         #Weibull + MC AGAN + MP GQ_RARA1-log + MultiSystems + Censorship
         simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
         #simData<-simData[1:12,]
         mle <- mle.vam(System & Time & Type ~ (AGAN() | Weibull(0.001,2.5)) & (GQR_ARA1(0.7,-1.3|log)),data=simData)
         theta<-c(0.3,2.2,0.4,-0.9)
         #L<-logLikelihood(mle,theta,c(TRUE,FALSE,FALSE))
         
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         rhoMP<-theta[3]
         rhoMP2<-theta[4]
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         Lcalc<-0
         T<-c(0,simData$Time[simData$System==1])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i]))+V; A<-rhoMP^(log(2)); C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i]))+V; A<-rhoMP^(log(3)); C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i]))+V; A<-rhoMP^(log(4)); C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==2])
         V<-0; C<-0; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-0; A<-1; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-0; A<-1; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-0; A<-1; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==3])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i]))+V; A<-rhoMP^(log(2)); C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-0; A<-1; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-0; A<-1; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==4])
         V<-0; C<-0; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-0; A<-1; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i]))+V; A<-rhoMP^(log(2));  C<-1; i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i]))+V; A<-rhoMP^(log(3)); C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==5])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i]))+V; A<-rhoMP^(log(2)); C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-0; A<-1; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i]))+V; A<-rhoMP^(log(3)); C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-0; A<-1; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-0; A<-1; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i]))+V; A<-rhoMP^(log(4)); C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         fix<-rep(TRUE,length(theta))
         fix[1]=FALSE
         alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
         Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
       },
       TGQRARA_3={
         #Weibull + MC ARA1 + MP GQ_RARA1-log + MultiSystems + Censorship
         simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
         #simData<-simData[1:12,]
         mle <- mle.vam(System & Time & Type ~ (ARA1(0.2) | Weibull(0.001,2.5)) & (GQR_ARA1(0.7,-1.3|log)),data=simData)
         theta<-c(0.3,2.2,0.7,0.4,-0.9)
         #L<-logLikelihood(mle,theta,c(TRUE,FALSE,FALSE))
         
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         rhoMC<-theta[3]
         rhoMP<-theta[4]
         rhoMP2<-theta[5]
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         Lcalc<-0
         T<-c(0,simData$Time[simData$System==1])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i]))+V; A<-rhoMP^(log(2)); C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i]))+V; A<-rhoMP^(log(3)); C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i]))+V; A<-rhoMP^(log(4)); C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==2])
         V<-0; C<-0; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==3])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i]))+V; A<-rhoMP^(log(2)); C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==4])
         V<-0; C<-0; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i]))+V; A<-rhoMP^(log(2));  C<-1; i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i]))+V; A<-rhoMP^(log(3)); C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==5])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i]))+V; A<-rhoMP^(log(2)); C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i]))+V; A<-rhoMP^(log(3)); C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-rhoMP2)*(A*(T[i+1]-T[i]))+V; A<-rhoMP^(log(4)); C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         fix<-rep(TRUE,length(theta))
         fix[1]=FALSE
         alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
         Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
       },
       TGQRARA_2={
         #Weibull + GQR_ARAInf
         simData<-data.frame(Time=c(18.09,52.07,95.71,145.75),Type=c(-1,-1,-1,-1),row.names=1:4)
         mle <- mle.vam(Time & Type ~ (GQR_ARAInf(0.1,0.5) | Weibull(0.001,2.5)),data=simData)
         theta<-c(0.03,2.4,0.7,-1.2)
         rho<-theta[3]
         rho2<-theta[4]
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         T<-simData$Time
         Lcalc<-log(h(T[1]))-H(T[1])
         V<-(1-rho2)*T[1]
         Lcalc<-Lcalc+log(rho^(1)*h(rho^(1)*(T[2]-T[1])+V))-(H(rho^(1)*(T[2]-T[1])+V)-H(V))
         V<-(1-rho2)*(rho^(1)*(T[2]-T[1])+V )
         Lcalc<-Lcalc+log(rho^(2)*h(rho^(2)*(T[3]-T[2])+V))-(H(rho^(2)*(T[3]-T[2])+V)-H(V))
         V<-(1-rho2)*(rho^(2)*(T[3]-T[2])+V)
         Lcalc<-Lcalc+log(rho^(3)*h(rho^(3)*(T[4]-T[3])+V))-(H(rho^(3)*(T[4]-T[3])+V)-H(V))
         fix<-rep(TRUE,length(theta))
         fix[1]=FALSE
         alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
         Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
       },
       TGQRARA_1={
         #Weibull + GQR_ARA1-log
         simData<-data.frame(Time=c(18.09,52.07,95.71,145.75),Type=c(-1,-1,-1,-1),row.names=1:4)
         mle <- mle.vam(Time & Type ~ (GQR_ARA1(0.7,-1.2|log) | Weibull(0.001,2.5)),data=simData)
         theta<-c(0.03,2.4,0.7,-1.2)
         rho<-theta[3]
         rho2<-theta[4]
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         T<-simData$Time
         Lcalc<-log(h(T[1]))-H(T[1])
         V<-(1-rho2)*T[1]
         Lcalc<-Lcalc+log(rho^(log(2))*h(rho^(log(2))*(T[2]-T[1])+V))-(H(rho^(log(2))*(T[2]-T[1])+V)-H(V))
         V<-rho^(log(2))*(T[2]-T[1])+V -rho2 *rho^(log(2))*(T[2]-T[1])
         Lcalc<-Lcalc+log(rho^(log(3))*h(rho^(log(3))*(T[3]-T[2])+V))-(H(rho^(log(3))*(T[3]-T[2])+V)-H(V))
         V<-rho^(log(3))*(T[3]-T[2])+V -rho2* rho^(log(3))*(T[3]-T[2])
         Lcalc<-Lcalc+log(rho^(log(4))*h(rho^(log(4))*(T[4]-T[3])+V))-(H(rho^(log(4))*(T[4]-T[3])+V)-H(V))
         fix<-rep(TRUE,length(theta))
         fix[1]=FALSE
         alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
         Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
       },
       TGQR3={
         #Weibull + MC ARA1 + MP GQR-log + MultiSystems + Censorship
         simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
         #simData<-simData[1:12,]
         mle <- mle.vam(System & Time & Type ~ (ARA1(0.2) | Weibull(0.001,2.5)) & (GQR(0.7|log)),data=simData)
         theta<-c(0.3,2.2,0.7,0.4)
         #L<-logLikelihood(mle,theta,c(TRUE,FALSE,FALSE))
         
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         rhoMC<-theta[3]
         rhoMP<-theta[4]
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         Lcalc<-0
         T<-c(0,simData$Time[simData$System==1])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(log(2)); C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(log(3)); C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(log(4)); C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==2])
         V<-0; C<-0; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==3])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(log(2)); C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==4])
         V<-0; C<-0; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(log(2));  C<-1; i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(log(3)); C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==5])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(log(2)); C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(log(3)); C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(log(4)); C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         fix<-rep(TRUE,length(theta))
         fix[1]=FALSE
         alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
         Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
       },
       TGQR2={
         #Weibull + GQR-sqrt
         simData<-data.frame(Time=c(18.09,52.07,95.71,145.75),Type=c(-1,-1,-1,-1),row.names=1:4)
         mle <- mle.vam(Time & Type ~ (GQR(0.7|sqrt) | Weibull(0.001,2.5)),data=simData)
         theta<-c(0.03,2.4,0.7)
         
         rho<-theta[3]
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         T<-simData$Time
         Lcalc<-log(h(T[1]))-H(T[1])
         Lcalc<-Lcalc+log(rho^(sqrt(1))*h(rho^(sqrt(1))*(T[2]-T[1])))-(H(rho^(sqrt(1))*(T[2]-T[1]))-H(rho^(sqrt(1))*(T[1]-T[1])))
         Lcalc<-Lcalc+log(rho^(sqrt(2))*h(rho^(sqrt(2))*(T[3]-T[2])))-(H(rho^(sqrt(2))*(T[3]-T[2]))-H(rho^(sqrt(2))*(T[2]-T[2])))
         Lcalc<-Lcalc+log(rho^(sqrt(3))*h(rho^(sqrt(3))*(T[4]-T[3])))-(H(rho^(sqrt(3))*(T[4]-T[3]))-H(rho^(sqrt(3))*(T[3]-T[3])))
         fix<-rep(TRUE,length(theta))
         fix[1]=FALSE
         alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
         Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
       },
       TGQR1={
         #Weibull + GQR-log
         simData<-data.frame(Time=c(18.09,52.07,95.71,145.75),Type=c(-1,-1,-1,-1),row.names=1:4)
         mle <- mle.vam(Time & Type ~ (GQR(0.7|log) | Weibull(0.001,2.5)),data=simData)
         theta<-c(0.03,2.4,0.7)
         
         rho<-theta[3]
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         T<-simData$Time
         Lcalc<-log(h(T[1]))-H(T[1])
         Lcalc<-Lcalc+log(rho^(log(2))*h(rho^(log(2))*(T[2]-T[1])))-(H(rho^(log(2))*(T[2]-T[1]))-H(rho^(log(2))*(T[1]-T[1])))
         Lcalc<-Lcalc+log(rho^(log(3))*h(rho^(log(3))*(T[3]-T[2])))-(H(rho^(log(3))*(T[3]-T[2]))-H(rho^(log(3))*(T[2]-T[2])))
         Lcalc<-Lcalc+log(rho^(log(4))*h(rho^(log(4))*(T[4]-T[3])))-(H(rho^(log(4))*(T[4]-T[3]))-H(rho^(log(4))*(T[3]-T[3])))
         fix<-rep(TRUE,length(theta))
         fix[1]=FALSE
         alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
         Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
       },
       TQR10={
         #Weibull + CM ARAInf + PM QAGAN + PM QR + mutlisystems
         simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,14)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,2,2,1, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,1,2,2,1,0),row.names=1:22)
         mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (QAGAN()+QR(0.7)),data=simData)
         theta<-c(0.3,1.8,0.3,0.7)
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         rhoMC<-theta[3]
         rhoMP<-theta[4]
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         Lcalc<-0
         T<-c(0,simData$Time[simData$System==1])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-0; A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==2])
         V<-0; C<-0; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==3])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-0; A<-A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-0; A<-A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-0; A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         
         V<-0; A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-0; A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         fix<-rep(TRUE,length(theta))
         fix[1]=FALSE
         alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
         Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
       },
       TQR9={
         #Weibull + CM ARAInf + PM AGAP + PM QR + mutlisystems
         simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,14)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,2,2,1, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,1,2,2,1,0),row.names=1:22)
         mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (AGAP()+QR(0.7)),data=simData)
         theta<-c(0.3,1.8,0.3,0.7)
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         rhoMC<-theta[3]
         rhoMP<-theta[4]
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         Lcalc<-0
         T<-c(0,simData$Time[simData$System==1])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-V; A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==2])
         V<-0; C<-0; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==3])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-V; A<-A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-V; A<-A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-V; A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         
         V<-V; A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-V; A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         fix<-rep(TRUE,length(theta))
         fix[1]=FALSE
         alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
         Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
       },
       TQR8={
         #Weibull + CM ARAInf + PM ABAO + PM QR + mutlisystems
         simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,14)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,2,2,1, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,1,2,2,1,0),row.names=1:22)
         mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (ABAO()+QR(0.7)),data=simData)
         theta<-c(0.3,1.8,0.3,0.7)
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         rhoMC<-theta[3]
         rhoMP<-theta[4]
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         Lcalc<-0
         T<-c(0,simData$Time[simData$System==1])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(A*(T[i+1]-T[i])+V); A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==2])
         V<-0; C<-0; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==3])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-A*(T[i+1]-T[i])+V; A<-A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-A*(T[i+1]-T[i])+V; A<-A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(A*(T[i+1]-T[i])+V); A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         
         V<-A*(T[i+1]-T[i])+V; A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-A*(T[i+1]-T[i])+V; A<-A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         fix<-rep(TRUE,length(theta))
         fix[1]=FALSE
         alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
         Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
       },
       TQR7={
         #Weibull + CM ARAInf + PM AGAN + PM QR + mutlisystems
         simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,14)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,2,2,1, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,1,2,2,1,0),row.names=1:22)
         mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (AGAN()+QR(0.7)),data=simData)
         theta<-c(0.3,1.8,0.3,0.7)
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         rhoMC<-theta[3]
         rhoMP<-theta[4]
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         Lcalc<-0
         T<-c(0,simData$Time[simData$System==1])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-0*(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-1+0*rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==2])
         V<-0; C<-0; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==3])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-0*(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-1+0*rhoMP^(C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-0*(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-1+0*rhoMP^(C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-0*(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-1+0*rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         
         V<-0*(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-1+0*rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-0*(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i])+V)); A<-1+0*rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         fix<-rep(TRUE,length(theta))
         fix[1]=FALSE
         alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
         Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
       },
       TQR6={
         #Weibull + MC ARA1 + MP QR + MultiSystems + Censorship
         simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
         #simData<-simData[1:12,]
         mle <- mle.vam(System & Time & Type ~ (ARA1(0.2) | Weibull(0.001,2.5)) & (QR(0.7)),data=simData)
         theta<-c(0.3,2.2,0.7,0.4)
         #L<-logLikelihood(mle,theta,c(TRUE,FALSE,FALSE))
         
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         rhoMC<-theta[3]
         rhoMP<-theta[4]
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         Lcalc<-0
         T<-c(0,simData$Time[simData$System==1])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==2])
         V<-0; C<-0; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==3])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==4])
         V<-0; C<-0; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(C)*A;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==5])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-(1-C)*((1-rhoMC)*(A*(T[i+1]-T[i]))+V); A<-rhoMP^(C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         fix<-rep(TRUE,length(theta))
         fix[1]=FALSE
         alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
         Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
       },
       TQR5={
         #Weibull + MC QR + MP ARA1 + MultiSystems + Censorship
         simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
         #simData<-simData[1:12,]
         mle <- mle.vam(System & Time & Type ~ (QR(0.2) | Weibull(0.001,2.5)) & (ARA1(0.7)),data=simData)
         theta<-c(0.3,2.2,0.7,0.4)
         #L<-logLikelihood(mle,theta,c(TRUE,FALSE,FALSE))
         
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         rhoMC<-theta[3]
         rhoMP<-theta[4]
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         Lcalc<-0
         T<-c(0,simData$Time[simData$System==1])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*((1-rhoMP)*(A*(T[i+1]-T[i]))+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*((1-rhoMP)*(A*(T[i+1]-T[i]))+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*((1-rhoMP)*(A*(T[i+1]-T[i]))+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==2])
         V<-0; C<-0; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*((1-rhoMP)*(A*(T[i+1]-T[i]))+V); A<-rhoMC^(1-C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*((1-rhoMP)*(A*(T[i+1]-T[i]))+V); A<-rhoMC^(1-C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*((1-rhoMP)*(A*(T[i+1]-T[i]))+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==3])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*((1-rhoMP)*(A*(T[i+1]-T[i]))+V); A<-rhoMC^(1-C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*((1-rhoMP)*(A*(T[i+1]-T[i]))+V); A<-rhoMC^(1-C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*((1-rhoMP)*(A*(T[i+1]-T[i]))+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==4])
         V<-0; C<-0; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*((1-rhoMP)*(A*(T[i+1]-T[i]))+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*((1-rhoMP)*(A*(T[i+1]-T[i]))+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*((1-rhoMP)*(A*(T[i+1]-T[i]))+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==5])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*((1-rhoMP)*(A*(T[i+1]-T[i]))+V); A<-rhoMC^(1-C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*((1-rhoMP)*(A*(T[i+1]-T[i]))+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*((1-rhoMP)*(A*(T[i+1]-T[i]))+V); A<-rhoMC^(1-C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*((1-rhoMP)*(A*(T[i+1]-T[i]))+V); A<-rhoMC^(1-C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*((1-rhoMP)*(A*(T[i+1]-T[i]))+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*((1-rhoMP)*(A*(T[i+1]-T[i]))+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         fix<-rep(TRUE,length(theta))
         fix[1]=FALSE
         alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
         Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
       }, 
      TQR4={
         #Weibull + MC QR + MP ARAinf + MultiSystems + Censorship
         simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
         #simData<-simData[1:12,]
         mle <- mle.vam(System & Time & Type ~ (QR(0.2) | Weibull(0.001,2.5)) & (ARAInf(0.7)),data=simData)
         theta<-c(0.3,2.2,0.7,0.4)
         #L<-logLikelihood(mle,theta,c(TRUE,FALSE,FALSE))
         
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         rhoMC<-theta[3]
         rhoMP<-theta[4]
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         Lcalc<-0
         T<-c(0,simData$Time[simData$System==1])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*(1-rhoMP)*(A*(T[i+1]-T[i])+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*(1-rhoMP)*(A*(T[i+1]-T[i])+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*(1-rhoMP)*(A*(T[i+1]-T[i])+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==2])
         V<-0; C<-0; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*(1-rhoMP)*(A*(T[i+1]-T[i])+V); A<-rhoMC^(1-C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*(1-rhoMP)*(A*(T[i+1]-T[i])+V); A<-rhoMC^(1-C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*(1-rhoMP)*(A*(T[i+1]-T[i])+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==3])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*(1-rhoMP)*(A*(T[i+1]-T[i])+V); A<-rhoMC^(1-C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*(1-rhoMP)*(A*(T[i+1]-T[i])+V); A<-rhoMC^(1-C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*(1-rhoMP)*(A*(T[i+1]-T[i])+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==4])
         V<-0; C<-0; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*(1-rhoMP)*(A*(T[i+1]-T[i])+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*(1-rhoMP)*(A*(T[i+1]-T[i])+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*(1-rhoMP)*(A*(T[i+1]-T[i])+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         T<-c(0,simData$Time[simData$System==5])
         V<-0; C<-1; A<-1; i<-1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*(1-rhoMP)*(A*(T[i+1]-T[i])+V); A<-rhoMC^(1-C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*(1-rhoMP)*(A*(T[i+1]-T[i])+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*(1-rhoMP)*(A*(T[i+1]-T[i])+V); A<-rhoMC^(1-C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*(1-rhoMP)*(A*(T[i+1]-T[i])+V); A<-rhoMC^(1-C)*A; C<-0;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*(1-rhoMP)*(A*(T[i+1]-T[i])+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         V<-C*(1-rhoMP)*(A*(T[i+1]-T[i])+V); A<-rhoMC^(1-C)*A; C<-1;  i<-i+1
         Lcalc<-Lcalc+(1-C)*log(A*h(A*(T[i+1]-T[i])+V))-(H(A*(T[i+1]-T[i])+V)-H(V))
         fix<-rep(TRUE,length(theta))
         fix[1]=FALSE
         alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
         Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
       }, 
    TQR3={
         #Weibull + MC ARAInf + MP QR + MultiSystems + Censorship
      simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
      #simData<-simData[1:12,]
      mle <- mle.vam(System & Time & Type ~ (ARAInf(0.2) | Weibull(0.001,2.5)) & (QR(0.7)),data=simData)
      theta<-c(0.3,2.2,0.3,0.4)
      #L<-logLikelihood(mle,theta,c(TRUE,FALSE,FALSE))
      
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         rhoMC<-theta[3]
         rhoMP<-theta[4]
         h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
         H<-function(t) theta[1]*t^(theta[2])
         T<-simData$Time[simData$System==1]
         rho<-rhoMP
         Lcalc<--H(T[1])
         Lcalc<-Lcalc-H(rho*(T[2]-T[1]))
         Lcalc<-Lcalc-H(rho^2*(T[3]-T[2]))
         Lcalc<-Lcalc-H(rho^3*(T[4]-T[3]))
          T<-simData$Time[simData$System==2]
          rho<-rhoMC
          Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
          Lcalc<-Lcalc+log(h(T[2]-rho*T[1]))-(H(T[2]-rho*T[1])-H(T[1]-rho*T[1]))
          Lcalc<-Lcalc+log(h(T[3]-rho*T[2]-rho*(1-rho)*T[1]))-(H(T[3]-rho*T[2]-rho*(1-rho)*T[1])-H(T[2]-rho*T[2]-rho*(1-rho)*T[1]))
          Lcalc<-Lcalc-(H(T[4]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1])-H(T[3]-rho*T[3]-rho*(1-rho)*T[2]-rho*(1-rho)^2*T[1]))
          T<-simData$Time[simData$System==3]
          Lcalc<-Lcalc-H(T[1])
          Lcalc<-Lcalc+log(rhoMP*h(rhoMP*(T[2]-T[1])))-H(rhoMP*(T[2]-T[1]))
          V<-(1-rhoMC)*rhoMP*(T[2]-T[1])
          Lcalc<-Lcalc+log(rhoMP*h(rhoMP*(T[3]-T[2])+V))-(H(rhoMP*(T[3]-T[2])+V)-H(V))
          V<-(1-rhoMC)*(rhoMP*(T[3]-T[2])+V)
          Lcalc<-Lcalc-(H(rhoMP*(T[4]-T[3])+V)-H(V))
         T<-simData$Time[simData$System==4]
         Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
         Lcalc<-Lcalc-(H(T[2]-rhoMC*T[1])-H(T[1]-rhoMC*T[1]))
         Lcalc<-Lcalc-H(rhoMP*(T[3]-T[2]))
         Lcalc<-Lcalc-H(rhoMP^2*(T[4]-T[3]))
         T<-simData$Time[simData$System==5]
         Lcalc<-Lcalc-H(T[1])
         Lcalc<-Lcalc+log(rhoMP*h(rhoMP*(T[2]-T[1])))-H(rhoMP*(T[2]-T[1]))
         V<-(1-rhoMC)*rhoMP*(T[2]-T[1])
         Lcalc<-Lcalc-(H(rhoMP*(T[3]-T[2])+V)-H(V))
         Lcalc<-Lcalc+log(rhoMP^2*h(rhoMP^2*(T[4]-T[3])))-H(rhoMP^2*(T[4]-T[3]))
         V<-(1-rhoMC)*rhoMP^2*(T[4]-T[3])
         Lcalc<-Lcalc+log(rhoMP^2*h(rhoMP^2*(T[5]-T[4])+V))-(H(rhoMP^2*(T[5]-T[4])+V)-H(V))
         V<-(1-rhoMC)*(rhoMP^2*(T[5]-T[4])+V)
         Lcalc<-Lcalc-(H(rhoMP^2*(T[6]-T[5])+V)-H(V))
         Lcalc<-Lcalc-H(rhoMP^3*(T[7]-T[6]))
         fix<-rep(TRUE,length(theta))
         fix[1]=FALSE
         alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
         Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
       },     
  TQR2={
    #Weibull + QR + MultiSystems + Censorship
    simData<-data.frame(System=c(1,1,1,1,2,2,2,3),Time=c(18.09,52.07,95.71,145.75,15.02,45.1,82,20.1),Type=c(-1,-1,-1,0,-1,-1,-1,-1),row.names=1:8)
    mle <- mle.vam(System & Time & Type ~ (QR(0.7) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.03,2.4,0.7)
    
    rho<-theta[3]
    h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
    H<-function(t) theta[1]*t^(theta[2])
    T<-simData$Time[1:4]
    Lcalc<-log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(rho*h(rho*(T[2]-T[1])))-(H(rho*(T[2]-T[1]))-H(rho*(T[1]-T[1])))
    Lcalc<-Lcalc+log(rho^2*h(rho^2*(T[3]-T[2])))-(H(rho^2*(T[3]-T[2]))-H(rho^2*(T[2]-T[2])))
    Lcalc<-Lcalc-(H(rho^3*(T[4]-T[3]))-H(rho^3*(T[3]-T[3])))
    T<-simData$Time[5:7]
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(rho*h(rho*(T[2]-T[1])))-(H(rho*(T[2]-T[1]))-H(rho*(T[1]-T[1])))
    Lcalc<-Lcalc+log(rho^2*h(rho^2*(T[3]-T[2])))-(H(rho^2*(T[3]-T[2]))-H(rho^2*(T[2]-T[2])))
    T<-simData$Time[8]
    Lcalc<-Lcalc+log(h(T[1]))-H(T[1])
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
  TQR={
    #Weibull + QR
    simData<-data.frame(Time=c(18.09,52.07,95.71,145.75),Type=c(-1,-1,-1,-1),row.names=1:4)
    mle <- mle.vam(Time & Type ~ (QR(0.7) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.03,2.4,0.7)
    
    rho<-theta[3]
    h<-function(t) theta[1]*theta[2]*t^(theta[2]-1)
    H<-function(t) theta[1]*t^(theta[2])
    T<-simData$Time
    Lcalc<-log(h(T[1]))-H(T[1])
    Lcalc<-Lcalc+log(rho*h(rho*(T[2]-T[1])))-(H(rho*(T[2]-T[1]))-H(rho*(T[1]-T[1])))
    Lcalc<-Lcalc+log(rho^2*h(rho^2*(T[3]-T[2])))-(H(rho^2*(T[3]-T[2]))-H(rho^2*(T[2]-T[2])))
    Lcalc<-Lcalc+log(rho^3*h(rho^3*(T[4]-T[3])))-(H(rho^3*(T[4]-T[3]))-H(rho^3*(T[3]-T[3])))
    fix<-rep(TRUE,length(theta))
    fix[1]=FALSE
    alpha_Est<-(run(mle,fixed=fix,verbose=FALSE))[1]
    Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))
  },
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
print((L-Lcalc)/L)
print((C-Ccalc)/C)
print((dL-EstdL)/dL)
print((dC-EstdC)/dC)
print((d2L-Estd2L)/d2L)
print((d2C-Estd2C)/d2C)

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
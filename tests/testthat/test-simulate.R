
context("simulate")

test_that("Weibull",{
	sim <- sim.vam( ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
	set.seed(0.5)
	U<-runif(1)
	Tcour<-(-(0.3)^(-1)*log(U[1]))^(1/2.4)
	Data<-data.frame(Time=Tcour,Type=c(-1))
	set.seed(0.5)
	expect_that(simulate(sim,1),equals(Data,tolerance=0.00000000000001))
	Data_list<-list(Data)
	names(Data_list)<-"System1"
	set.seed(0.5)
	expect_that(simulate(sim,1,as.list=TRUE),equals(Data_list,tolerance=0.00000000000001))
	sim2 <- sim.vam( T & C ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
	Data2<-data.frame(T=Tcour,C=c(-1))
	set.seed(0.5)
	expect_that(simulate(sim2,1),equals(Data2,tolerance=0.00000000000001))
	Data2_list<-list(Data2)
	names(Data2_list)<-"System1"
	set.seed(0.5)
	expect_that(simulate(sim2,1,as.list=TRUE),equals(Data2_list,tolerance=0.00000000000001))
	names(Data2_list)<-"S1"
	sim2 <- sim.vam(S & T & C ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
	set.seed(0.5)
	expect_that(simulate(sim2,1,as.list=TRUE),equals(Data2_list,tolerance=0.00000000000001))
}
)

test_that("Weibull + MultiSystems",{
	sim <- sim.vam( ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
	set.seed(0.5)
	U<-runif(1)
	Tcour1<-(-(0.3)^(-1)*log(U[1]))^(1/2.4)
	U<-runif(1)
	Tcour2<-(-(0.3)^(-1)*log(U[1]))^(1/2.4)
	U<-runif(1)
	Tcour3<-(-(0.3)^(-1)*log(U[1]))^(1/2.4)
	Data<-data.frame(System=1:3,Time=c(Tcour1,Tcour2,Tcour3),Type=c(-1,-1,-1))
	set.seed(0.5)
	expect_that(simulate(sim,1,nb.system=3),equals(Data,tolerance=0.00000000000001))
	Data1<-data.frame(Time=Tcour1,Type=c(-1))
	Data2<-data.frame(Time=Tcour2,Type=c(-1))
	Data3<-data.frame(Time=Tcour3,Type=c(-1))
	Data_list<-list(Data1,Data2,Data3)
	names(Data_list)<-c("System1","System2","System3")
	set.seed(0.5)
	expect_that(simulate(sim,1,as.list=TRUE,nb.system=3),equals(Data_list,tolerance=0.00000000000001))
	sim2 <- sim.vam( T & C ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
	Data1<-data.frame(T=Tcour1,C=c(-1))
	Data2<-data.frame(T=Tcour2,C=c(-1))
	Data3<-data.frame(T=Tcour3,C=c(-1))
	Data_list<-list(Data1,Data2,Data3)
	names(Data_list)<-c("System1","System2","System3")
	set.seed(0.5)
	expect_that(simulate(sim2,1,nb.system=3,as.list=TRUE),equals(Data_list,tolerance=0.00000000000001))
	Data<-data.frame(System=1:3,T=c(Tcour1,Tcour2,Tcour3),C=c(-1,-1,-1))
	set.seed(0.5)
	expect_that(simulate(sim2,1,nb.system=3),equals(Data,tolerance=0.00000000000001))
	names(Data_list)<-c("S1","S2","S3")
	sim2 <- sim.vam(S & T & C ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
	set.seed(0.5)
	expect_that(simulate(sim2,1,nb.system=3,as.list=TRUE),equals(Data_list,tolerance=0.00000000000001))
	Data<-data.frame(S=1:3,T=c(Tcour1,Tcour2,Tcour3),C=c(-1,-1,-1))
	set.seed(0.5)
	expect_that(simulate(sim2,1,nb.system=3),equals(Data,tolerance=0.00000000000001))
}
)

test_that("loglinear",{
	sim <- sim.vam( ~ (ARAInf(0.6) | LogLinear(0.01,0.8)))
	set.seed(0.5)
	U<-runif(1)
	Tcour<-0.8^(-1)*log(1+0.8*(-(0.01)^(-1)*log(U[1])))
	Data<-data.frame(Time=Tcour,Type=c(-1))
	set.seed(0.5)
	expect_that(simulate(sim,1),equals(Data,tolerance=0.00000000000001))
}
)

test_that("Weibull+ARAInf",{
	sim <- sim.vam( ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
	rho<-0.6
	set.seed(0.5)
	T<-c()
	Vr<-0
	Tprec<-0
	for(i in 1:5){
		U<-runif(1)
		Tcour<-Tprec+(Vr^2.4-(0.3)^(-1)*log(U[1]))^(1/2.4)-Vr
		T<-c(T,Tcour)
		Vr<-(1-rho)*(Vr+Tcour-Tprec)
		Tprec<-Tcour
	}
	Data<-data.frame(Time=T,Type=rep(-1,5))
	set.seed(0.5)
	expect_that(simulate(sim,5),equals(Data,tolerance=0.00000000000001))
	Data_list<-list(Data)
	names(Data_list)<-"System1"
	set.seed(0.5)
	expect_that(simulate(sim,5,as.list=TRUE),equals(Data_list,tolerance=0.00000000000001))
	sim2 <- sim.vam( T & C ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
	Data2<-data.frame(T=T,C=rep(-1,5))
	set.seed(0.5)
	expect_that(simulate(sim2,5),equals(Data2,tolerance=0.00000000000001))
	Data2_list<-list(Data2)
	names(Data2_list)<-"System1"
	set.seed(0.5)
	expect_that(simulate(sim2,5,as.list=TRUE),equals(Data2_list,tolerance=0.00000000000001))
	names(Data2_list)<-"S1"
	sim2 <- sim.vam(S & T & C ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
	set.seed(0.5)
	expect_that(simulate(sim2,5,as.list=TRUE),equals(Data2_list,tolerance=0.00000000000001))
}
)

test_that("Weibull+ARAInf+Mulisystems",{
	sim <- sim.vam( ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
	rho<-0.6
	set.seed(0.5)
	Cens<-4.1
	T<-c()
	Vr<-0
	Tprec<-0
	while(1){
		U<-runif(1)
		Tcour<-Tprec+(Vr^2.4-(0.3)^(-1)*log(U[1]))^(1/2.4)-Vr
		if(Tcour>Cens) break
		T<-c(T,Tcour)
		Vr<-(1-rho)*(Vr+Tcour-Tprec)
		Tprec<-Tcour
	}
	T1<-c(T,Cens)
	T<-c()
	Vr<-0
	Tprec<-0
	while(1){
		U<-runif(1)
		Tcour<-Tprec+(Vr^2.4-(0.3)^(-1)*log(U[1]))^(1/2.4)-Vr
		if(Tcour>Cens) break
		T<-c(T,Tcour)
		Vr<-(1-rho)*(Vr+Tcour-Tprec)
		Tprec<-Tcour
	}
	T2<-c(T,Cens)
	Data<-data.frame(System=c(rep(1,length(T1)),rep(2,length(T2))),Time=c(T1,T2),Type=c(rep(-1,length(T1)-1),0,rep(-1,length(T2)-1),0))
	set.seed(0.5)
	expect_that(simulate(sim,T>(RC=4.1),nb.system=2),equals(Data,tolerance=0.00000000000001))
	Data1<-data.frame(Time=T1,Type=c(rep(-1,length(T1)-1),0))
	Data2<-data.frame(Time=T2,Type=c(rep(-1,length(T2)-1),0))
	Data_list<-list(Data1,Data2)
	names(Data_list)<-c("System1","System2")
	set.seed(0.5)
	expect_that(simulate(sim,T>(RC=4.1),nb.system=2,as.list=TRUE),equals(Data_list,tolerance=0.00000000000001))
	sim2 <- sim.vam( T & C ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
	Data2<-data.frame(System=c(rep(1,length(T1)),rep(2,length(T2))),T=c(T1,T2),C=c(rep(-1,length(T1)-1),0,rep(-1,length(T2)-1),0))
	set.seed(0.5)
	expect_that(simulate(sim2,T>(RC=4.1),nb.system=2),equals(Data2,tolerance=0.00000000000001))
	Data1<-data.frame(T=T1,C=c(rep(-1,length(T1)-1),0))
	Data2<-data.frame(T=T2,C=c(rep(-1,length(T2)-1),0))
	Data2_list<-list(Data1,Data2)
	names(Data2_list)<-c("System1","System2")
	set.seed(0.5)
	expect_that(simulate(sim2,T>(RC=4.1),nb.system=2,as.list=TRUE),equals(Data2_list,tolerance=0.00000000000001))
	names(Data2_list)<-c("S1","S2")
	sim2 <- sim.vam(S & T & C ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
	set.seed(0.5)
	expect_that(simulate(sim2,T>(RC=4.1),nb.system=2,as.list=TRUE),equals(Data2_list,tolerance=0.00000000000001))
}
)

test_that("Weibull+QR",{
	sim <- sim.vam( ~ (QR(0.9) | Weibull(0.3,2.4)))
	rho<-0.9
	set.seed(0.5)
	Cens<-4.5
	T<-c()
	A<-1
	Tprec<-0
	while(1){
		U<-runif(1)
		Tcour<-Tprec+((-(0.3)^(-1)*log(U[1]))^(1/2.4))/A
		if(Tcour>Cens) break
		T<-c(T,Tcour)
		A<-A*rho
		Tprec<-Tcour
	}
	T1<-c(T,Cens)
	T<-c()
	A<-1
	Tprec<-0
	while(1){
		U<-runif(1)
		Tcour<-Tprec+((-(0.3)^(-1)*log(U[1]))^(1/2.4))/A
		if(Tcour>Cens) break
		T<-c(T,Tcour)
		A<-A*rho
		Tprec<-Tcour
	}
	T2<-c(T,Cens)
	Data<-data.frame(System=c(rep(1,length(T1)),rep(2,length(T2))),Time=c(T1,T2),Type=c(rep(-1,length(T1)-1),0,rep(-1,length(T2)-1),0))
	set.seed(0.5)
	expect_that(simulate(sim,T>(RC=4.5),nb.system=2),equals(Data,tolerance=0.00000000000001))
}
)

test_that("Weibull+GQR_ARA1",{
	sim <- sim.vam( ~ (GQR_ARA1(1.2,0.6|sqrt) | Weibull(0.3,2.4)))
	rho1<-1.2
	rho2<-0.6
	set.seed(0.5)
	Cens<-4.5
	T<-c()
	Vr<-0
	A<-1
	Tprec<-0
	i<-0
	while(1){
		U<-runif(1)
		Tcour<-Tprec+((Vr^2.4-(0.3)^(-1)*log(U[1]))^(1/2.4)-Vr)/A
		if(Tcour>Cens) break
		T<-c(T,Tcour)
		i<-i+1
		Vr<-Vr+(1-rho2)*A*(Tcour-Tprec)
		A<-rho1^(sqrt(i))
		Tprec<-Tcour
	}
	T1<-c(T,Cens)
	T<-c()
	Vr<-0
	A<-1
	Tprec<-0
	i<-0
	while(1){
		U<-runif(1)
		Tcour<-Tprec+((Vr^2.4-(0.3)^(-1)*log(U[1]))^(1/2.4)-Vr)/A
		if(Tcour>Cens) break
		T<-c(T,Tcour)
		i<-i+1
		Vr<-Vr+(1-rho2)*A*(Tcour-Tprec)
		A<-rho1^(sqrt(i))
		Tprec<-Tcour
	}
	T2<-c(T,Cens)
	Data<-data.frame(System=c(rep(1,length(T1)),rep(2,length(T2))),Time=c(T1,T2),Type=c(rep(-1,length(T1)-1),0,rep(-1,length(T2)-1),0))
	T<-c()
	Vr<-0
	A<-1
	Tprec<-0
	i<-0
	while(1){
		U<-runif(1)
		Tcour<-Tprec+((Vr^2.4-(0.3)^(-1)*log(U[1]))^(1/2.4)-Vr)/A
		if(Tcour>Cens) break
		T<-c(T,Tcour)
		i<-i+1
		Vr<-Vr+(1-rho2)*A*(Tcour-Tprec)
		A<-rho1^(sqrt(i))
		Tprec<-Tcour
	}
	T1<-c(T,Cens)
	T<-c()
	Vr<-0
	A<-1
	Tprec<-0
	i<-0
	while(1){
		U<-runif(1)
		Tcour<-Tprec+((Vr^2.4-(0.3)^(-1)*log(U[1]))^(1/2.4)-Vr)/A
		if(Tcour>Cens) break
		T<-c(T,Tcour)
		i<-i+1
		Vr<-Vr+(1-rho2)*A*(Tcour-Tprec)
		A<-rho1^(sqrt(i))
		Tprec<-Tcour
	}
	T2<-c(T,Cens)
	Data2<-data.frame(System=c(rep(1,length(T1)),rep(2,length(T2))),Time=c(T1,T2),Type=c(rep(-1,length(T1)-1),0,rep(-1,length(T2)-1),0))
	set.seed(0.5)
	expect_that(simulate(sim,T>(RC=4.5),nb.system=2),equals(Data,tolerance=0.00000000000001))
	expect_that(simulate(sim,T>(RC=4.5),nb.system=2),equals(Data2,tolerance=0.00000000000001))
}
)

test_that("Weibull+GQR_ARAm",{
	sim <- sim.vam( ~ (GQR_ARAm(1.2,0.6|sqrt,2) | Weibull(0.3,2.4)))
	rho1<-1.2
	rho2<-0.6
	set.seed(0.5)
	Cens<-4.5
	T<-c()
	Vr<-0
	A<-1;A1<-1
	Tprec<-0;Tprec1<-0
	i<-0
	while(1){
		U<-runif(1)
		Tcour<-Tprec+((Vr^2.4-(0.3)^(-1)*log(U[1]))^(1/2.4)-Vr)/A
		if(Tcour>Cens) break
		T<-c(T,Tcour)
		i<-i+1
		Vr<-Vr+(1-rho2)*A*(Tcour-Tprec)-rho2*(1-rho2)*A1*(Tprec-Tprec1)
		A1<-A
		A<-rho1^(sqrt(i))
		Tprec1<-Tprec
		Tprec<-Tcour
	}
	T1<-c(T,Cens)
	T<-c()
	Vr<-0
	A<-1;A1<-1
	Tprec<-0;Tprec1<-0
	i<-0
	while(1){
		U<-runif(1)
		Tcour<-Tprec+((Vr^2.4-(0.3)^(-1)*log(U[1]))^(1/2.4)-Vr)/A
		if(Tcour>Cens) break
		T<-c(T,Tcour)
		i<-i+1
		Vr<-Vr+(1-rho2)*A*(Tcour-Tprec)-rho2*(1-rho2)*A1*(Tprec-Tprec1)
		A1<-A
		A<-rho1^(sqrt(i))
		Tprec1<-Tprec
		Tprec<-Tcour
	}
	T2<-c(T,Cens)
	Data<-data.frame(System=c(rep(1,length(T1)),rep(2,length(T2))),Time=c(T1,T2),Type=c(rep(-1,length(T1)-1),0,rep(-1,length(T2)-1),0))
	T<-c()
	Vr<-0
	A<-1;A1<-1
	Tprec<-0;Tprec1<-0
	i<-0
	while(1){
		U<-runif(1)
		Tcour<-Tprec+((Vr^2.4-(0.3)^(-1)*log(U[1]))^(1/2.4)-Vr)/A
		if(Tcour>Cens) break
		T<-c(T,Tcour)
		i<-i+1
		Vr<-Vr+(1-rho2)*A*(Tcour-Tprec)-rho2*(1-rho2)*A1*(Tprec-Tprec1)
		A1<-A
		A<-rho1^(sqrt(i))
		Tprec1<-Tprec
		Tprec<-Tcour
	}
	T1<-c(T,Cens)
	T<-c()
	Vr<-0
	A<-1;A1<-1
	Tprec<-0;Tprec1<-0
	i<-0
	while(1){
		U<-runif(1)
		Tcour<-Tprec+((Vr^2.4-(0.3)^(-1)*log(U[1]))^(1/2.4)-Vr)/A
		if(Tcour>Cens) break
		T<-c(T,Tcour)
		i<-i+1
		Vr<-Vr+(1-rho2)*A*(Tcour-Tprec)-rho2*(1-rho2)*A1*(Tprec-Tprec1)
		A1<-A
		A<-rho1^(sqrt(i))
		Tprec1<-Tprec
		Tprec<-Tcour
	}
	T2<-c(T,Cens)
	Data2<-data.frame(System=c(rep(1,length(T1)),rep(2,length(T2))),Time=c(T1,T2),Type=c(rep(-1,length(T1)-1),0,rep(-1,length(T2)-1),0))
	set.seed(0.5)
	expect_that(simulate(sim,T>(RC=4.5),nb.system=2),equals(Data,tolerance=0.00000000000001))
	expect_that(simulate(sim,T>(RC=4.5),nb.system=2),equals(Data2,tolerance=0.00000000000001))
}
)

test_that("Weibull+CM ARA1+PM GQR AtTimes_notcycle",{
set.seed(0.5)
PM<-c(10,15,20)
simCMPM<-sim.vam(  ~ (ARA1(-1.2) | Weibull(.003,2.6)) & (GQR(.9|log) | AtTimes(c(10,15,20),FALSE)))
Data<-simulate(simCMPM,16,nb.system=2)
Data1<-Data[Data$System==1,c(2,3)]
Data2<-Data[Data$System==2,c(2,3)]
Data3<-simulate(simCMPM,12)
H<-function(t) 0.003*t^(2.6)

test<-function(Data,PM){
  n<-dim(Data)[1]
  PMnotOK<-0
  U<-c()
  V<-0; A<-1;Tprec<-0;iPM<-0
  for(i in 1:n){
    Tcour<-Data$Time[i]
    if(Data$Type[i]<0){
      U<-c(U,exp(-(H(A*(Tcour-Tprec)+V)-H(V))))
      V<-V+(1-(-1.2))*A*(Tcour-Tprec)
    } else {
      iPM=iPM+1
      if((iPM>length(PM))||(Tcour!=PM[iPM])){PMnotOK<-PMnotOK+1}
      A<-0.9^(log(iPM+1))
      V<-0
    }
    Tprec<-Tcour
  }
  list(PMnotOK,U)
}
set.seed(0.5)
res1<-test(Data1,PM);U1<-res1[[2]]
U1bis<-runif(dim(Data1)[1])[Data1$Type==-1]
res2<-test(Data2,PM);U2<-res2[[2]]
U2bis<-runif(dim(Data2)[1])[Data2$Type==-1]
res3<-test(Data3,PM);U3<-res3[[2]]
U3bis<-runif(dim(Data3)[1])[Data3$Type==-1]

	expect_that(res1[[1]],equals(0,tolerance=0.00000000000001))
	expect_that(sort(Data1$Time),equals(Data1$Time,tolerance=0.00000000000001))
	expect_that(U1,equals(U1bis,tolerance=0.00000000000001))
	expect_that(res2[[1]],equals(0,tolerance=0.00000000000001))
	expect_that(sort(Data2$Time),equals(Data2$Time,tolerance=0.00000000000001))
	expect_that(U2,equals(U2bis,tolerance=0.00000000000001))
	expect_that(res3[[1]],equals(0,tolerance=0.00000000000001))
	expect_that(sort(Data3$Time),equals(Data3$Time,tolerance=0.00000000000001))
	expect_that(U3,equals(U3bis,tolerance=0.0000000000001))
}
)

test_that("LogLinear+CM ARA3+PM GQR-ARAInf AtTimes_cycle",{
set.seed(0.5)
PMbis<-c(10,15,20)
PM<-c(10,15,20,30,35,40,50,55,60)
simCMPM<-sim.vam(  ~ (ARAm(0.6|3) | LogLinear(.005,0.8)) & (GQR_ARAInf(1.2,0.2|sqrt) | AtTimes(c(10,15,20))))
Data<-simulate(simCMPM,20,nb.system=2)
Data1<-Data[Data$System==1,c(2,3)]
Data2<-Data[Data$System==2,c(2,3)]
Data3<-simulate(simCMPM,40)
H<-function(t) 0.005/0.8*exp(t*0.8)

test<-function(Data,PM){
  n<-dim(Data)[1]
  PMnotOK<-0
  U<-c()
  V<-0; A<-1;V2<-0;V3<-0
  Tprec<-0;Tprec2<-0;Tprec3<-0;iPM<-0
  for(i in 1:n){
    Tcour<-Data$Time[i]
    if(Data$Type[i]<0){
      U<-c(U,exp(-(H(A*(Tcour-Tprec)+V)-H(V))))
      V<-V+(1-0.6)*A*(Tcour-Tprec)-0.6*V2-0.6*V3
      V3<-(1-0.6)*V2
      V2<-(1-0.6)*A*(Tcour-Tprec)
    } else {
      iPM=iPM+1
      if((iPM>length(PM))||(Tcour!=PM[iPM])){PMnotOK<-PMnotOK+1}
      V<-(1-(0.2))*(A*(Tcour-Tprec)+V)
      V3<-(1-0.2)*V2
      V2<-(1-0.2)*A*(Tcour-Tprec)
      A<-1.2^(sqrt(iPM))
    }
    Tprec3<-Tprec2
    Tprec2<-Tprec
    Tprec<-Tcour
  }
  list(PMnotOK,U)
}
set.seed(0.5)
res1<-test(Data1,PM);U1<-res1[[2]]
U1bis<-runif(dim(Data1)[1])[Data1$Type==-1]
res2<-test(Data2,PM);U2<-res2[[2]]
U2bis<-runif(dim(Data2)[1])[Data2$Type==-1]
res3<-test(Data3,PM);U3<-res3[[2]]
U3bis<-runif(dim(Data3)[1])[Data3$Type==-1]

expect_that(res1[[1]],equals(0,tolerance=0.00000000000001))
expect_that(sort(Data1$Time),equals(Data1$Time,tolerance=0.00000000000001))
expect_that(U1,equals(U1bis,tolerance=0.00000000000001))
expect_that(res2[[1]],equals(0,tolerance=0.00000000000001))
expect_that(sort(Data2$Time),equals(Data2$Time,tolerance=0.00000000000001))
expect_that(U2,equals(U2bis,tolerance=0.00000000000001))
expect_that(res3[[1]],equals(0,tolerance=0.00000000000001))
expect_that(sort(Data3$Time),equals(Data3$Time,tolerance=0.00000000000001))
expect_that(U3,equals(U3bis,tolerance=0.0000000000001))
}
)

test_that("LogLinear+CM ARA3+PM GQR-ARA1 AtVirtualAge",{
set.seed(0.5)
simCMPM<-sim.vam(  ~ (ARAm(0.6|3) | LogLinear(.005,0.8)) & (GQR_ARA1(1.2,0.8|sqrt) | AtVirtualAge(6)))
Data<-simulate(simCMPM,20,nb.system=2)
Data1<-Data[Data$System==1,c(2,3)]
Data2<-Data[Data$System==2,c(2,3)]
Data3<-simulate(simCMPM,40)
H<-function(t) 0.005/0.8*exp(t*0.8)

test<-function(Data,n){
  PMnotOK<-0
  U<-c()
  V<-0; A<-1;V2<-0;V3<-0
  Tprec<-0;Tprec2<-0;Tprec3<-0;iPM<-0
  for(i in 1:n){
    Tcour<-Data$Time[i]
    if(Data$Type[i]<0){
      U<-c(U,exp(-(H(A*(Tcour-Tprec)+V)-H(V))))
      V<-V+(1-0.6)*A*(Tcour-Tprec)-0.6*V2-0.6*V3
      V3<-(1-0.6)*V2
      V2<-(1-0.6)*A*(Tcour-Tprec)
    } else {
      iPM=iPM+1
      if(abs(A*(Tcour-Tprec)+V-6)>0.0000001){PMnotOK<-PMnotOK+1}
      V<-(1-(0.8))*A*(Tcour-Tprec)+V
      V3<-V2
      V2<-(1-0.8)*A*(Tcour-Tprec)
      A<-1.2^(sqrt(iPM))
    }
    Tprec3<-Tprec2
    Tprec2<-Tprec
    Tprec<-Tcour
  }
  list(PMnotOK,U)
}
set.seed(0.5)
res1<-test(Data1,dim(Data1)[1]);U1<-res1[[2]]
U1bis<-runif(dim(Data1)[1])[Data1$Type==-1]
res2<-test(Data2,dim(Data2)[1]);U2<-res2[[2]]
U2bis<-runif(dim(Data2)[1])[Data2$Type==-1]
res3<-test(Data3,dim(Data3)[1]);U3<-res3[[2]]
U3bis<-runif(dim(Data3)[1])[Data3$Type==-1]

expect_that(res1[[1]],equals(0,tolerance=0.00000000000001))
expect_that(sort(Data1$Time),equals(Data1$Time,tolerance=0.00000000000001))
expect_that(U1,equals(U1bis,tolerance=0.00000000000001))
expect_that(res2[[1]],equals(0,tolerance=0.00000000000001))
expect_that(sort(Data2$Time),equals(Data2$Time,tolerance=0.00000000000001))
expect_that(U2,equals(U2bis,tolerance=0.00000000000001))
expect_that(res3[[1]],equals(0,tolerance=0.00000000000001))
expect_that(sort(Data3$Time),equals(Data3$Time,tolerance=0.00000000000001))
expect_that(U3,equals(U3bis,tolerance=0.0000000000001))
}
)

test_that("LogLinear+PM ARA3+CM GQR-ARA1 AtFailureProbability",{
set.seed(0.5)
simCMPM<-sim.vam(  ~ (GQR_ARA1(1.2,0.8|sqrt) | LogLinear(.005,0.8)) & (ARAm(0.6|3) | AtFailureProbability(0.7)))
Data<-simulate(simCMPM,20,nb.system=2)
Data1<-Data[Data$System==1,c(2,3)]
Data2<-Data[Data$System==2,c(2,3)]
Data3<-simulate(simCMPM,40)
H<-function(t) 0.005/0.8*exp(t*0.8)

test<-function(Data,n){
  PMnotOK<-0
  U<-c()
  V<-0; A<-1;V2<-0;V3<-0
  Tprec<-0;Tprec2<-0;Tprec3<-0;iPM<-0;iCM<-0
  for(i in 1:n){
    Tcour<-Data$Time[i]
    if(Data$Type[i]<0){
      iCM=iCM+1
      U<-c(U,exp(-(H(A*(Tcour-Tprec)+V)-H(V))))
      V<-(1-(0.8))*A*(Tcour-Tprec)+V
      V3<-V2
      V2<-(1-0.8)*A*(Tcour-Tprec)
      A<-1.2^(sqrt(iCM))
    } else {
      iPM=iPM+1
      if(abs(exp(-(H(A*(Tcour-Tprec)+V)-H(V)))-(1-0.7))>0.0000001){PMnotOK<-PMnotOK+1}
      V<-V+(1-0.6)*A*(Tcour-Tprec)-0.6*V2-0.6*V3
      V3<-(1-0.6)*V2
      V2<-(1-0.6)*A*(Tcour-Tprec)
    }
    Tprec3<-Tprec2
    Tprec2<-Tprec
    Tprec<-Tcour
  }
  list(PMnotOK,U)
}
set.seed(0.5)
res1<-test(Data1,dim(Data1)[1]);U1<-res1[[2]]
U1bis<-runif(dim(Data1)[1])[Data1$Type==-1]
res2<-test(Data2,dim(Data2)[1]);U2<-res2[[2]]
U2bis<-runif(dim(Data2)[1])[Data2$Type==-1]
res3<-test(Data3,dim(Data3)[1]);U3<-res3[[2]]
U3bis<-runif(dim(Data3)[1])[Data3$Type==-1]

expect_that(res1[[1]],equals(0,tolerance=0.00000000000001))
expect_that(sort(Data1$Time),equals(Data1$Time,tolerance=0.00000000000001))
expect_that(U1,equals(U1bis,tolerance=0.00000000000001))
expect_that(res2[[1]],equals(0,tolerance=0.00000000000001))
expect_that(sort(Data2$Time),equals(Data2$Time,tolerance=0.00000000000001))
expect_that(U2,equals(U2bis,tolerance=0.00000000000001))
expect_that(res3[[1]],equals(0,tolerance=0.00000000000001))
expect_that(sort(Data3$Time),equals(Data3$Time,tolerance=0.00000000000001))
expect_that(U3,equals(U3bis,tolerance=0.0000000000001))
}
)

test_that("LogLinear+CM GQR_ARA2-sqrt+PM GQR-ARA4-log AtIntensity",{
set.seed(0.5)
simCMPM<-sim.vam(  ~ (GQR_ARAm(1.5,0.4|sqrt,2) | LogLinear(.005,0.8)) & (GQR_ARAm(1.2,0.8|log,4) | AtIntensity(1.6)))
Data<-simulate(simCMPM,20,nb.system=2)
Data1<-Data[Data$System==1,c(2,3)]
Data2<-Data[Data$System==2,c(2,3)]
Data3<-simulate(simCMPM,40)
H<-function(t) 0.005/0.8*exp(t*0.8)
h<-function(t) 0.005*exp(t*0.8)

test<-function(Data,n){
  PMnotOK<-0
  U<-c()
  V<-0; A<-1;V2<-0;V3<-0;V4<-0
  Tprec<-0;Tprec2<-0;Tprec3<-0;iPM<-0;iCM<-0
  for(i in 1:n){
    Tcour<-Data$Time[i]
    if(Data$Type[i]<0){
      iCM=iCM+1
      U<-c(U,exp(-(H(A*(Tcour-Tprec)+V)-H(V))))
      V<-V+(1-(0.4))*A*(Tcour-Tprec)-0.4*V2
      V4<-V3
      V3<-(1-0.4)*V2
      V2<-(1-0.4)*A*(Tcour-Tprec)
      A<-1.5^(sqrt(iCM))*1.2^(log(iPM+1))
    } else {
      iPM=iPM+1
      if(abs(h(A*(Tcour-Tprec)+V)*A-1.6)>0.0000001){PMnotOK<-PMnotOK+1;print(h(A*(Tcour-Tprec)+V)*A)}
      V<-V+(1-0.8)*A*(Tcour-Tprec)-0.8*V2-0.8*V3-0.8*V4
      V4<-(1-0.8)*V3
      V3<-(1-0.8)*V2
      V2<-(1-0.8)*A*(Tcour-Tprec)
      A<-1.5^(sqrt(iCM))*1.2^(log(iPM+1))
    }
    Tprec3<-Tprec2
    Tprec2<-Tprec
    Tprec<-Tcour
  }
  list(PMnotOK,U)
}
set.seed(0.5)
res1<-test(Data1,dim(Data1)[1]);U1<-res1[[2]]
U1bis<-runif(dim(Data1)[1])[Data1$Type==-1]
res2<-test(Data2,dim(Data2)[1]);U2<-res2[[2]]
U2bis<-runif(dim(Data2)[1])[Data2$Type==-1]
res3<-test(Data3,dim(Data3)[1]);U3<-res3[[2]]
U3bis<-runif(dim(Data3)[1])[Data3$Type==-1]

expect_that(res1[[1]],equals(0,tolerance=0.00000000000001))
expect_that(sort(Data1$Time),equals(Data1$Time,tolerance=0.00000000000001))
expect_that(U1,equals(U1bis,tolerance=0.00000000000001))
expect_that(res2[[1]],equals(0,tolerance=0.00000000000001))
expect_that(sort(Data2$Time),equals(Data2$Time,tolerance=0.00000000000001))
expect_that(U2,equals(U2bis,tolerance=0.00000000000001))
expect_that(res3[[1]],equals(0,tolerance=0.00000000000001))
expect_that(sort(Data3$Time),equals(Data3$Time,tolerance=0.00000000000001))
expect_that(U3,equals(U3bis,tolerance=0.0000000000001))
}
)

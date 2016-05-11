
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


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

}
)
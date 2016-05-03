
context("simulate")

test_that("Weibull",{
	theta<-c(0.3,2.4,0.6)	
	sim <- sim.vam( ~ (ARAInf(theta[3]) | Weibull(theta[1],theta[2])))
	set.seed(0.5)
	U<-runif(1)
	Tcour<-(-(theta[1])^(-1)*log(U[1]))^(1/theta[2])
	Data<-data.frame(Time=Tcour,Type=c(-1))
	set.seed(0.5)
	expect_that(simulate(sim,1),equals(Data,tolerance=0.00000000000001))
	Data_list<-list(Data)
	names(Data_list)<-"System1"
	set.seed(0.5)
	expect_that(simulate(sim,1,as.list=TRUE),equals(Data_list,tolerance=0.00000000000001))
	sim2 <- sim.vam( T & C ~ (ARAInf(theta[3]) | Weibull(theta[1],theta[2])))
	Data2<-data.frame(T=Tcour,C=c(-1))
	set.seed(0.5)
	expect_that(simulate(sim2,1),equals(Data2,tolerance=0.00000000000001))
	Data2_list<-list(Data2)
	names(Data2_list)<-"System1"
	set.seed(0.5)
	expect_that(simulate(sim2,1,as.list=TRUE),equals(Data2_list,tolerance=0.00000000000001))

}
)
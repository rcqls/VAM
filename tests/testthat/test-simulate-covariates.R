context("simulate with covariates")

test_that("Weibull+ARAInf+2cov",{
	set.seed(3)
	n<-5
	m<-6
    SimCov<-data.frame(cov1=runif(m,4,6),cov2=rbinom(m,1,c(0.3,0.5)))
    set.seed(4)
	simARAC<-sim.vam(~(ARAInf(0.6)|Weibull(0.1,2.2|0.5*cov1+(-1)*cov2)),data.covariates=SimCov)
	(simDataC<-simulate(simARAC,n))
	set.seed(4)
	simData<-simDataC
	for(i in 1:m){
  		a<-0.1*exp(0.5*SimCov[i,1]-SimCov[i,2])
  		simARA<-sim.vam(~(ARAInf(0.6)|Weibull(a,2.2)))
  		Data<-simulate(simARA,n)
  		simData[(n*(i-1)+1):(n*i),2:3]<-Data
	}
	expect_that(simData,equals(simDataC,tolerance=0.00000000000001))
}
)

test_that("Weibull+ARAInf+1cov",{
	set.seed(3)
	n<-5
	m<-6
    SimCov<-data.frame(cov1=runif(m,4,6))
    set.seed(4)
	simARAC<-sim.vam(~(ARAInf(0.6)|Weibull(0.1,2.2|0.5*cov1)),data.covariates=SimCov)
	(simDataC<-simulate(simARAC,n))
	set.seed(4)
	simData<-simDataC
	for(i in 1:m){
  		a<-0.1*exp(0.5*SimCov[i,1])
  		simARA<-sim.vam(~(ARAInf(0.6)|Weibull(a,2.2)))
  		Data<-simulate(simARA,n)
  		simData[(n*(i-1)+1):(n*i),2:3]<-Data
	}
	expect_that(simData,equals(simDataC,tolerance=0.00000000000001))
}
)

test_that("Weibull+ARA1+1cov+1system",{
	set.seed(3)
	n<-5
	m<-1
    SimCov<-data.frame(cov1=runif(m,4,6))
    set.seed(4)
	simARAC<-sim.vam(~(ARA1(0.6)|Weibull(0.1,2.2|0.5*cov1)),data.covariates=SimCov)
	simDataC<-simulate(simARAC,n)
	set.seed(4)
	simData<-simDataC
  		a<-0.1*exp(0.5*SimCov[1,1])
  		simARA<-sim.vam(~(ARA1(0.6)|Weibull(a,2.2)))
  		Data<-simulate(simARA,n)
  		simData[1:n,1:2]<-Data
	expect_that(simData,equals(simDataC,tolerance=0.00000000000001))
}
)

test_that("Weibull+ARA1+2cov+1system",{
	set.seed(3)
    SimCov<-data.frame(cov1=-0.4,cov2=0.3)
    set.seed(4)
	simARAC<-sim.vam(~(ARA1(0.6)|Weibull(0.1,2.2|-0.2*cov1-0.1*cov2)),data.covariates=SimCov)
	simDataC<-simulate(simARAC,n)
	set.seed(4)
	simData<-simDataC
  		a<-0.1*exp(-0.2*SimCov[1,1]-0.1*SimCov[1,2])
  		simARA<-sim.vam(~(ARA1(0.6)|Weibull(a,2.2)))
  		Data<-simulate(simARA,n)
  		simData[1:n,1:2]<-Data
	expect_that(simData,equals(simDataC,tolerance=0.00000000000001))
}
)

test_that("Weibull+GQR_ARAm+2cov",{
	set.seed(3)
	n<-5
	m<-6
    SimCov<-data.frame(cov1=runif(m,4,6),cov2=rbinom(m,1,c(0.3,0.5)))
    set.seed(4)
	simARAC<-sim.vam(~(GQR_ARAm(0.9,0.6|log,2)|Weibull(0.1,2.2|0.5*cov1+(-1)*cov2)),data.covariates=SimCov)
	(simDataC<-simulate(simARAC,n))
	set.seed(4)
	simData<-simDataC
	for(i in 1:m){
  		a<-0.1*exp(0.5*SimCov[i,1]-SimCov[i,2])
  		simARA<-sim.vam(~(GQR_ARAm(0.9,0.6|log,2)|Weibull(a,2.2)))
  		Data<-simulate(simARA,n)
  		simData[(n*(i-1)+1):(n*i),2:3]<-Data
	}
	expect_that(simData,equals(simDataC,tolerance=0.00000000000001))
}
)

test_that("loglinear+QR+2cov",{
	set.seed(3)
	n<-5
	m<-6
    SimCov<-data.frame(cov1=runif(m,4,6),cov2=rbinom(m,1,c(0.3,0.5)))
    set.seed(4)
	simARAC<-sim.vam(~(QR(0.9)|LogLinear(0.01,0.8|0.5*cov1+(-1)*cov2)),data.covariates=SimCov)
	(simDataC<-simulate(simARAC,n))
	set.seed(4)
	simData<-simDataC
	for(i in 1:m){
  		a<-0.01*exp(0.5*SimCov[i,1]-SimCov[i,2])
  		simARA<-sim.vam(~(QR(0.9)|LogLinear(a,0.8)))
  		Data<-simulate(simARA,n)
  		simData[(n*(i-1)+1):(n*i),2:3]<-Data
	}
	expect_that(simData,equals(simDataC,tolerance=0.00000000000001))
}
)

test_that("loglinear+QR+1cov",{
	set.seed(3)
	n<-5
	m<-6
	SimCov<-data.frame(cov1=runif(m,4,6))
	set.seed(4)
	simARAC<-sim.vam(~(QR(0.9)|LogLinear(0.01,0.8|(-0.2)*cov1)),data.covariates=SimCov)
	(simDataC<-simulate(simARAC,n))
	set.seed(4)
	simData<-simDataC
	for(i in 1:m){
		a<-0.01*exp(-0.2*SimCov[i,1])
		simARA<-sim.vam(~(QR(0.9)|LogLinear(a,0.8)))
		Data<-simulate(simARA,n)
		simData[(n*(i-1)+1):(n*i),2:3]<-Data
		}
	expect_that(simData,equals(simDataC,tolerance=0.00000000000001))
}
)

test_that("Weibull+ARA1+PM_ARAInf+1cov",{
	set.seed(3)
	n<-8
	m<-6
	SimCov<-data.frame(cov1=runif(m,3,6))
	set.seed(4)
	simARAC<-sim.vam(~(ARA1(0.6)|Weibull(0.1,2.5|(-0.2)*cov1))&(ARAInf(0.7)|Periodic(3)),data.covariates=SimCov)
	(simDataC<-simulate(simARAC,n))
	set.seed(4)
	simData<-simDataC
	for(i in 1:m){
		a<-0.1*exp(-0.2*SimCov[i,1])
		simARA<-sim.vam(~(ARA1(0.6)|Weibull(a,2.5))&(ARAInf(0.7)|Periodic(3)))
		Data<-simulate(simARA,n)
		simData[(n*(i-1)+1):(n*i),2:3]<-Data
		}
	expect_that(simData,equals(simDataC,tolerance=0.00000000000001))
}
)

test_that("Weibull+ARA1+PM_GQR+2cov",{
	set.seed(3)
	n<-8
	m<-6
	SimCov<-data.frame(cov1=runif(m,3,6),cov2=rbinom(m,1,c(0.3,0.5)))
	set.seed(4)
	simARAC<-sim.vam(~(ARA1(0.6)|Weibull(0.1,2.5|(-0.2)*cov1+1.1*cov2))&(GQR(0.9|sqrt)|Periodic(3)),data.covariates=SimCov)
	(simDataC<-simulate(simARAC,n))
	set.seed(4)
	simData<-simDataC
	for(i in 1:m){
		a<-0.1*exp(-0.2*SimCov[i,1]+1.1*SimCov[i,2])
		simARA<-sim.vam(~(ARA1(0.6)|Weibull(a,2.5))&(GQR(0.9|sqrt)|Periodic(3)))
		Data<-simulate(simARA,n)
		simData[(n*(i-1)+1):(n*i),2:3]<-Data
		}
	expect_that(simData,equals(simDataC,tolerance=0.00000000000001))
}
)

test_that("Weibull+GQR_ARAm+PM_GQR+2cov",{
	set.seed(3)
	n<-8
	m<-6
	SimCov<-data.frame(cov1=runif(m,3,6),cov2=rbinom(m,1,c(0.3,0.5)))
	set.seed(4)
	simARAC<-sim.vam(~(GQR_ARAm(1.1,0.7|log,2)|Weibull(0.1,2.5|-0.2*cov1+1.1*cov2))&(GQR(0.9|sqrt)|Periodic(3)),data.covariates=SimCov)
	(simDataC<-simulate(simARAC,n))
	set.seed(4)
	simData<-simDataC
	for(i in 1:m){
		a<-0.1*exp(-0.2*SimCov[i,1]+1.1*SimCov[i,2])
		simARA<-sim.vam(~(GQR_ARAm(1.1,0.7|log,2)|Weibull(a,2.5))&(GQR(0.9|sqrt)|Periodic(3)))
		Data<-simulate(simARA,n)
		simData[(n*(i-1)+1):(n*i),2:3]<-Data
		}
	expect_that(simData,equals(simDataC,tolerance=0.00000000000001))
}
)

test_that("Weibull+ARA1+PM_GQR+2cov",{
	set.seed(3)
	n<-8
	m<-6
	SimCov<-data.frame(cov1=runif(m,3,6),cov2=rbinom(m,1,c(0.3,0.5)))
	set.seed(4)
	simARAC<-sim.vam(~(ARA1(0.6)|Weibull(0.1,2.5|(-0.2)*cov1+1.1*cov2))&(GQR(0.9|sqrt)|Periodic(3)),data.covariates=SimCov)
	(simDataC<-simulate(simARAC,n))
	set.seed(4)
	simData<-simDataC
	for(i in 1:m){
		a<-0.1*exp(-0.2*SimCov[i,1]+1.1*SimCov[i,2])
		simARA<-sim.vam(~(ARA1(0.6)|Weibull(a,2.5))&(GQR(0.9|sqrt)|Periodic(3)))
		Data<-simulate(simARA,n)
		simData[(n*(i-1)+1):(n*i),2:3]<-Data
		}
	expect_that(simData,equals(simDataC,tolerance=0.00000000000001))
}
)

test_that("cov+AtIntensity",{
	set.seed(3)
	n<-5
	m<-6
	SimCov<-data.frame(cov1=runif(m,3,6),cov2=rbinom(m,1,c(0.3,0.5)))
	set.seed(4)
	simARAC<-sim.vam(~(ARAm(0.6|2)|Weibull(0.1,2.5|(-0.2)*cov1+1.1*cov2))&(AGAN()|AtIntensity(0.75)),data.covariates=SimCov)
	(simDataC<-simulate(simARAC,n))
	set.seed(4)
	simData<-simDataC
	for(i in 1:m){
		a<-0.1*exp(-0.2*SimCov[i,1]+1.1*SimCov[i,2])
		simARA<-sim.vam(~(ARAm(0.6|2)|Weibull(a,2.5))&(AGAN()|AtIntensity(0.75)))
		Data<-simulate(simARA,n)
		simData[(n*(i-1)+1):(n*i),2:3]<-Data
		}
	expect_that(simData,equals(simDataC,tolerance=0.00000000000001))
}
)

test_that("cov+AtVirtualAge",{
	set.seed(3)
	n<-5
	m<-6
	SimCov<-data.frame(cov1=runif(m,3,6),cov2=rbinom(m,1,c(0.3,0.5)))
	set.seed(4)
	simARAC<-sim.vam(~(ARAm(0.6|2)|Weibull(0.1,2.5|(-0.2)*cov1+1.1*cov2))&(AGAN()|AtVirtualAge(1.7)),data.covariates=SimCov)
	(simDataC<-simulate(simARAC,n))
	set.seed(4)
	simData<-simDataC
	for(i in 1:m){
		a<-0.1*exp(-0.2*SimCov[i,1]+1.1*SimCov[i,2])
		simARA<-sim.vam(~(ARAm(0.6|2)|Weibull(a,2.5))&(AGAN()|AtVirtualAge(1.7)))
		Data<-simulate(simARA,n)
		simData[(n*(i-1)+1):(n*i),2:3]<-Data
		}
	expect_that(simData,equals(simDataC,tolerance=0.00000000000001))
}
)

test_that("cov+AtFailureProbability",{
	set.seed(3)
	n<-5
	m<-6
	SimCov<-data.frame(cov1=runif(m,3,6),cov2=rbinom(m,1,c(0.3,0.5)))
	set.seed(4)
	simARAC<-sim.vam(~(ARAm(0.6|2)|Weibull(0.1,2.5|(-0.2)*cov1+1.1*cov2))&(ARAInf(0.8)|AtFailureProbability(0.5)),data.covariates=SimCov)
	(simDataC<-simulate(simARAC,n))
	set.seed(4)
	simData<-simDataC
	for(i in 1:m){
		a<-0.1*exp(-0.2*SimCov[i,1]+1.1*SimCov[i,2])
		simARA<-sim.vam(~(ARAm(0.6|2)|Weibull(a,2.5))&(ARAInf(0.8)|AtFailureProbability(0.5)))
		Data<-simulate(simARA,n)
		simData[(n*(i-1)+1):(n*i),2:3]<-Data
		}
	expect_that(simData,equals(simDataC,tolerance=0.00000000000001))
}
)


test_that("Weibull+GQR_ARAm+PM_GQR+PM_GQR_ARAm+2cov+Periodic+AtIntensity",{
	set.seed(3)
	n<-10
	m<-6
	SimCov<-data.frame(cov1=runif(m,3,6),cov2=rbinom(m,1,c(0.3,0.5)))
	set.seed(4)
	simARAC<-sim.vam(~(GQR_ARAm(1.1,0.7|log,2)|Weibull(0.1,2.5|(-0.2)*cov1+1.1*cov2))&(GQR(0.9|sqrt)+GQR_ARAm(1.2,0.6|log,3)|Periodic(3)*AtIntensity(0.75)),data.covariates=SimCov)
	(simDataC<-simulate(simARAC,n))
	set.seed(4)
	simData<-simDataC
	for(i in 1:m){
		a<-0.1*exp(-0.2*SimCov[i,1]+1.1*SimCov[i,2])
		simARA<-sim.vam(~(GQR_ARAm(1.1,0.7|log,2)|Weibull(a,2.5))&(GQR(0.9|sqrt)+GQR_ARAm(1.2,0.6|log,3)|Periodic(3)*AtIntensity(0.75)))
		Data<-simulate(simARA,n)
		simData[(n*(i-1)+1):(n*i),2:3]<-Data
		}
	expect_that(simData,equals(simDataC,tolerance=0.00000000000001))
}
)

test_that("Weibull+GQR_ARAm+PM_GQR+PM_GQR_ARAm+2cov+PeriodicProb",{
	set.seed(3)
	n<-10
	m<-6
	SimCov<-data.frame(cov1=runif(m,3,6),cov2=rbinom(m,1,c(0.3,0.5)))
	set.seed(4)
	simARAC<-sim.vam(~(GQR_ARAm(1.1,0.7|log,2)|Weibull(0.1,2.5|(-0.2)*cov1+1.1*cov2))&(GQR(0.9|sqrt)+GQR_ARAm(1.2,0.6|log,3)|Periodic(3,prob=c(0.55,0.45))),data.covariates=SimCov)
	(simDataC<-simulate(simARAC,n))
	set.seed(4)
	simData<-simDataC
	for(i in 1:m){
		a<-0.1*exp(-0.2*SimCov[i,1]+1.1*SimCov[i,2])
		simARA<-sim.vam(~(GQR_ARAm(1.1,0.7|log,2)|Weibull(a,2.5))&(GQR(0.9|sqrt)+GQR_ARAm(1.2,0.6|log,3)|Periodic(3,prob=c(0.55,0.45))))
		Data<-simulate(simARA,n)
		simData[(n*(i-1)+1):(n*i),2:3]<-Data
		}
	expect_that(simData,equals(simDataC,tolerance=0.00000000000001))
}
)

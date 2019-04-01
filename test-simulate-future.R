context("simulate future of a data set")

test_that("Initialization data set empty",{
	sim1 <- sim.vam( ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
  sim2 <- sim.vam( ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
	set.seed(0.5)
	Data<-simulate(sim1, 5)
	set.seed(0.5)
	expect_that(simulate(sim2,5,data=c()),equals(Data,tolerance=0.00000000000001))
}
)

test_that("One system only CM",{
	sim1 <- sim.vam( ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
  sim2 <- sim.vam( ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
	set.seed(0.5)
	Data<-simulate(sim1, 10)
	set.seed(0.5)
  InitData<-simulate(sim2, 4)
  Data2<-simulate(sim2,10,data=InitData)
  expect_that(Data2[1:4,],equals(InitData,tolerance=0.00000000000001))
	expect_that(Data2,equals(Data,tolerance=0.00000000000001))

  set.seed(0.5)
  Data<-simulate(sim1, 10)
  set.seed(0.5)
  InitData<-simulate(sim2, 1)
  Data2<-simulate(sim2,10,data=InitData)
  expect_that(Data2[1,],equals(InitData,tolerance=0.00000000000001))
  expect_that(Data2,equals(Data,tolerance=0.00000000000001))

  sim1 <- sim.vam(T&U ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
  sim2 <- sim.vam( T&U ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
	set.seed(0.5)
	Data<-simulate(sim1, 10)
	set.seed(0.5)
  InitData<-simulate(sim2, 4)
  Data2<-simulate(sim2,10,data=InitData)
  expect_that(Data2[1:4,],equals(InitData,tolerance=0.00000000000001))
	expect_that(Data2,equals(Data,tolerance=0.00000000000001))

  sim1 <- sim.vam(S&T&U ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
  sim2 <- sim.vam( S&T&U ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
  set.seed(0.5)
  Data<-simulate(sim1, 10)
  set.seed(0.5)
  InitData<-simulate(sim2, 4)
  Data2<-simulate(sim2,10,data=InitData)
  expect_that(Data2[1:4,],equals(InitData,tolerance=0.00000000000001))
  expect_that(Data2,equals(Data,tolerance=0.00000000000001))

  sim1 <- sim.vam(T&U ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
  sim2 <- sim.vam( T&U ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
  set.seed(0.5)
  Data<-simulate(sim1, 10,nb.system=1)
  set.seed(0.5)
  InitData<-simulate(sim2, 4)
  Data2<-simulate(sim2,10,nb.system=1,data=InitData)
  expect_that(Data2[1:4,],equals(InitData,tolerance=0.00000000000001))
  expect_that(Data2,equals(Data,tolerance=0.00000000000001))

  sim1 <- sim.vam( ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
  sim2 <- sim.vam( ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
  set.seed(0.5)
  Data<-simulate(sim1, Time>10)
  set.seed(0.5)
  InitData<-simulate(sim2, 4)
  Data2<-simulate(sim2,Time>10,data=InitData)
  expect_that(Data2[1:4,],equals(InitData,tolerance=0.00000000000001))
  expect_that(Data2,equals(Data,tolerance=0.00000000000001))

  sim1 <- sim.vam( ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
  sim2 <- sim.vam( ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
  set.seed(0.5)
  Data<-simulate(sim1, Time>(RC=10))
  set.seed(0.5)
  InitData<-simulate(sim2, 4)
  Data2<-simulate(sim2,Time>(RC=10),data=InitData)
  expect_that(Data2[1:4,],equals(InitData,tolerance=0.00000000000001))
  expect_that(Data2,equals(Data,tolerance=0.00000000000001))

	sim1 <- sim.vam( ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
	sim2 <- sim.vam( ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
	set.seed(0.5)
	Data<-simulate(sim1, Time>(RC=10))
	set.seed(0.5)
	InitData<-simulate(sim2, Time>(RC=3.8))
  set.seed(0.5)
  runif(length(InitData$Time)-1)
	Data2<-simulate(sim2,Time>(RC=10),data=InitData)
	expect_that(Data2,equals(Data,tolerance=0.00000000000001))

  sim1 <- sim.vam( ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
  sim2 <- sim.vam( ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
  set.seed(0.5)
  Data<-simulate(sim1, Time>3)
  set.seed(0.5)
  InitData<-simulate(sim2, 4)
  Data2<-simulate(sim2,Time>3,data=InitData)
  expect_that(Data2,equals(Data,tolerance=0.00000000000001))

  sim1 <- sim.vam( ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
  sim2 <- sim.vam( ~ (ARAInf(0.6) | Weibull(0.3,2.4)))
  set.seed(0.5)
  Data<-simulate(sim1, 3)
  set.seed(0.5)
  InitData<-simulate(sim2, 4)
  Data2<-simulate(sim2,3,data=InitData)
  expect_that(Data2,equals(Data,tolerance=0.00000000000001))
}
)


test_that("One system CM and PM",{
	sim1 <- sim.vam( ~ (ARA1(0.6) | Weibull(0.3,2.4))&(AGAN()|AtFailureProbability(0.6)))
  sim2 <- sim.vam( ~ (ARA1(0.6) | Weibull(0.3,2.4))&(AGAN()|AtFailureProbability(0.6)))
	set.seed(0.5)
	Data<-simulate(sim1, 15)
	set.seed(0.5)
  InitData<-simulate(sim2, 7)
  Data2<-simulate(sim2,15,data=InitData)
  expect_that(Data2[1:7,],equals(InitData,tolerance=0.00000000000001))
	expect_that(Data2,equals(Data,tolerance=0.00000000000001))

  sim1 <- sim.vam(S&T ~ (ARA1(0.6) | Weibull(0.3,2.4))&(AGAN()|AtFailureProbability(0.6)))
  sim2 <- sim.vam(S&T ~ (ARA1(0.6) | Weibull(0.3,2.4))&(AGAN()|AtFailureProbability(0.6)))
  set.seed(0.5)
  Data<-simulate(sim1, 15)
  set.seed(0.5)
  InitData<-simulate(sim2, 7)
  Data2<-simulate(sim2,15,data=InitData)
  expect_that(Data2[1:7,],equals(InitData,tolerance=0.00000000000001))
  expect_that(Data2,equals(Data,tolerance=0.00000000000001))
}
)

test_that("Several systems only CM",{
	sim1 <- sim.vam( ~ (AGAN() | Weibull(0.3,2.4)))
  sim2 <- sim.vam( ~ (AGAN() | Weibull(0.3,2.4)))
  n1_init<-5
	n_tot<-15
	n2_init<-6
	set.seed(0.5)
	Data_1<-simulate(sim1, n_tot)
	set.seed(0.5)
	runif(n_tot-n2_init)
	Data_2<-simulate(sim1, n_tot)
	set.seed(0.5)
  InitData<-rbind(cbind(System=1,Data_1[1:n1_init,]),cbind(System=2,Data_2[1:n2_init,]))
	Data<-rbind(cbind(System=1,Data_1),cbind(System=2,Data_2))
	runif(n1_init)
  Data2<-simulate(sim2,n_tot,data=InitData,nb.syst=2)
	expect_that(Data2,equals(Data,tolerance=0.00000000000001))

	sim1 <- sim.vam( ~ (AGAN() | Weibull(0.3,2.4)))
	sim2 <- sim.vam( ~ (AGAN() | Weibull(0.3,2.4)))
	n1_init<-5
	n_tot<-15
	n2_init<-6
	n3_init<-3
	set.seed(0.5)
	Data_1<-simulate(sim1, n_tot)
	set.seed(0.5)
	runif(n_tot-n2_init)
	Data_2<-simulate(sim1, n_tot)
	set.seed(0.5)
	runif(2*n_tot-n2_init-n3_init)
	Data_3<-simulate(sim1, n_tot)
	set.seed(0.5)
	InitData<-rbind(cbind(System=1,Data_1[1:n1_init,]),cbind(System=2,Data_2[1:n2_init,]),cbind(System=3,Data_3[1:n3_init,]))
	Data<-rbind(cbind(System=1,Data_1),cbind(System=2,Data_2),cbind(System=3,Data_3))
	runif(n1_init)
	Data2<-simulate(sim2,n_tot,data=InitData,nb.syst=3)
	expect_that(Data2,equals(Data,tolerance=0.00000000000001))

	sim1 <- sim.vam( ~ (AGAN() | Weibull(0.3,2.4)))
	sim2 <- sim.vam( ~ (AGAN() | Weibull(0.3,2.4)))
	n1_init<-5
	n_tot<-15
	n2_init<-0
	n3_init<-3
	set.seed(0.5)
	Data_1<-simulate(sim1, n_tot)
	set.seed(0.5)
	runif(n_tot-n2_init)
	Data_2<-simulate(sim1, n_tot)
	set.seed(0.5)
	runif(2*n_tot-n2_init-n3_init)
	Data_3<-simulate(sim1, n_tot)
	set.seed(0.5)
	InitData<-rbind(cbind(System=1,Data_1[1:n1_init,]),cbind(System=3,Data_3[1:n3_init,]))
	Data<-rbind(cbind(System=1,Data_1),cbind(System=2,Data_2),cbind(System=3,Data_3))
	runif(n1_init)
	Data2<-simulate(sim2,n_tot,data=InitData,nb.syst=3)
	expect_that(Data2,equals(Data,tolerance=0.00000000000001))

	sim1 <- sim.vam(T&U ~ (AGAN() | Weibull(0.3,2.4)))
	sim2 <- sim.vam(T&U ~ (AGAN() | Weibull(0.3,2.4)))
	n1_init<-5
	n_tot<-15
	n2_init<-6
	set.seed(0.5)
	Data_1<-simulate(sim1, n_tot)
	set.seed(0.5)
	runif(n_tot-n2_init)
	Data_2<-simulate(sim1, n_tot)
	set.seed(0.5)
	InitData<-rbind(cbind(System=1,Data_1[1:n1_init,]),cbind(System=2,Data_2[1:n2_init,]))
	Data<-rbind(cbind(System=1,Data_1),cbind(System=2,Data_2))
	runif(n1_init)
	Data2<-simulate(sim2,n_tot,data=InitData,nb.syst=2)
	expect_that(Data2,equals(Data,tolerance=0.00000000000001))

	sim1 <- sim.vam(S&T&U ~ (AGAN() | Weibull(0.3,2.4)))
	sim2 <- sim.vam(S&T&U ~ (AGAN() | Weibull(0.3,2.4)))
	n1_init<-5
	n_tot<-15
	n2_init<-6
	set.seed(0.5)
	Data_1<-simulate(sim1, n_tot)
	set.seed(0.5)
	runif(n_tot-n2_init)
	Data_2<-simulate(sim1, n_tot)
	set.seed(0.5)
	InitData<-rbind(cbind(S=1,Data_1[1:n1_init,]),cbind(S=2,Data_2[1:n2_init,]))
	Data<-rbind(cbind(S=1,Data_1),cbind(S=2,Data_2))
	runif(n1_init)
	Data2<-simulate(sim2,n_tot,data=InitData,nb.syst=2)
	expect_that(Data2,equals(Data,tolerance=0.00000000000001))
}
)

test_that("Several systems CM and PM",{
	sim1 <- sim.vam( ~ (AGAN() | Weibull(0.3,2.4))&(AGAN()|AtFailureProbability(0.6)))
  sim2 <- sim.vam( ~ (AGAN() | Weibull(0.3,2.4))&(AGAN()|AtFailureProbability(0.6)))
  n1_init<-5
	n_tot<-15
	n2_init<-6
	set.seed(0.5)
	Data_1<-simulate(sim1, n_tot)
	set.seed(0.5)
	runif(n_tot-n2_init)
	Data_2<-simulate(sim1, n_tot)
	set.seed(0.5)
  InitData<-rbind(cbind(System=1,Data_1[1:n1_init,]),cbind(System=2,Data_2[1:n2_init,]))
	Data<-rbind(cbind(System=1,Data_1),cbind(System=2,Data_2))
	runif(n1_init)
  Data2<-simulate(sim2,n_tot,data=InitData,nb.syst=2)
	expect_that(Data2,equals(Data,tolerance=0.00000000000001))
}
)

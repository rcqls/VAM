require(VAM)

set.seed(10)
n<-1000
sim <- sim.vam( ~ (ARA1(0.8)|Weibull(0.005,3))) 
df <- simulate(sim,100,nb.syst=10)

mle <- mle.vam(System&Time&Type~(ARA1(0.8)|Weibull(0.005,3)),data=df) 
coef(mle)
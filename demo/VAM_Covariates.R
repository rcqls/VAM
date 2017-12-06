require(VAM)

set.seed(100)
n<-100
df.cov <- data.frame(cov1=runif(n),cov2=runif(n))
sim <- sim.vam( ~ (ARA1(0.8)|Weibull(0.005,3|1*cov1-2*cov2)),data.covariates = df.cov) 
df <- simulate(sim)
df2 <- simulate(sim,5)

mle <- mle.vam(System&Time&Type~(ARA1(0.8)|Weibull(0.005,3|1*cov1-2*cov2)),data=df,data.covariates = df.cov) 
coef(mle)

mle2 <- mle.vam(System&Time&Type~(ARA1(0.8)|Weibull(0.005,3|1*cov1-2*cov2)),data=df2,data.covariates = df.cov) 
coef(mle2)
require(VAM)
#gctorture(TRUE)
#Rprofmem("boot.memprof",threshold=1000)
set.seed(10)
n<-10
df.cov <- data.frame(cov1=runif(n),cov2=runif(n))
sim <- sim.vam( ~ (ARA1(0.8)|Weibull(0.005,3|1*cov1+2*cov2)),data.covariates = df.cov) 
df <- simulate(sim)

mle <- mle.vam(System&Time&Type~(ARA1(0.8)|Weibull(0.005,3|1*cov1+2*cov2)),data=df,data.covariates = df.cov) 
coef(mle)

mle.0 <- mle.vam(System&Time&Type~(ARA1(0.8)|Weibull(0.005,3|0*cov1+0*cov2)),data=df,data.covariates = df.cov) 
coef(mle.0)

df2 <- simulate(sim,50)

 mle2 <- mle.vam(System&Time&Type~(ARA1(0.8)|Weibull(0.005,3|1*cov1+2*cov2)),data=df2,data.covariates = df.cov) 
 coef(mle2)

 mle2.0 <- mle.vam(System&Time&Type~(ARA1(0.8)|Weibull(0.005,3|0*cov1+0*cov2)),data=df2,data.covariates = df.cov) 
 coef(mle2.0)
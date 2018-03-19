require(VAM)
require(pryr)
gcinfo(TRUE)
#gctorture(TRUE)
#Rprofmem("boot.memprof",threshold=1000)
set.seed(10)
n<-10
df.cov <- data.frame(cov1=runif(n),cov2=runif(n))
sim <- sim.vam( ~ (ARA1(0.8)|Weibull(0.005,3|1*cov1+2*cov2)),data.covariates = df.cov) 
df <- simulate(sim)
print("ici0AV");gc(TRUE); print("ici0AP")
mle <- mle.vam(System&Time&Type~(ARA1(0.8)|Weibull(0.005,3|1*cov1+2*cov2)),data=df,data.covariates = df.cov) 
print("ici1AV");gc(TRUE); print("ici1AP")
#print(mem_change(
	coef(mle)
#))
print("ici2AV");gc(TRUE); print("ici2AP")


print("ici00AV");gc(TRUE); print("ici00AP")
mle <- mle.vam(System&Time&Type~(ARA1(0.8)|Weibull(0.005,3|1*cov1+2*cov2)),data=df,data.covariates = df.cov) 
print("ici11AV");gc(TRUE); print("ici11AP")
#print(mem_change(
	coef(mle)
#))
print("ici22AV");gc(TRUE); print("ici22AP")

# mle.0 <- mle.vam(System&Time&Type~(ARA1(0.8)|Weibull(0.005,3|0*cov1+0*cov2)),data=df,data.covariates = df.cov) 
#  gc(TRUE); print("ici")
# coef(mle.0)
#  gc(TRUE); print("ici")

# df2 <- simulate(sim,50)
#  gc(TRUE); print("ici")

#  mle2 <- mle.vam(System&Time&Type~(ARA1(0.8)|Weibull(0.005,3|1*cov1+2*cov2)),data=df2,data.covariates = df.cov) 
#   gc(TRUE); print("ici")
#  coef(mle2)
#  gc(TRUE); print("ici")

#  mle2.0 <- mle.vam(System&Time&Type~(ARA1(0.8)|Weibull(0.005,3|0*cov1+0*cov2)),data=df2,data.covariates = df.cov) 
#   gc(TRUE); print("ici")
#  coef(mle2.0)
#  gc(TRUE); print("ici")
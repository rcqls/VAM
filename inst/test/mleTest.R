data.frame(Time=c(0,2.5,3.8,4.9,5.35),Type=rep(-1,5)) -> tmp
tmp1<-tmp[1:4,]
tmp2<-tmp[2:5,]
require(VAM)

mleCpp <- mle.vam(Time & Type ~ (ARA1(.4) | Weibull(.001,2.5)) & (ARA1(.7) + ARA1(.7) | Periodic(1000,prob=c(.5,.5)) ),data=tmp2)
print(mleCpp$rcpp$contrast(c(1,1.75,.65,.7,.7)->par0))

mle <- mle.vam(Time & Type ~ (ARA1(.4) | Weibull(.001,2.5)) & (ARA1(.7) + ARA1(.7) | Periodic(1000,prob=c(.5,.5)) ),data=tmp2)
print(contrast.mle.vam(mle,par0))

beta<-par0[2]
rho <- par0[3]
sum(log(beta*(tmp2$Time-rho*tmp1$Time)^(beta-1)))->part1
sum((tmp2$Time-rho*tmp1$Time)^beta) -> part2
sum(((1-rho)*tmp1$Time)^beta)->part3
print(part1 - length(tmp1$Time)*log(part2-part3))

#script used to develop tests for logLik and contrast in the presence of covariates
#The considered test name
nbtest<-"WA8cov"

switch(nbtest,
       WA8cov={
         #Weibull+GQR_ARAm+PM_ARA1+PM_GQR_ARAm+2cov
         DataC<-data.frame(System=c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4),Time=c(0.79,1.583,1.761,2,2.652,2.705,3.445,4,4.57,4.901,0.069,0.477,0.682,1.004,1.327,1.704,1.83,1.848,0.735,1.004,1.408,1.874,2,2.466,2.514,2.919,3.376,1.02,1.401,2,2.216,2.752,3.426,3.592,4,4.35,4.82),Type=c(-1,-1,-1,2,1,-1,1,2,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,-1,1,2,-1,-1,1,1,-1,-1,2,-1,1,1,-1,2,-1,1))
         Cov<-data.frame(cov1=c(4.336,5.615,4.770,4.655),cov2=c(1,0,0,1))
         mle<-mle.vam(System&Time&Type~(GQR_ARAm(1.1,0.6|sqrt,3)|Weibull(0.1,2.1|0.6*cov1+(-0.7)*cov2))&(ARA1(0.9)+GQR_ARAm(1.3,0.8|log,2)),data=DataC,data.covariates = Cov)
         theta<-c(0.11,2.5,1.15,0.62,0.92,1.25,0.75,0.65,-0.56)
         
         n<-4
         mle2<-list()
         length(mle2)<-n
         for(i in 1:n){
           a<-0.11*exp(0.65*Cov$cov1[i]+(-0.56)*Cov$cov2[i])
           mle2[[i]]<-mle.vam(System&Time&Type~(GQR_ARAm(1.15,0.62|sqrt,3)|Weibull(a,2.5))&(ARA1(0.92)+GQR_ARAm(1.25,0.75|log,2)),data=DataC[DataC$System==i,])
         }
         Lcalc<-sum(unlist(lapply(mle2,logLik.mle.vam)))
         
       },
       WA7cov={
         #Weibull+ABAO+PM_AGAN+2cov
         DataC<-data.frame(System=c(rep(1,5),rep(2,4),rep(3,1),rep(4,5)),Time=c(0.800,2.646,3.190,3.916,4.109,0.910,1.127,1.245, 1.349, 0.541,1.397,1.463,2.406,2.506,3.159),Type=c(1,-1,1,1,0, -1,-1,1,-1, 0, 1,-1,-1,1,1))
         Cov<-data.frame(cov1=c(4.336,5.615,4.770,4.655),cov2=c(1,0,0,1))
         mle<-mle.vam(System&Time&Type~(ABAO()|Weibull(0.1,2.1|0.6*cov1+0.7*cov2))&(AGAN()),data=DataC,data.covariates = Cov)
         theta<-c(0.1,2.5,0.6,-0.9)
         
         n<-4
         mle2<-list()
         length(mle2)<-n
         for(i in 1:n){
           a<-0.1*exp(0.6*Cov$cov1[i]-0.9*Cov$cov2[i])
           mle2[[i]]<-mle.vam(System&Time&Type~(ABAO()|Weibull(a,2.5))&(AGAN()),data=DataC[DataC$System==i,])
         }
         Lcalc<-sum(unlist(lapply(mle2,logLik.mle.vam)))
         
       },
       WA6cov={
         #Weibull+GQRARA3+PM_AGAN+2cov
         DataC<-data.frame(System=c(rep(1,5),rep(2,4),rep(3,1),rep(4,5)),Time=c(0.800,2.646,3.190,3.916,4.109,0.910,1.127,1.245, 1.349, 0.541,1.397,1.463,2.406,2.506,3.159),Type=c(1,-1,1,1,0, -1,-1,1,-1, 0, 1,-1,-1,1,1))
         Cov<-data.frame(cov1=c(4.336,5.615,4.770,4.655),cov2=c(1,0,0,1))
         mle<-mle.vam(System&Time&Type~(GQR_ARAm(0.9,0.6|log,3)|Weibull(0.1,2.1|0.6*cov1+0.7*cov2))&(AGAN()),data=DataC,data.covariates = Cov)
         theta<-c(0.1,2.5,0.9,0.3,0.6,-0.9)
         
         n<-4
         mle2<-list()
         length(mle2)<-n
         for(i in 1:n){
           a<-0.1*exp(0.6*Cov$cov1[i]-0.9*Cov$cov2[i])
           mle2[[i]]<-mle.vam(System&Time&Type~(GQR_ARAm(0.9,0.3|log,3)|Weibull(a,2.5))&(AGAN()),data=DataC[DataC$System==i,])
         }
         Lcalc<-sum(unlist(lapply(mle2,logLik.mle.vam)))
         
       },
       WA5cov={
         #Weibull+GQRARA3+PM_ARA2+2cov
         DataC<-data.frame(System=c(rep(1,5),rep(2,4),rep(3,1),rep(4,5)),Time=c(0.800,2.646,3.190,3.916,4.109,0.910,1.127,1.245, 1.349, 0.541,1.397,1.463,2.406,2.506,3.159),Type=c(1,-1,1,1,0, -1,-1,1,-1, 0, 1,-1,-1,1,1))
         Cov<-data.frame(cov1=c(4.336,5.615,4.770,4.655),cov2=c(1,0,0,1))
         mle<-mle.vam(System&Time&Type~(GQR_ARAm(0.9,0.6|log,3)|Weibull(0.1,2.1|0.6*cov1+0.7*cov2))&(ARAm(0.5|2)),data=DataC,data.covariates = Cov)
         theta<-c(0.1,2.5,0.9,0.3,0.7,0.6,-0.9)
         
         n<-4
         mle2<-list()
         length(mle2)<-n
         for(i in 1:n){
           a<-0.1*exp(0.6*Cov$cov1[i]-0.9*Cov$cov2[i])
           mle2[[i]]<-mle.vam(System&Time&Type~(GQR_ARAm(0.9,0.3|log,3)|Weibull(a,2.5))&(ARAm(0.7|2)),data=DataC[DataC$System==i,])
         }
         Lcalc<-sum(unlist(lapply(mle2,logLik.mle.vam)))
         
       },
       WA2cov={
         #Weibull+ARAInf+2cov
         DataC<-data.frame(System=c(rep(1,5),rep(2,4),rep(3,1),rep(4,5)),Time=c(0.800,2.646,3.190,3.916,4.109,0.910,1.127,1.245, 1.349, 0.541,1.397,1.463,2.406,2.506,3.159),Type=rep(-1,15))
         Cov<-data.frame(cov1=c(4.336,5.615,4.770,4.655),cov2=c(0,0,0,1))
         mle<-mle.vam(System&Time&Type~(ARAInf(0.8)|Weibull(0.15,2.3|0.6*cov1-0.9*cov2)),data=DataC,data.covariates = Cov)
         theta<-c(0.15,2.3,0.8,0.6,-0.9)
         
         n<-4
         mle2<-list()
         length(mle2)<-n
         for(i in 1:n){
           a<-0.15*exp(0.6*Cov$cov1[i]-0.9*Cov$cov2[i])
           mle2[[i]]<-mle.vam(System&Time&Type~(ARAInf(0.8)|Weibull(a,2.3)),data=DataC[DataC$System==i,])
         }
         Lcalc<-sum(unlist(lapply(mle2,logLik.mle.vam)))
                    
       },
       WA1cov={
         #Weibull+ARAInf+cov
         DataC<-data.frame(System=c(rep(1,5),rep(2,4),rep(3,1),rep(4,5)),Time=c(0.800,2.646,3.190,3.916,4.109,0.910,1.127,1.245, 1.349, 0.541,1.397,1.463,2.406,2.506,3.159),Type=rep(-1,15))
         Cov<-data.frame(cov1=c(4.336,5.615,4.770,4.655))
         mle<-mle.vam(System&Time&Type~(ARAInf(0.8)|Weibull(0.15,2.3|0.6*cov1)),data=DataC,data.covariates = Cov)
         theta<-c(0.15,2.3,0.8,0.6)
         
         n<-4
         mle2<-list()
         length(mle2)<-n
         for(i in 1:n){
           a<-0.15*exp(0.6*Cov$cov1[i])
           mle2[[i]]<-mle.vam(System&Time&Type~(ARAInf(0.8)|Weibull(a,2.3)),data=DataC[DataC$System==i,])
         }
         Lcalc<-sum(unlist(lapply(mle2,logLik.mle.vam)))
         
       },
       WA3cov={
         #Weibull+ARA1+cov+1syst
         DataC<-data.frame(System=c(rep(1,5)),Time=c(0.800,2.646,3.190,3.916,4.109),Type=rep(-1,5))
         Cov<-data.frame(cov1=c(4.336))
         mle<-mle.vam(System&Time&Type~(ARA1(0.8)|Weibull(0.15,2.3|0.6*cov1)),data=DataC,data.covariates = Cov)
         theta<-c(0.15,2.3,0.8,0.6)
         
         n<-1
         mle2<-list()
         length(mle2)<-n
         for(i in 1:n){
           a<-0.15*exp(0.6*Cov$cov1[i])
           mle2[[i]]<-mle.vam(System&Time&Type~(ARA1(0.8)|Weibull(a,2.3)),data=DataC[DataC$System==i,])
         }
         Lcalc<-sum(unlist(lapply(mle2,logLik.mle.vam)))
         
       },
       WA4cov={
         #loglinear+GQRARA3+1cov
         DataC<-data.frame(System=c(rep(1,5),rep(2,4),rep(3,1),rep(4,5)),Time=c(0.800,2.646,3.190,3.916,4.109,0.910,1.127,1.245, 1.349, 0.541,1.397,1.463,2.406,2.506,3.159),Type=rep(-1,15))
         Cov<-data.frame(cov1=c(4.336,5.615,4.770,4.655))
         mle<-mle.vam(System&Time&Type~(GQR_ARAm(0.9,0.6|log,3)|LogLinear(0.01,0.8|0.6*cov1)),data=DataC,data.covariates = Cov)
         theta<-c(10,0.8,0.9,0.6,0.05)
         
         n<-4
         mle2<-list()
         length(mle2)<-n
         for(i in 1:n){
           a<-10*exp(0.05*Cov$cov1[i])
           mle2[[i]]<-mle.vam(System&Time&Type~(GQR_ARAm(0.9,0.6|log,3)|LogLinear(a,0.8)),data=DataC[DataC$System==i,])
         }
         Lcalc<-sum(unlist(lapply(mle2,logLik.mle.vam)))
         
       },
)



#########
L<-logLik(mle,theta,c(TRUE,FALSE,FALSE))
dL<-logLik(mle,theta,c(FALSE,TRUE,FALSE))
d2L<-logLik(mle,theta,c(FALSE,FALSE,TRUE))
C<-contrast(mle,theta,c(TRUE,FALSE,FALSE))
dC<-contrast(mle,theta,c(FALSE,TRUE,FALSE))
d2C<-contrast(mle,theta,c(FALSE,FALSE,TRUE))


fix<-rep(TRUE,length(theta))
fix[1]=FALSE
run(mle,fixed=fix,verbose=FALSE);alpha_Est<-coef(mle)[1]
Ccalc<-contrast(mle,c(alpha_Est,theta[2:length(theta)]))

epsilon<-1e-7
EstdL<-dL
EstdC<-dC
for(i in (1:length(dL))){
  theta1<-theta
  theta1[i]<-theta1[i]+epsilon
  EstdL[i]<-(logLik(mle,theta1,c(TRUE,FALSE,FALSE))-logLik(mle,theta,c(TRUE,FALSE,FALSE)))/epsilon
}
for(i in (1:length(dC))){
  theta1<-theta
  theta1[i+1]<-theta1[i+1]+epsilon
  EstdC[i]<-(contrast(mle,theta1,c(TRUE,FALSE,FALSE))-contrast(mle,theta,c(TRUE,FALSE,FALSE)))/epsilon
}
Estd2L<-d2L
Estd2C<-d2C
for(i in (1:(dim(d2L)[1]))){
  theta1<-theta
  theta1[i]<-theta1[i]+epsilon
  Estd2L[i,]<-(logLik(mle,theta1,c(FALSE,TRUE,FALSE))-logLik(mle,theta,c(FALSE,TRUE,FALSE)))/epsilon
}
if(!(is.null(dim(d2C)))){
  for(i in (1:(dim(d2C)[1]))){
    theta1<-theta
    theta1[i+1]<-theta1[i+1]+epsilon
    Estd2C[i,]<-(contrast(mle,theta1,c(FALSE,TRUE,FALSE))-contrast(mle,theta,c(FALSE,TRUE,FALSE)))/epsilon
  }
} else {
  theta1<-theta
  theta1[2]<-theta1[2]+epsilon
  Estd2C<-(contrast(mle,theta1,c(FALSE,TRUE,FALSE))-contrast(mle,theta,c(FALSE,TRUE,FALSE)))/epsilon
}


ecrit<-function(aprinter,nl,nc){
  if((nl>0)&&(nc>0)){
    for(i in 1:nc){
      for(j in 1:nl){
        if((i==1)&&(j==1)) res<-paste("matrix(c(",aprinter[1,1],sep="")
        else res<-paste(res,aprinter[j,i],sep=',')
      }
    }
    res<-paste(res,"),nrow=",nl,',byrow=TRUE)',sep="")
  }
  else if((nl>0)||(nc>0)){
    for(j in 1:max(nl,nc)){
      if(j==1) res<-paste("c(",aprinter[1],sep="")
      else res<-paste(res,aprinter[j],sep=',')
    }
    res<-paste(res,")",sep="")
  }
  else res<-paste(aprinter)
  return(res)
  
}
print(max(abs((L-Lcalc)),abs((C-Ccalc)),abs((dL-EstdL)),abs((dC-EstdC)),abs((d2L-Estd2L)),abs((d2C-Estd2C))))
print(max(abs((L-Lcalc)/L),abs((C-Ccalc)/C),abs((dL-EstdL)/dL),abs((dC-EstdC)/dC),abs((d2L-Estd2L)/d2L),abs((d2C-Estd2C)/d2C)))

res<-paste(" lnL=",ecrit(L,0,0),"\n",sep="")
res<-paste(res," dlnL=",ecrit(dL,0,length(dL)),"\n",sep="")
res<-paste(res," d2lnL=",ecrit(d2L,dim(d2L)[1],dim(d2L)[2]),"\n",sep="")

res<-paste(res," C=",ecrit(C,0,0),"\n",sep="")
res<-paste(res," dC=",ecrit(dC,0,length(dC)),"\n",sep="")
if(!(is.null(dim(d2C)))){
  res<-paste(res," d2C=",ecrit(d2C,dim(d2C)[1],dim(d2C)[2]),"\n",sep="")
} else {
  res<-paste(res," d2C=",ecrit(d2C,0,0),"\n",sep="")
}
res<-paste(res,"expect_that(logLik.mle.vam(mle,par0=theta) ,equals(lnL,tolerance=0.00000000000001))\n expect_that(contrast.mle.vam(mle,par0=theta) ,equals(C,tolerance=0.00000000000001))\n expect_that(logLik.mle.vam(mle,par0=theta,with_value=FALSE,with_gradient=TRUE) ,equals(dlnL,tolerance=0.00000000000001))\n expect_that(contrast(mle,par0=theta,with_value=FALSE,with_gradient=TRUE) ,equals(dC,tolerance=0.00000000000001))\n expect_that(logLik.mle.vam(mle,par0=theta,with_value=FALSE,with_hessian=TRUE) ,equals(d2lnL,tolerance=0.00000000000001))\n expect_that(contrast(mle,par0=theta,with_value=FALSE,with_hessian=TRUE) ,equals(d2C,tolerance=0.00000000000001))")

message(res)

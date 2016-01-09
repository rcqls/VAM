require(VAM)
##modMP <- model.vam( ~ (ARA1(.8) | Weibull(.0012,2.4)) & (ARA1(.7) + ARA1(.6) + ARA1(.5)))
modMP <- model.vam( ~ (ARA1(.87) | Weibull(.0015,2.7)) & (ARAInf(.44)))
# Simulated model
## form0 <- ~ (ARA1(.4) | Weibull(.001,2.5)) & (ARA1(.7) + ARA1(.6) + ARA1(.5) | AtIntensity(1.5,modMP) )
form0 <- ~ (ARA1(.9) | Weibull(.001,2.5)) & (ARAInf(.4) | AtIntensity(0.2,modMP))
print(parse.vam.formula(NULL,form0) ->model0)

sim0<-sim.vam( form0 )

(tmp<-simulate(sim0,100))

sim0bis <- sim.vam( ~ (ARA1(.9) | Weibull(.001,2.5)) & (ARAInf(.4) | AtIntensity(0.2)) )

(tmp2<-simulate(sim0bis,100))

#(tmp3<-simulate(sim0,TimeGreaterThanCensorship(1000) | SizeGreaterThan(1000)))
#(tmp3<-simulate(sim0,Time > (RightCensorship=1000) & Size > 1000))
(tmp3<-simulate(sim0, T > (RC=1000) | S > 1000))

#(tmp3<-simulate(sim0,EndAt(time=1000),nb.system=3))
(tmp3bis<-simulate(sim0,Time > (RightCensorship=1000),nb.system=3))

#(tmp3<-simulate(sim0,Time > (RightCensorship=Unif(50,100))))
#(tmp3<-simulate(sim0,T > (RC=Unif(50,100))))

#(tmp4<-simulate(sim0,SizeOfTypeGreaterThan(size=10,type=1) & SizeOfTypeGreaterThan(size=10,type=-1)))
#(tmp3<-simulate(sim0,Size[Type=1]>10 & Size[Type=-1]>10))
(tmp4<-simulate(sim0,S[1]>10 & S[-1]>10))


context("Contrast and log-likelihood")

test_that("Weibull",{
	simData<-data.frame(Time=c(3.36),Type=c(-1),row.names=1:1)
	mle <- mle.vam(Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)),data=simData)
	theta<-c(0.3,1.8,0.6)
	lnL<- -2.30449245951301
	dlnL<- c(-5.52619699756427,-1.45367181592636,0)
	d2lnL<- matrix(c(-11.1111111111111,-10.7372278181901,0,-10.7372278181901,-4.21250787723964,0,0,0,0),nrow=3,byrow=TRUE)
	C<- -1.62415430907299
	dC<-c(0.555555555555556,0)
	d2C<-matrix(c(-0.308641975308642,0,0,0),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("LogLinear",{
	simData<-data.frame(Time=c(3.36),Type=c(-1),row.names=1:1)
	mle <- mle.vam(Time & Type ~ (ARAInf(0.4) | LogLinear(0.001,2.5)),data=simData)
	theta<-c(0.3,0.8,0.6)
	lnL<- -3.65431355894635
	dlnL<- c(-13.7944691820681,-8.74189899224908,0)
	d2lnL<- matrix(c(-11.1111111111111,-40.3396633074969,0,-40.3396633074969,-31.9886643027399,0,0,0,0),nrow=3,byrow=TRUE)
	C<- -1.15270302129339
	dC<-c(1.00478465516967,0)
	d2C<-matrix(c(-0.678445876029819,0,0,0),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + MC ARAInf",{
	simData<-data.frame(Time=c(3.36,4.04,4.97,5.16),Type=c(-1,-1,-1,-1),row.names=1:4)
	mle <- mle.vam(Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)),data=simData)
	theta<-c(0.3,0.8,0.6)
	lnL<- -7.37830963135462
	dlnL<- c(9.33348796771948,5.77076155284033,1.19923836457015)
	d2lnL<- matrix(c(-44.4444444444444,-5.26230480903023,-0.723044448779292,-5.26230480903023,-7.7471885008781,-6.30726171420585,-0.723044448779292,-6.30726171420585,0.435684125802398),nrow=3,byrow=TRUE)
	C<- -5.36231016699152
	dC<-c(2.08694474533596,0.693079297460398)
	d2C<-matrix(c(-4.31732300351524,-3.55104259186756,-3.55104259186756,-0.351517905180765),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARA1",{
	simData<-data.frame(Time=c(3.36,4.04,4.97,5.16),Type=c(-1,-1,-1,-1),row.names=1:4)
	mle <- mle.vam(Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5)),data=simData)
	theta<-c(0.3,0.8,0.6)
	lnL<- -7.60531410020218
	dlnL<- c(9.43046924122556,6.95299388770745,0.950641033424141)
	d2lnL<- matrix(c(-44.4444444444444,-5.58984454112956,-0.51705781427391,-5.58984454112956,-8.18607969766148,-5.04449078731948,-0.51705781427391,-5.04449078731948,1.70955451742894),nrow=3,byrow=TRUE)
	C<- -5.5202288755881
	dC<-c(2.90098061507367,0.575831839583211)
	d2C<-matrix(c(-4.65895384634265,-3.11529381117442,-3.11529381117442,0.845549840945982),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("LogLinear + CM ARAInf",{
	simData<-data.frame(Time=c(3.36,4.04,4.97,5.16),Type=c(-1,-1,-1,-1),row.names=1:4)
	mle <- mle.vam(Time & Type ~ (ARAInf(0.4) | LogLinear(0.001,2.5)),data=simData)
	theta<-c(0.3,0.8,0.6)
	lnL<- -5.23995692659076
	dlnL<- c( -9.43872569762339,-6.59383234684053,-3.31872264053507)
	d2lnL<- matrix(c(-44.4444444444444,-48.6775744894684,15.1103911982164,-48.6775744894684,-35.8815691840954,2.5634113080498,15.1103911982164,2.5634113080498,-4.2900710638037),nrow=3,byrow=TRUE)
	C<- -4.54940776890445
	dC<-c(-0.540963707156266,-5.19764150566532)
	d2C<-matrix(c(-2.73176918257204,-8.24081529909406,-8.24081529909406,4.60248500309959),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARA1 + Censorship",{
	simData<-data.frame(Time=c(3.36,4.04,4.97,5.16),Type=c(-1,-1,-1,0),row.names=1:4)
	mle <- mle.vam(Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5)),data=simData)
	theta<-c(0.3,0.8,0.6)
	lnL<- -6.02251633965988
	dlnL<- c(6.09713590789223,4.92458686319668,0.494259031587594)
	d2lnL<- matrix(c(-33.3333333333333,-5.58984454112956,-0.51705781427391,-5.58984454112956,-6.62357969766148,-2.76258077813674,-0.51705781427391,-2.76258077813674,0.66813185942727),nrow=3,byrow=TRUE)
	C<- -5.02903383164658
	dC<-c(2.30481524930606,0.25193147227744)
	d2C<-matrix(c(-3.49421538475699,-1.43675832133097,-1.43675832133097,0.116785670792289),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + multiSystems",{
	simData<-data.frame(System=c(1,2),Time=c(3.36,2.34),Type=c(-1,-1),row.names=1:2)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.3,0.8,0.6)
    logLik(mle,theta,TRUE,FALSE,FALSE)
    logLik(mle,theta,FALSE,TRUE,FALSE)
    logLik(mle,theta,FALSE,FALSE,TRUE)
    contrast(mle,theta,TRUE,FALSE,FALSE)
    contrast(mle,theta,FALSE,TRUE,FALSE)
    contrast(mle,theta,FALSE,FALSE,TRUE)
    theta<-c(0.3,1.8,0.6)
	lnL<- -3.62638734326211
	dlnL<- c(-6.81229529030084,-1.22612955802,0)
	d2lnL<- matrix(c(-22.2222222222222,-14.6644419082528,0,-14.6644419082528,-5.52276726509835,0,0,0,0),nrow=3,byrow=TRUE)
	C<- -2.9907189791188
	dC<-c(0.99730223588401,0)
	d2C<-matrix(c(-0.67625373876488,0,0,0),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARAInf + Censorship + multiSystems",{
	simData<-data.frame(System=c(1,1,1,1,2,2,2,3),Time=c(3.36,4.04,4.97,5.16,2.34,3.46,5.02,4),Type=c(-1,-1,-1,0,-1,-1,-1,0),row.names=1:8)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.3,0.8,0.6)
    logLik(mle,theta,TRUE,FALSE,FALSE)
    logLik(mle,theta,FALSE,TRUE,FALSE)
    logLik(mle,theta,FALSE,FALSE,TRUE)
    contrast(mle,theta,TRUE,FALSE,FALSE)
    contrast(mle,theta,FALSE,TRUE,FALSE)
    contrast(mle,theta,FALSE,FALSE,TRUE)
    theta<-c(0.3,1.8,0.6)
	lnL<- -10.7518360550876
	dlnL<- c(-16.6082469606206,-4.41424143581748,0.565291823677669)
	d2lnL<- matrix(c(-66.6666666666667,-42.1908945206754,17.8580795012604,-42.1908945206754,-16.0970944508072,5.49620639347231,17.8580795012604,5.49620639347231,-6.72216164111632),nrow=3,byrow=TRUE)
	C<- -9.39660957648529
	dC<-c(1.32804496573541,-1.86523751785326)
	d2C<-matrix(c(-1.66488790656225,-3.08811044563958,-3.08811044563958,-3.53294273656389),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARAInf + PM ARAInf + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (ARAInf(0.5)),data=simData)
    theta<-c(0.3,1.8,0.3,0.8)
	lnL<- -20.7592907965932
	dlnL<- c(-31.3904767850125,-9.15468712792062,6.37367814316912,3.41532811689743)
	d2lnL<- matrix(c(-100,-60.5303216502559,30.8177749831935,29.6544563881339,-60.5303216502559,-20.7750026445059,19.227056494657,8.24534883396676,30.8177749831935,19.227056494657,-6.01682216912879,-4.05851834664032,29.6544563881339,8.24534883396676,-4.05851834664032,-11.086407735654),nrow=4,byrow=TRUE)
	C<- -17.7866638232513
	dC<-c(0.130510297108522,1.64630982399109,-1.13359008787765)
	d2C<-matrix(c(-2.82300923948606,3.10568724954822,-3.76042913146761,3.10568724954822,-1.35443592315717,1.21888246564082,-3.76042913146761,1.21888246564082,-6.87069738841726),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARA1 + PM ARA1 + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5)) & (ARA1(0.5)),data=simData)
    theta<-c(0.3,1.8,0.3,0.8)
	lnL<- -22.803615512535
	dlnL<- c(-40.8686832560709,-14.5723320558868,5.97450010229359,6.11227548582455)
	d2lnL<- matrix(c(-100,-81.9188765158751,27.6301819586786,41.5815929285388,-81.9188765158751,-30.6353597460098,19.0183236155591,19.4470265766715,27.6301819586786,19.0183236155591,-0.112216256019809,-0.0688331609468675,41.5815929285388,19.4470265766715,-0.0688331609468675,-7.61532027419852),nrow=4,byrow=TRUE)
	C<- -18.2796917972255
	dC<-c(-0.399993316686276,1.19435278095851,-1.08152947541163)
	d2C<-matrix(c(-2.54493624659463,2.32631632432271,-2.45799033270139,2.32631632432271,0.400805866820642,1.41358355278611,-2.45799033270139,1.41358355278611,-7.64153503772091),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARA1 + PM ARAInf + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5)) & (ARAInf(0.5)),data=simData)
    theta<-c(0.3,1.8,0.3,0.8)
	lnL<- -21.4669497529826
	dlnL<- c(-34.4483795857384,-11.2730770047724,4.88510828172436,4.51171083061658)
	d2lnL<- matrix(c(-100,-68.4654207557133,23.5723992902919,35.4708804565318,-68.4654207557133,-25.4038060829933,15.7318257510918,11.6570381491589,23.5723992902919,15.7318257510918,-0.209594625288678,-1.22056348858044,35.4708804565318,11.6570381491589,-1.22056348858044,-13.2267740291002),nrow=4,byrow=TRUE)
	C<- -18.0144407241976
	dC<-c(-0.294417101217524,1.10519497567318,-1.17616384704568)
	d2C<-matrix(c(-3.15304348912119,2.36504515416587,-3.93130478045259,2.36504515416587,0.276480347185599,1.07610314939065,-3.93130478045259,1.07610314939065,-6.78594073027517),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARAInf + PM ARA1 + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (ARA1(0.5)),data=simData)
    theta<-c(0.3,1.8,0.3,0.8)
	lnL<- -21.8391429413397
	dlnL<- c(-36.4999918106382,-11.5827201058525,7.96675197477725,4.26686686122415)
	d2lnL<- matrix(c(-100,-70.5112747556172,38.1582103953299,32.8735045622425,-70.5112747556172,-24.0010530944472,24.0292065539748,13.1106219399306,38.1582103953299,24.0292065539748,-8.15842021737515,-4.8218313691535,32.8735045622425,13.1106219399306,-4.8218313691535,-6.78161720542848),nrow=4,byrow=TRUE)
	C<- -18.0531853838104
	dC<-c(0.0277818554092377,1.68355856820899,-1.14613822932782)
	d2C<-matrix(c(-2.23368136554917,2.97638451294678,-2.64164697486646,2.97638451294678,-1.56051218100989,1.37245448857941,-2.64164697486646,1.37245448857941,-6.95365901039246),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARAInf + PM ARAInf + PM ARA1 + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,10)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43),Type=c(1,2,2,1, -1,-1,-1,0, 1,-1,-1,2, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,0),row.names=1:26)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (ARAInf(0.5)+ARA1(-1)),data=simData)
    theta<-c(0.3,1.8,0.3,0.8,-1)
	lnL<- -28.5661244738475
	dlnL<- c(-54.3256561040677,-25.9832001219506,8.64422921282077,2.49039875914388,2.46171865807292)
	d2lnL<- matrix(c(-111.111111111111,-121.351826223399,38.9463329627851,31.2547997200674,9.05229027861169,-121.351826223399,-57.7998939468941,26.3496292986742,11.5491413087507,8.57176749488547,38.9463329627851,26.3496292986742,-6.6714102959566,-6.57594695903925,-0.289361831349399,31.2547997200674,11.5491413087507,-6.57594695903925,-9.72658096422371,-0.679084831246593,9.05229027861169,8.57176749488547,-0.289361831349399,-0.679084831246593,0.112472116953841),nrow=5,byrow=TRUE)
	C<- -21.9373903366026
	dC<-c(-3.4212782198478,1.40326589132925,-3.32054263577063,0.778702648695225)
	d2C<-matrix(c(-4.727246115034,1.51436129328023,-5.87868127434093,1.63317888623568,1.51436129328023,-1.31359453402465,0.78363997271008,0.315720104487723,-5.87868127434093,0.78363997271008,-7.84461370458548,0.065862948983942,1.63317888623568,0.315720104487723,0.065862948983942,0.0994433721114185),nrow=4,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + AGAN",{
	simData<-data.frame(Time=c(3.36,4.04,4.97,5.16),Type=c(-1,-1,-1,-1),row.names=1:4)
    mle <- mle.vam(Time & Type ~ (AGAN() | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.3,0.8)
	lnL<- -6.90098251895707
	dlnL<- c(8.75359407434948,3.3717767832448)
	d2lnL<- matrix(c(-44.4444444444444,-2.40399936753948,-2.40399936753948,-7.6652713149061),nrow=2,byrow=TRUE)
	C<- -5.25256034408912
	dC<-1.99329429144225
	d2C<--9.26821752088532
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("LogLinear + CM ABAO",{
    simData<-data.frame(Time=c(3.36,4.04,4.97,5.16),Type=c(-1,-1,-1,-1),row.names=1:4)
    mle <- mle.vam(Time & Type ~ (ABAO() | LogLinear(0.001,2.5)),data=simData)
    theta<-c(0.3,0.8)
	lnL<- -13.6870254780068
	dlnL<- c(-62.9837808690102,-73.9249749593491)
	d2lnL<- matrix(c(-44.4444444444444,-304.849916531164,-304.849916531164,-390.943849373404),nrow=2,byrow=TRUE)
	C<- -1.77041141396439
	dC<-1.55193690275558
	d2C<--4.47702275378921
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARAInf + PM AGAN + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (AGAN()),data=simData)
    theta<-c(0.3,1.8,0.3)
	lnL<- -20.373593780469
	dlnL<- c(-25.6662791768442,-8.48131134572237,5.63552763127959)
	d2lnL<- matrix(c(-100,-52.7375595999022,27.2348357573035,-52.7375595999022,-19.8545917834264,17.1554852531126,27.2348357573035,17.1554852531126,-5.23495386861318),nrow=3,byrow=TRUE)
	C<- -18.237304657921
	dC<-c(-1.18653460936235,1.86834447361728)
	d2C<-matrix(c(-3.90302067583932,3.61294880573467,3.61294880573467,-1.23797415150246),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARAInf + PM ABAO + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (ABAO()),data=simData)
    theta<-c(0.3,1.8,0.3)
	lnL<- -26.1289527701662
	dlnL<- c(-60.0989490737245,-26.356099335404,12.0916680482488)
	d2lnL<- matrix(c(-100,-131.380361146203,54.6108671722709,-131.380361146203,-59.0094678187948,40.3465131542417,54.6108671722709,40.3465131542417,-14.5659541685396),nrow=3,byrow=TRUE)
	C<- -17.9966681180183
	dC<-c(-0.06559856379452,1.16349710241685)
	d2C<-matrix(c(-2.36452544203233,1.90128283795096,1.90128283795096,-3.25686180452323),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM AGAN + PM ARAInf + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (AGAN() | Weibull(0.001,2.5)) & (ARAInf(0.5)),data=simData)
    theta<-c(0.3,1.8,0.6)
	lnL<- -18.3086928056166
	dlnL<- c(-17.4175161247657,-3.39893433955917,3.29451466868132)
	d2lnL<- matrix(c(-100,-34.0891031594782,20.1682479366685,-34.0891031594782,-11.4701363569356,8.57719280033067,20.1682479366685,8.57719280033067,-5.85778413970232),nrow=3,byrow=TRUE)
	C<- -17.2035868236485
	dC<-c(0.357573081011327,1.07203985966202)
	d2C<-matrix(c(-3.62570578781634,1.40918907729665,1.40918907729665,-3.47013127037913),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ABAO + PM ARAInf + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ABAO() | Weibull(0.001,2.5)) & (ARAInf(0.5)),data=simData)
    theta<-c(0.3,1.8,0.6)
	lnL<- -24.1615521488875
	dlnL<- c(-49.0982333806125,-19.7432656164282,7.2827784977942)
	d2lnL<- matrix(c(-100,-103.784190055221,41.227797314041,-103.784190055221,-44.1308967618595,22.3471560737357,41.227797314041,22.3471560737357,-11.2234514062957),nrow=3,byrow=TRUE)
	C<- -18.1575205644636
	dC<-c(-0.41684010691203,-0.394556106708563)
	d2C<-matrix(c(-2.96771319633972,-1.62522576816206,-1.62522576816206,-4.17301130756893),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARAInf + PM AGAN + PM ARA1 + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,10)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43),Type=c(1,2,2,1, -1,-1,-1,0, 1,-1,-1,2, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,0),row.names=1:26)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (AGAN()+ARA1(-1)),data=simData)
    theta<-c(0.3,1.8,0.3,-1)
	lnL<- -28.3428237855895
	dlnL<- c(-48.0459603807939,-24.5636656174879,7.49234651428247,2.31473950208313)
	d2lnL<- matrix(c(-111.111111111111,-109.700844422173,33.4242320339798,8.61601065450215,-109.700844422173,-53.6196344211062,22.5095841815981,8.03660251513632,33.4242320339798,22.5095841815981,-5.26135655808678,-0.180393079852264,8.61601065450215,8.03660251513632,-0.180393079852264,0.103529366551017),nrow=4,byrow=TRUE)
	C<- -22.8546653348801
	dC<-c(-5.13360355515621,1.57229269676842,0.788683594049934)
	d2C<-matrix(c(-5.61350145753183,1.812665200363,1.66531233950821,1.812665200363,-1.19934416903005,0.360960522359852,1.66531233950821,0.360960522359852,0.100675520443334),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull3",{
	simData<-data.frame(Time=c(3.36),Type=c(-1),row.names=1:1)
    mle <- mle.vam(Time & Type ~ (ARAInf(0.4) | Weibull3(0.001,2.5,5)),data=simData)
    theta<-c(0.3,1.8,4,0.6)
	lnL<- -6.28349650594271
	dlnL<- c(-20.8805277090384,-14.1662361334174,-0.920550636729322,0)
	d2lnL<- matrix(c(-11.1111111111111,-55.726172072379,-3.43082096301078,0,-55.726172072379,-36.7534932861936,-3.48854152770058,0,-3.43082096301078,-3.48854152770058,0.0228198059801746,0,0,0,0,0),nrow=4,byrow=TRUE)
	C<- -2.00229062844351
	dC<-c(0.250199288093325,-0.0329926542472398,0)
	d2C<-matrix(c(-0.0292037862530474,-0.0369910685383343,0,-0.0369910685383343,0.01048162449677,0,0,0,0),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull3 + CM ARAInf + PM AGAN + PM ARA1 + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,10)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43),Type=c(1,2,2,1, -1,-1,-1,0, 1,-1,-1,2, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,0),row.names=1:26)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull3(0.001,2.5,3)) & (AGAN()+ARA1(-1)),data=simData)
    theta<-c(0.3,1.8,4,0.3,-1)
	lnL<- -57.8500874222302
	dlnL<- c(-184.624600600652,-129.704071132207,-7.67058446328368,7.10281102985035,2.17137698160561)
	d2lnL<- matrix(c(-111.111111111111,-507.945231260056,-30.444281557097,26.9768005426141,7.49166354730297,-507.945231260056,-342.217351176241,-30.3516547909226,28.3089982845618,8.81234540697149,-30.444281557097,-30.3516547909226,0.061326762184357,-0.111853699101766,-0.0401271189392565,26.9768005426141,28.3089982845618,-0.111853699101766,-4.59887941278902,-0.184645955524138,7.49166354730297,8.81234540697149,-0.0401271189392565,-0.184645955524138,0.0830597880298306),nrow=5,byrow=TRUE)
	C<- -21.2401490744713
	dC<-c(-0.625238686313029,0.0659037960818534,0.247477614221491,0.267598534144922)
	d2C<-matrix(c(-0.640233922422266,0.162128309735713,0.396499507783223,0.466081163203727,0.162128309735713,-0.02785719138912,-0.0653908606723684,-0.042571225920877,0.396499507783223,-0.0653908606723684,-0.474210863914348,0.0143037570576406,0.466081163203727,-0.042571225920877,0.0143037570576406,0.0183816353566123),nrow=4,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + QR",{
	simData<-data.frame(Time=c(18.09,52.07,95.71,145.75),Type=c(-1,-1,-1,-1),row.names=1:4)
    mle <- mle.vam(Time & Type ~ (QR(0.7) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.03,2.4,0.7)
	lnL<- -161.716137729918
	dlnL<- c(-5393.65204881456,-489.310119943957,-789.802230021644)
	d2lnL<- matrix(c(-4444.44444444444,-16764.8888403069,-27012.4552864358,-16764.8888403069,-1529.01885118223,-2770.66341478734,-27012.4552864358,-2770.66341478734,-4691.475862947),nrow=3,byrow=TRUE)
	C<- -14.8047584598105
	dC<-c(1.5034291586224,1.02192487729127)
	d2C<-matrix(c(-0.760709275195016,0.824096797385021,0.824096797385021,-46.3104800299686),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + QR + MultiSystems + Censorship",{
	simData<-data.frame(System=c(1,1,1,1,2,2,2,3),Time=c(18.09,52.07,95.71,145.75,15.02,45.1,82,20.1),Type=c(-1,-1,-1,0,-1,-1,-1,-1),row.names=1:8)
    mle <- mle.vam(System & Time & Type ~ (QR(0.7) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.03,2.4,0.7)
	lnL<- -293.755998370921
	dlnL<- c(-9842.97490606146,-881.745942144996,-1158.18274977465)
	d2lnL<- matrix(c(-7777.77777777778,-30181.3946412058,-39291.8059448694,-30181.3946412058,-2718.32845215212,-4014.06792884068,-39291.8059448694,-4014.06792884068,-6162.26303978957),nrow=3,byrow=TRUE)
	C<- -24.8250704635939
	dC<-c(2.72891618219407,-6.72454482275228)
	d2C<-matrix(c(-1.33241877266666,-2.82017460407049,-2.82017460407049,-64.9658790434754),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + GQR + MultiSystems + Censorship",{
	simData<-data.frame(System=c(1,1,1,1,2,2,2,3),Time=c(18.09,52.07,95.71,145.75,15.02,45.1,82,20.1),Type=c(-1,-1,-1,0,-1,-1,-1,-1),row.names=1:8)
    mle <- mle.vam(System & Time & Type ~ (GQR(0.7) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.03,2.4,0.7)
	lnL<- -293.755998370921
	dlnL<- c(-9842.97490606146,-881.745942144996,-1158.18274977465)
	d2lnL<- matrix(c(-7777.77777777778,-30181.3946412058,-39291.8059448694,-30181.3946412058,-2718.32845215212,-4014.06792884068,-39291.8059448694,-4014.06792884068,-6162.26303978957),nrow=3,byrow=TRUE)
	C<- -24.8250704635939
	dC<-c(2.72891618219407,-6.72454482275228)
	d2C<-matrix(c(-1.33241877266666,-2.82017460407049,-2.82017460407049,-64.9658790434754),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + MC ARAInf + MP QR + MultiSystems + Censorship",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.2) | Weibull(0.001,2.5)) & (QR(0.7)),data=simData)
    theta<-c(0.3,2.2,0.3,0.4)
	lnL<- -32.5275523694307
	dlnL<- c(-28.2832798133722,-18.6640539449763,9.63883415728792,34.2616205525264)
	d2lnL<- matrix(c(-100,-62.2696148233821,44.804062670517,-14.1279314915788,-62.2696148233821,-23.5335804516524,25.1942219659137,16.6907921309135,44.804062670517,25.1942219659137,-18.6975158606326,3.8646832453501,-14.1279314915788,16.6907921309135,3.8646832453501,-111.910753346101),nrow=4,byrow=TRUE)
	C<- -30.0196292995781
	dC<-c(-9.59873223233373,3.1161787129969,36.318390045458)
	d2C<-matrix(c(-2.74251727027051,4.03873049585152,19.4143015947087,4.03873049585152,-5.20711892675317,0.312190876917756,19.4143015947087,0.312190876917756,-103.782193895224),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + MC ARAInf + MP GQR + MultiSystems + Censorship",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.2) | Weibull(0.001,2.5)) & (GQR(0.7)),data=simData)
    theta<-c(0.3,2.2,0.3,0.4)
	lnL<- -32.5275523694307
	dlnL<- c(-28.2832798133722,-18.6640539449763,9.63883415728792,34.2616205525264)
	d2lnL<- matrix(c(-100,-62.2696148233821,44.804062670517,-14.1279314915788,-62.2696148233821,-23.5335804516524,25.1942219659137,16.6907921309135,44.804062670517,25.1942219659137,-18.6975158606326,3.8646832453501,-14.1279314915788,16.6907921309135,3.8646832453501,-111.910753346101),nrow=4,byrow=TRUE)
	C<- -30.0196292995781
	dC<-c(-9.59873223233373,3.1161787129969,36.318390045458)
	d2C<-matrix(c(-2.74251727027051,4.03873049585152,19.4143015947087,4.03873049585152,-5.20711892675317,0.312190876917756,19.4143015947087,0.312190876917756,-103.782193895224),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + MC QR + MP ARAinf + MultiSystems + Censorship",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (QR(0.2) | Weibull(0.001,2.5)) & (ARAInf(0.7)),data=simData)
    theta<-c(0.3,2.2,0.7,0.4)
	lnL<- -22.6714198028668
	dlnL<- c(-24.5487480023091,-10.4391063192085,10.4179303806011,8.17557943541913)
	d2lnL<- matrix(c(-100,-48.2377408144253,-38.6068987313296,38.1382567498831,-48.2377408144253,-16.749614379971,2.43752197636647,18.2846101648813,-38.6068987313296,2.43752197636647,-64.2647694835314,8.69514292754855,38.1382567498831,18.2846101648813,8.69514292754855,-16.7016066567203),nrow=4,byrow=TRUE)
	C<- -20.6878718129505
	dC<-c(-3.92653075979929,15.6302463153278,3.02653483684862)
	d2C<-matrix(c(-3.01060933513669,11.473694438973,3.26669580677216,11.473694438973,-44.9791941884165,0.328568244428881,3.26669580677216,0.328568244428881,-6.28621100077208),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + MC GQR + MP ARAinf + MultiSystems + Censorship",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (GQR(0.2) | Weibull(0.001,2.5)) & (ARAInf(0.7)),data=simData)
    theta<-c(0.3,2.2,0.7,0.4)
	lnL<- -22.6714198028668
	dlnL<- c(-24.5487480023091,-10.4391063192085,10.4179303806011,8.17557943541913)
	d2lnL<- matrix(c(-100,-48.2377408144253,-38.6068987313296,38.1382567498831,-48.2377408144253,-16.749614379971,2.43752197636647,18.2846101648813,-38.6068987313296,2.43752197636647,-64.2647694835314,8.69514292754855,38.1382567498831,18.2846101648813,8.69514292754855,-16.7016066567203),nrow=4,byrow=TRUE)
	C<- -20.6878718129505
	dC<-c(-3.92653075979929,15.6302463153278,3.02653483684862)
	d2C<-matrix(c(-3.01060933513669,11.473694438973,3.26669580677216,11.473694438973,-44.9791941884165,0.328568244428881,3.26669580677216,0.328568244428881,-6.28621100077208),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + MC QR + MP ARA1 + MultiSystems + Censorship",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (QR(0.2) | Weibull(0.001,2.5)) & (ARA1(0.7)),data=simData)
    theta<-c(0.3,2.2,0.7,0.4)
	lnL<- -23.7894585920001
	dlnL<- c(-28.2755439660868,-12.8061999010132,9.95331693250107,7.80657442228432)
	d2lnL<- matrix(c(-100,-56.1280527537744,-40.1556102249964,36.9082400394338,-56.1280527537744,-20.8063184146786,1.42290080543924,19.5922436011137,-40.1556102249964,1.42290080543924,-65.0612496802743,8.45275962596154,36.9082400394338,19.5922436011137,8.45275962596154,-6.30749913607539),nrow=4,byrow=TRUE)
	C<- -21.2826616413254
	dC<-c(-4.63612829622253,15.7984192436662,2.43416349522436)
	d2C<-matrix(c(-3.26433266163505,11.5575946586495,3.275469834776,11.5575946586495,-44.4692405703843,0.4237344848212,3.275469834776,0.4237344848212,-1.25454107692245),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + MC GQR + MP ARA1 + MultiSystems + Censorship",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (GQR(0.2) | Weibull(0.001,2.5)) & (ARA1(0.7)),data=simData)
    theta<-c(0.3,2.2,0.7,0.4)
	lnL<- -23.7894585920001
	dlnL<- c(-28.2755439660868,-12.8061999010132,9.95331693250107,7.80657442228432)
	d2lnL<- matrix(c(-100,-56.1280527537744,-40.1556102249964,36.9082400394338,-56.1280527537744,-20.8063184146786,1.42290080543924,19.5922436011137,-40.1556102249964,1.42290080543924,-65.0612496802743,8.45275962596154,36.9082400394338,19.5922436011137,8.45275962596154,-6.30749913607539),nrow=4,byrow=TRUE)
	C<- -21.2826616413254
	dC<-c(-4.63612829622253,15.7984192436662,2.43416349522436)
	d2C<-matrix(c(-3.26433266163505,11.5575946586495,3.275469834776,11.5575946586495,-44.4692405703843,0.4237344848212,3.275469834776,0.4237344848212,-1.25454107692245),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + MC ARA1 + MP QR + MultiSystems + Censorship",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ARA1(0.2) | Weibull(0.001,2.5)) & (QR(0.7)),data=simData)
    theta<-c(0.3,2.2,0.7,0.4)
	lnL<- -30.7955085887019
	dlnL<- c(-17.7519855519537,-14.0471499881583,5.57391234995491,35.2921127397085)
	d2lnL<- matrix(c(-100,-42.9150891403866,34.496627961891,-10.692957534305,-42.9150891403866,-14.992375117352,13.1220381095624,17.6992813536295,34.496627961891,13.1220381095624,-9.51515958904598,3.10193258543223,-10.692957534305,17.6992813536295,3.10193258543223,-108.017802166014),nrow=4,byrow=TRUE)
	C<- -29.6533223228752
	dC<-c(-9.26099481424148,1.7266358392996,36.4846572514971)
	d2C<-matrix(c(-2.84109562339117,0.921435285537838,19.436402292816,0.921435285537838,-3.24605586075682,0.492868480339134,19.436402292816,0.492868480339134,-103.191786539807),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + MC ARA1 + MP GQR + MultiSystems + Censorship",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ARA1(0.2) | Weibull(0.001,2.5)) & (GQR(0.7)),data=simData)
    theta<-c(0.3,2.2,0.7,0.4)
	lnL<- -30.7955085887019
	dlnL<- c(-17.7519855519537,-14.0471499881583,5.57391234995491,35.2921127397085)
	d2lnL<- matrix(c(-100,-42.9150891403866,34.496627961891,-10.692957534305,-42.9150891403866,-14.992375117352,13.1220381095624,17.6992813536295,34.496627961891,13.1220381095624,-9.51515958904598,3.10193258543223,-10.692957534305,17.6992813536295,3.10193258543223,-108.017802166014),nrow=4,byrow=TRUE)
	C<- -29.6533223228752
	dC<-c(-9.26099481424148,1.7266358392996,36.4846572514971)
	d2C<-matrix(c(-2.84109562339117,0.921435285537838,19.436402292816,0.921435285537838,-3.24605586075682,0.492868480339134,19.436402292816,0.492868480339134,-103.191786539807),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARAInf + PM AGAN + PM QR + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,14)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,2,2,1, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,1,2,2,1,0),row.names=1:22)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (AGAN()+QR(0.7)),data=simData)
    theta<-c(0.3,1.8,0.3,0.7)
	lnL<- -15.3776604100342
	dlnL<- c(-16.8839786764703,-6.49659205082436,1.65928758620766,-1.31833774282327)
	d2lnL<- matrix(c(-55.5555555555556,-23.1026445701935,6.6769924030313,-12.9658877141728,-23.1026445701935,-8.76119478441277,3.54105200760483,-0.249063896896798,6.6769924030313,3.54105200760483,-0.348675974659114,0.424155922414777,-12.9658877141728,-0.249063896896798,0.424155922414777,-9.43522280386179),nrow=4,byrow=TRUE)
	C<- -13.8106937157808
	dC<-c(-3.00874952823924,0.65125156005001,0.639142682205718)
	d2C<-matrix(c(-2.75804429545594,0.857595219934551,1.92573996820204,0.857595219934551,-0.0495361222215684,-0.173844430400474,1.92573996820204,-0.173844430400474,-5.78894052115768),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARAInf + PM ABAO + PM QR + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,14)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,2,2,1, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,1,2,2,1,0),row.names=1:22)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (ABAO()+QR(0.7)),data=simData)
    theta<-c(0.3,1.8,0.3,0.7)
	lnL<- -15.1273150096585
	dlnL<- c(-26.0520900755714,-10.5065811467075,4.99930295489381,-1.03970110052812)
	d2lnL<- matrix(c(-55.5555555555556,-50.458665573324,22.9559022056105,-20.6085274779509,-50.458665573324,-21.9851818182436,16.6728931405193,0.502472049276604,22.9559022056105,16.6728931405193,-6.3852128088804,0.653945073228676,-20.6085274779509,0.502472049276604,0.653945073228676,-32.2667327917739),nrow=4,byrow=TRUE)
	C<- -12.0178248619506
	dC<-c(-1.27489604495264,0.799396794455955,2.73074018890396)
	d2C<-matrix(c(-2.54265291780059,1.89238551886214,4.78762247611762,1.89238551886214,-1.75005615781763,-1.04107047445181,4.78762247611762,-1.04107047445181,-15.9057019425099),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARAInf + PM AGAP + PM QR + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,14)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,2,2,1, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,1,2,2,1,0),row.names=1:22)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (AGAP()+QR(0.7)),data=simData)
    theta<-c(0.3,1.8,0.3,0.7)
	lnL<- -14.5656646781159
	dlnL<- c(-13.0903809462426,-6.11418462574372,1.6634538982099,0.0952034871029612)
	d2lnL<- matrix(c(-55.5555555555556,-21.9554171398778,8.27586368148283,-16.8255121858473,-21.9554171398778,-8.76274038345714,3.93173424400979,2.21273613208318,8.27586368148283,3.93173424400979,-1.07031782202512,0.705936713280325,-16.8255121858473,2.21273613208318,0.705936713280325,-26.9820311867999),nrow=4,byrow=TRUE)
	C<- -13.5368268028835
	dC<-c(-3.21667162714252,0.571266833897918,2.31570970458978)
	d2C<-matrix(c(-2.86489719501896,0.725611074306498,4.58214897692001,0.725611074306498,-0.277272865854988,-0.39088134918004,4.58214897692001,-0.39088134918004,-16.745832783835),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARAInf + PM QAGAN + PM QR + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,14)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,2,2,1, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,1,2,2,1,0),row.names=1:22)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (QAGAN()+QR(0.7)),data=simData)
    theta<-c(0.3,1.8,0.3,0.7)
    logLik(mle,theta,TRUE,FALSE,FALSE)
    logLik(mle,theta,FALSE,TRUE,FALSE)
    logLik(mle,theta,FALSE,FALSE,TRUE)
    contrast(mle,theta,TRUE,FALSE,FALSE)
    contrast(mle,theta,FALSE,TRUE,FALSE)
    contrast(mle,theta,FALSE,FALSE,TRUE)
	lnL<- -14.6375440080987
	dlnL<- c(-12.2768743397197,-6.18133058378875,1.58113962320165,0.29122303249598)
	d2lnL<- matrix(c(-55.5555555555556,-20.8628565336123,6.41649919301129,-16.1721137012039,-20.8628565336123,-8.59055133935837,3.40189211724887,2.38382967201624,6.41649919301129,3.40189211724887,-0.358098930352653,0.446408035084388,-16.1721137012039,2.38382967201624,0.446408035084388,-26.2539585896259),nrow=4,byrow=TRUE)
	C<- -13.7141634750771
	dC<-c(-3.5265344621852,0.764640830505925,2.34912266205614)
	d2C<-matrix(c(-3.00345576787943,0.977648388311007,4.59835097571389,0.977648388311007,-0.0231461647505786,-0.362286740283802,4.59835097571389,-0.362286740283802,-16.6732477149685),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + GQR-log",{
	simData<-data.frame(Time=c(18.09,52.07,95.71,145.75),Type=c(-1,-1,-1,-1),row.names=1:4)
    mle <- mle.vam(Time & Type ~ (GQR(0.7|log) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.03,2.4,0.7)
	lnL<- -313.857271999251
	dlnL<- c(-10545.5445904744,-1049.71538021225,-1077.25045015678)
	d2lnL<- matrix(c(-4444.44444444444,-35478.6147656454,-36271.5544905992,-35478.6147656454,-3544.07657292827,-4129.45460979288,-36271.5544905992,-4129.45460979288,-2814.88047641441),nrow=3,byrow=TRUE)
	C<- -15.0236174259842
	dC<-c(1.353795843943,-2.69009472747077)
	d2C<-matrix(c(-0.784919224791155,-1.93784363307562,-1.93784363307562,-4.37065037672597),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + GQR-sqrt",{
	simData<-data.frame(Time=c(18.09,52.07,95.71,145.75),Type=c(-1,-1,-1,-1),row.names=1:4)
    mle <- mle.vam(Time & Type ~ (GQR(0.7|sqrt) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.03,2.4,0.7)
	lnL<- -244.558483645129
	dlnL<- c(-8207.95806347584,-788.351238674561,-1050.82462756074)
	d2lnL<- matrix(c(-4444.44444444444,-26754.9654997246,-35501.3463704942,-26754.9654997246,-2579.29345017032,-3912.38114820971,-35501.3463704942,-3912.38114820971,-3897.23540002875),nrow=3,byrow=TRUE)
	C<- -14.8642260155259
	dC<-c(1.46759523701393,-2.80862496401469)
	d2C<-matrix(c(-0.759600585249061,-2.10353871070821,-2.10353871070821,-9.8224362844605),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + MC ARA1 + MP GQR-log + MultiSystems + Censorship",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ARA1(0.2) | Weibull(0.001,2.5)) & (GQR(0.7|log)),data=simData)
    theta<-c(0.3,2.2,0.7,0.4)
	lnL<- -25.9295776644323
	dlnL<- c(-19.8315635006478,-11.5441301568133,6.18485740297678,18.5557013220856)
	d2lnL<- matrix(c(-100,-42.8895750190514,36.5331114719639,-16.5532077750091,-42.8895750190514,-15.0267599068787,13.6707008274522,9.78922374712259,36.5331114719639,13.6707008274522,-9.70090792671755,4.84792055163571,-16.5532077750091,9.78922374712259,4.84792055163571,-68.3647543081143),nrow=4,byrow=TRUE)
	C<- -24.5471694757317
	dC<-c(-6.42347607281054,1.82311153869292,20.5320149357969)
	d2C<-matrix(c(-3.11946173172753,0.967523538961955,12.7215337059711,0.967523538961955,-3.10657798789795,0.726777249182269,12.7215337059711,0.726777249182269,-63.5667947445182),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + GQR_ARA1-log",{
	simData<-data.frame(Time=c(18.09,52.07,95.71,145.75),Type=c(-1,-1,-1,-1),row.names=1:4)
    mle <- mle.vam(Time & Type ~ (GQR_ARA1(0.7,-1.2|log) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.03,2.4,0.7,-1.2)
	lnL<- -5227.28210215578
	dlnL<- c(-174523.718461965,-27932.4679482722,-15588.9239900853,2925.83535839578)
	d2lnL<- matrix(c(-4444.44444444444,-931711.328605459,-519891.725970003,97574.7515774708,-931711.328605459,-149124.678752068,-89338.3200091043,17794.0729848851,-519891.725970003,-89338.3200091043,-29458.0754936798,8538.61396441627,97574.7515774708,17794.0729848851,8538.61396441627,-474.24229613857),nrow=4,byrow=TRUE)
	C<- -20.2814577676388
	dC<-c(-2.46616536000733,-4.07878379460467,0.827470378552931)
	d2C<-matrix(c(-0.707449731980772,-2.33853540176398,0.658809838099249,-2.33853540176398,3.22719216928062,0.213056341848732,0.658809838099249,0.213056341848732,0.406344931337982),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + GQR_ARAInf",{
	simData<-data.frame(Time=c(18.09,52.07,95.71,145.75),Type=c(-1,-1,-1,-1),row.names=1:4)
    mle <- mle.vam(Time & Type ~ (GQR_ARAInf(0.1,0.5) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.03,2.4,0.7,-1.2)
	lnL<- -6920.95828002956
	dlnL<- c(-230985.44201529,-41289.7322048457,-33910.4905043345,8972.05361316284)
	d2lnL<- matrix(c(-4444.44444444444,-1376981.61444881,-1130746.37557052,299160.41537896,-1376981.61444881,-247372.240704713,-212484.726352308,61190.0304852329,-1130746.37557052,-212484.726352308,-136997.192671136,43010.4052016972,299160.41537896,61190.0304852329,43010.4052016972,-9344.88981359903),nrow=4,byrow=TRUE)
	C<- -21.2263586825304
	dC<-c(-4.1153551612053,-7.66919857501101,2.41875651139994)
	d2C<-matrix(c(-1.41796907601273,-3.61075961185818,2.48367083926731,-3.61075961185818,3.94864357549099,0.8268593639824,2.48367083926731,0.8268593639824,0.485567016552806),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.0000000000001))
}
)

test_that("Weibull + MC ARA1 + MP GQ_RARA1-log + MultiSystems + Censorship",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ARA1(0.2) | Weibull(0.001,2.5)) & (GQR_ARA1(0.7,-1.3|log)),data=simData)
    theta<-c(0.3,2.2,0.7,0.4,-0.9)
	lnL<- -33.2961906877108
	dlnL<- c(-92.8760104131109,-46.6801033336665,9.34931942198527,-47.9119117936478,10.470539758891)
	d2lnL<- matrix(c(-100,-200.417044082063,42.793609075184,-201.850249846753,44.1676540638691,-200.417044082063,-108.270211854658,21.9813689091466,-131.618058974325,35.9779791929647,42.793609075184,21.9813689091466,-8.36124676707121,10.1156562971593,-0.419667547402598,-201.850249846753,-131.618058974325,10.1156562971593,-61.6035552690088,33.7481794183083,44.1676540638691,35.9779791929647,-0.419667547402598,33.7481794183083,-2.57572452700716),nrow=5,byrow=TRUE)
	C<- -18.1231933425977
	dC<-c(-1.23444943291652,-0.35436399527184,-2.14127069707828,0.455284173279421)
	d2C<-matrix(c(-3.89657698655457,-1.94313484129869,-6.7909902653984,1.75655268681334,-1.94313484129869,-4.09281601550727,-3.06754108013221,0.909378943225631,-6.7909902653984,-3.06754108013221,-12.3460043040192,3.37568923785665,1.75655268681334,0.909378943225631,3.37568923785665,-0.446933256627635),nrow=4,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + MC AGAN + MP GQ_RARA1-log + MultiSystems + Censorship",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,8)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.52),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,0),row.names=1:24)
    mle <- mle.vam(System & Time & Type ~ (AGAN() | Weibull(0.001,2.5)) & (GQR_ARA1(0.7,-1.3|log)),data=simData)
    theta<-c(0.3,2.2,0.4,-0.9)
	lnL<- -26.0110174352452
	dlnL<- c(-55.3770729201669,-26.0264586110469,-24.1768593541568,5.75212431635257)
	d2lnL<- matrix(c(-100,-114.966549141451,-100.876259060329,24.5911031510572,-114.966549141451,-59.7228171559624,-63.1479745895377,19.3873280431689,-100.876259060329,-63.1479745895377,-23.7034729937036,17.679086547779,24.5911031510572,19.3873280431689,17.679086547779,-1.47232544091911),nrow=4,byrow=TRUE)
	C<- -18.8108175042004
	dC<-c(-3.65566013752983,-4.5478239603539,0.967057556134173)
	d2C<-matrix(c(-5.87229883911108,-7.39016980734201,2.44324201647881,-7.39016980734201,-4.81734728102459,3.39712910973041,2.44324201647881,3.39712910973041,-0.250226518299337),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + MC ARA1 + MP GQ_RARAInf-log + MultiSystems + Censorship",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,7)),Time=c(3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0),row.names=1:23)
    mle <- mle.vam(System & Time & Type ~ (ARA1(0.2) | Weibull(0.001,2.5)) & (GQR_ARAInf(0.7,-1.3|log)),data=simData)
    theta<-c(0.3,2.2,0.7,0.4,-0.9)
	lnL<- -45.1767974450898
	dlnL<- c(-136.836651977583,-92.1950025296914,11.4275781876523,-88.5648534727329,34.3648437225561)
	d2lnL<- matrix(c(-100,-355.76555726871,49.6772787472738,-336.030699309883,127.652490859454,-355.76555726871,-262.782172493883,29.4935147997201,-279.640452671131,136.238226683968,49.6772787472738,29.4935147997201,-8.50644897780949,16.5536939999478,-3.34889216065876,-336.030699309883,-279.640452671131,16.5536939999478,-88.2082103062016,110.873727273422,127.652490859454,136.238226683968,-3.34889216065876,110.873727273422,-30.0687703335966),nrow=5,byrow=TRUE)
	C<- -19.568162233621
	dC<-c(-4.65710153153946,-0.795765584639633,-5.88281332709037,2.95530764264158)
	d2C<-matrix(c(-7.85285540201067,-2.78664457162768,-10.5680421360252,7.12689306214278,-2.78664457162768,-4.14065180208066,-2.82691291028266,1.32051846411931,-10.5680421360252,-2.82691291028266,-2.44964958457988,6.71488713174694,7.12689306214278,1.32051846411931,6.71488713174694,-1.62344453972361),nrow=4,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.0000000000001))
}
)

test_that("Weibull + ARAm-m=2",{
	simData<-data.frame(Time=c(18.09,52.07,95.71,145.75,198.7,220.9),Type=c(-1,-1,-1,-1,-1,-1),row.names=1:6)
    mle <- mle.vam(Time & Type ~ (ARAm(0.5|2) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.03,2.4,0.7)
	lnL<- -2600.46156675627
	dlnL<- c(-87032.7099834185,-10989.6661601188,6745.47891744838)
	d2lnL<- matrix(c(-6666.66666666667,-367174.735520771,225153.425274492,-367174.735520771,-46307.9946188132,32956.7724964216,225153.425274492,32956.7724964216,-16880.4479240297),nrow=3,byrow=TRUE)
	C<- -25.948370978711
	dC<-c(0.32106229767103,6.36255820206048)
	d2C<-matrix(c(-0.909264674311954,3.87406673594902,3.87406673594902,-0.394244003389716),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + ARAm-m=3",{
	simData<-data.frame(Time=c(18.09,52.07,95.71,145.75,198.7,220.9),Type=c(-1,-1,-1,-1,-1,-1),row.names=1:6)
    mle <- mle.vam(Time & Type ~ (ARAm(0.5|3) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.03,2.4,0.7)
	lnL<- -2436.77330082115
	dlnL<- c(-81567.7532962045,-10188.2094022023,5719.57195795803)
	d2lnL<- matrix(c(-6666.66666666667,-340453.309431142,190929.671896146,-340453.309431142,-42483.9898494684,27534.7616151142,190929.671896146,27534.7616151142,-14554.2780650748),nrow=3,byrow=TRUE)
	C<- -25.8206271210325
	dC<-c(0.407907102520102,5.69194547894592)
	d2C<-matrix(c(-0.936417820176478,3.08815779480951,3.08815779480951,-0.920858811778887),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + ARAm-m=1",{
	simData<-data.frame(Time=c(18.09,52.07,95.71,145.75,198.7,220.9),Type=c(-1,-1,-1,-1,-1,-1),row.names=1:6)
    mle <- mle.vam(Time & Type ~ (ARAm(0.5|1) | Weibull(0.001,2.5)),data=simData)
    mle2 <- mle.vam(Time & Type ~ (ARA1(0.5) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.03,2.4,0.7)
	lnL<- logLik(mle2,theta,TRUE,FALSE,FALSE)
	dlnL<- logLik(mle2,theta,FALSE,TRUE,FALSE)
	d2lnL<- logLik(mle2,theta,FALSE,FALSE,TRUE)
	C<-contrast(mle2,theta,TRUE,FALSE,FALSE)
	dC<-contrast(mle2,theta,FALSE,TRUE,FALSE)
	d2C<-contrast(mle2,theta,FALSE,FALSE,TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + ARAm-m=Inf",{
	simData<-data.frame(Time=c(18.09,52.07,95.71,145.75,198.7,220.9),Type=c(-1,-1,-1,-1,-1,-1),row.names=1:6)
    mle <- mle.vam(Time & Type ~ (ARAm(0.5|8) | Weibull(0.001,2.5)),data=simData)
    mle2 <- mle.vam(Time & Type ~ (ARAInf(0.5) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.03,2.4,0.7)
	lnL<- logLik(mle2,theta,TRUE,FALSE,FALSE)
	dlnL<- logLik(mle2,theta,FALSE,TRUE,FALSE)
	d2lnL<- logLik(mle2,theta,FALSE,FALSE,TRUE)
	C<-contrast(mle2,theta,TRUE,FALSE,FALSE)
	dC<-contrast(mle2,theta,FALSE,TRUE,FALSE)
	d2C<-contrast(mle2,theta,FALSE,FALSE,TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARAm-m=2 + PM ARAm-m=4 + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,14)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0),row.names=1:22)
    theta<-c(0.3,1.8,0.3,0.7)
    mle <- mle.vam(System & Time & Type ~ (ARAm(0.3|2) | Weibull(0.001,2.5)) & (ARAm(0.6|4)),data=simData)
	lnL<- -20.1504767593027
	dlnL<- c(-29.4476398159832,-10.1914497132837,4.64183577052888,7.06399852898271)
	d2lnL<- matrix(c(-88.8888888888889,-55.5084151127721,21.8799936040503,65.7589094670528,-55.5084151127721,-17.226152459112,12.4350883648324,26.1191546857653,21.8799936040503,12.4350883648324,-0.565415653412342,-11.0052964447512,65.7589094670528,26.1191546857653,-11.0052964447512,-34.6307206131453),nrow=4,byrow=TRUE)
	C<- -17.2679964267276
	dC<-c(-1.45254558638836,1.19718389181072,-3.28868050453445)
	d2C<-matrix(c(-1.65379858663581,1.56285390950733,-5.16845746266277,1.56285390950733,0.537921529475571,0.93175674258043,-5.16845746266277,0.93175674258043,-8.3953177399354),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARAm-m=4 + PM ARAm-m=2 + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,14)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0),row.names=1:22)
    theta<-c(0.3,1.8,0.3,0.7)
    mle <- mle.vam(System & Time & Type ~ (ARAm(0.3|4) | Weibull(0.001,2.5)) & (ARAm(0.6|2)),data=simData)
	lnL<- -21.2311122143017
	dlnL<- c(-38.0598977198794,-14.5704485906476,9.67269834242892,9.68850554138202)
	d2lnL<- matrix(c(-88.8888888888889,-76.3677526880262,47.5797934336762,68.2956526339584,-76.3677526880262,-27.8962044371839,31.1680564206178,36.300015570557,47.5797934336762,31.1680564206178,-10.4594408536748,-23.6602841209597,68.2956526339584,36.300015570557,-23.6602841209597,-17.2696560246775),nrow=4,byrow=TRUE)
	C<- -16.9072016909912
	dC<-c(-1.098937083971,1.27947329434443,-2.35906215480096)
	d2C<-matrix(c(-1.80841986731932,2.5205717042266,-2.94233502985972,2.5205717042266,-0.663494335954779,-2.0138163778875,-2.94233502985972,-2.0138163778875,-7.30859939934244),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM QR + PM ARAm-m=3 + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,14)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0),row.names=1:22)
    theta<-c(0.3,1.8,1.3,0.7)
    mle <- mle.vam(System & Time & Type ~ (QR(0.3) | Weibull(0.001,2.5)) & (ARAm(0.6|2)),data=simData)
	lnL<- -39.8773909185144
	dlnL<-c(-113.496292381551,-54.536361547164,-205.338866836277,9.25106291122166)
	d2lnL<-matrix(c(-88.8888888888889,-209.792361406732,-753.693658685026,56.2966050971949,-209.792361406732,-107.518863114582,-513.017025854981,28.1201421190358,-753.693658685026,-513.017025854981,-1439.71841158142,63.1526616736015,56.2966050971949,28.1201421190358,63.1526616736015,-9.17483031281349),nrow=4,byrow=TRUE)
	C<- -19.103634347171
	dC<- c(-3.57285017870748,-22.2489054803471,-4.42470990833311)
	d2C<- matrix(c(-4.53266947644281,-23.8722488024133,-7.19041963030991,-23.8722488024133,-55.5299981102586,-5.26318833097213,-7.19041963030991,-5.26318833097213,-5.52825758512622),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARAInf + PM ARAm-m=3 + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,14)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0),row.names=1:22)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.3) | Weibull(0.001,2.5)) & (ARAm(0.6|3)),data=simData)
    mle2 <- mle.vam(System & Time & Type ~ (ARAm(0.3|22) | Weibull(0.001,2.5)) & (ARAm(0.6|3)),data=simData)
    theta<-c(0.03,2.4,0.3,0.7)
	lnL<- logLik(mle2,theta,TRUE,FALSE,FALSE)
	dlnL<- logLik(mle2,theta,FALSE,TRUE,FALSE)
	d2lnL<- logLik(mle2,theta,FALSE,FALSE,TRUE)
	C<-contrast(mle2,theta,TRUE,FALSE,FALSE)
	dC<-contrast(mle2,theta,FALSE,TRUE,FALSE)
	d2C<-contrast(mle2,theta,FALSE,FALSE,TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARA1 + PM ARAm-m=3 + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,14)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0),row.names=1:22)
    mle <- mle.vam(System & Time & Type ~ (ARA1(0.3) | Weibull(0.001,2.5)) & (ARAm(0.6|3)),data=simData)
    mle2 <- mle.vam(System & Time & Type ~ (ARAm(0.3|1) | Weibull(0.001,2.5)) & (ARAm(0.6|3)),data=simData)
    theta<-c(0.03,2.4,0.3,0.7)
	lnL<- logLik(mle2,theta,TRUE,FALSE,FALSE)
	dlnL<- logLik(mle2,theta,FALSE,TRUE,FALSE)
	d2lnL<- logLik(mle2,theta,FALSE,FALSE,TRUE)
	C<-contrast(mle2,theta,TRUE,FALSE,FALSE)
	dC<-contrast(mle2,theta,FALSE,TRUE,FALSE)
	d2C<-contrast(mle2,theta,FALSE,FALSE,TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARAm-m3 + PM ARAInf + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,14)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0),row.names=1:22)
    mle <- mle.vam(System & Time & Type ~ (ARAm(0.3|3) | Weibull(0.001,2.5)) & (ARAInf(0.6)),data=simData)
    mle2 <- mle.vam(System & Time & Type ~ (ARAm(0.3|3) | Weibull(0.001,2.5)) & (ARAm(0.6|22)),data=simData)
    theta<-c(0.03,2.4,0.3,0.7)
	lnL<- logLik(mle2,theta,TRUE,FALSE,FALSE)
	dlnL<- logLik(mle2,theta,FALSE,TRUE,FALSE)
	d2lnL<- logLik(mle2,theta,FALSE,FALSE,TRUE)
	C<-contrast(mle2,theta,TRUE,FALSE,FALSE)
	dC<-contrast(mle2,theta,FALSE,TRUE,FALSE)
	d2C<-contrast(mle2,theta,FALSE,FALSE,TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARAm-m=3 + PM ARA1 + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,14)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0),row.names=1:22)
    mle <- mle.vam(System & Time & Type ~ (ARAm(0.3|3) | Weibull(0.001,2.5)) & (ARA1(0.6)),data=simData)
    mle2 <- mle.vam(System & Time & Type ~ (ARAm(0.3|3) | Weibull(0.001,2.5)) & (ARAm(0.6|1)),data=simData)
    theta<-c(0.03,2.4,0.3,0.7)
	lnL<- logLik(mle2,theta,TRUE,FALSE,FALSE)
	dlnL<- logLik(mle2,theta,FALSE,TRUE,FALSE)
	d2lnL<- logLik(mle2,theta,FALSE,FALSE,TRUE)
	C<-contrast(mle2,theta,TRUE,FALSE,FALSE)
	dC<-contrast(mle2,theta,FALSE,TRUE,FALSE)
	d2C<-contrast(mle2,theta,FALSE,FALSE,TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARAm-m=3 + PM QR + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,14)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0),row.names=1:22)
    theta<-c(0.3,1.8,0.3,0.7)
    mle <- mle.vam(System & Time & Type ~ (ARAm(0.3|2) | Weibull(0.001,2.5)) & (QR(0.6)),data=simData)
	lnL<- -30.1766718759053
	dlnL<- c(4.04792345878618,-12.6501273642652,0.613707394378638,58.2110225533236)
	d2lnL<- matrix(c(-88.8888888888889,-17.2123801970225,3.19172509693457,-20.2489724413023,-17.2123801970225,-8.59612513999474,1.13460477244041,34.8653501316305,3.19172509693457,1.13460477244041,-0.14487695793176,2.37457191752127,-20.2489724413023,34.8653501316305,2.37457191752127,-120.718396691945),nrow=4,byrow=TRUE)
	C<- -30.0739654098369
	dC<-c(-13.5742423375074,0.78506782906508,57.1238760685502)
	d2C<-matrix(c(-5.05992318248314,0.555517025337616,40.1634277012813,0.555517025337616,0.0149341030595374,1.78892782197386,40.1634277012813,1.78892782197386,-119.4756615276),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + GQR_ARAm-m=3",{
	simData<-data.frame(Time=c(18.09,52.07,95.71,145.75,198.7,220.9,230),Type=c(-1,-1,-1,-1,-1,-1,0),row.names=1:7)
    mle <- mle.vam(Time & Type ~ (GQR_ARAm(1.2,0.5|3) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.03,2.4,1.3,0.7)
	lnL<- -19524.0996600096
	dlnL<- c(-651440.845101926,-98931.83785181,-130451.573242292,44948.9374958108)
	d2lnL<- matrix(c(-6666.66666666667,-3298691.66087083,-4349260.39370564,1498505.11901772,-3298691.66087083,-501269.384415686,-717595.854835935,254816.142137408,-4349260.39370564,-717595.854835935,-837390.260758198,292604.491724178,1498505.11901772,254816.142137408,292604.491724178,-115028.435473886),nrow=4,byrow=TRUE)
	C<- -29.4078957809114
	dC<-c(-1.46081479729706,-13.8073591192265,7.58145006849591)
	d2C<-matrix(c(-1.13865627103432,-7.02814050027214,3.92396572758071,-7.02814050027214,-8.65037145925458,2.91612651943195,3.92396572758071,2.91612651943195,-2.45231214454718),nrow=3,byrow=TRUE)
	expect_that(logLik(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM GQR_ARAm-m=3 vs ARAm + PM ARA1 + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,14)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0),row.names=1:22)
    mle <- mle.vam(System & Time & Type ~ (GQR_ARAm(1,0.3|3) | Weibull(0.001,2.5)) & (ARA1(0.6)),data=simData)
    mle2 <- mle.vam(System & Time & Type ~ (ARAm(0.3|3) | Weibull(0.001,2.5)) & (ARAm(0.6|1)),data=simData)
    theta<-c(0.03,2.4,1,0.3,0.7)
    theta2<-c(0.03,2.4,0.3,0.7)
	lnL<- logLik(mle2,theta2,TRUE,FALSE,FALSE)
	dlnL<- logLik(mle2,theta2,FALSE,TRUE,FALSE)
	d2lnL<- logLik(mle2,theta2,FALSE,FALSE,TRUE)
	C<-contrast(mle2,theta2,TRUE,FALSE,FALSE)
	dC<-contrast(mle2,theta2,FALSE,TRUE,FALSE)
	d2C<-contrast(mle2,theta2,FALSE,FALSE,TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE)[c(1,2,4,5)],equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE)[c(1,2,4,5),c(1,2,4,5)],equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE)[c(1,3,4)],equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE)[c(1,3,4),c(1,3,4)],equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM GQR_ARAm-m=22-fun=log vs GQR + PM ARA1 + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,14)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0),row.names=1:22)
    mle <- mle.vam(System & Time & Type ~ (GQR_ARAm(0.7,1|22,log) | Weibull(0.001,2.5)) & (ARA1(0.6)),data=simData)
    mle2 <- mle.vam(System & Time & Type ~ (GQR(0.7|log) | Weibull(0.001,2.5)) & (ARAm(0.6|1)),data=simData)
    theta<-c(0.03,2.4,0.7,1,0.8)
    theta2<-c(0.03,2.4,0.7,0.8)
	lnL<- logLik(mle2,theta2,TRUE,FALSE,FALSE)
	dlnL<- logLik(mle2,theta2,FALSE,TRUE,FALSE)
	d2lnL<- logLik(mle2,theta2,FALSE,FALSE,TRUE)
	C<-contrast(mle2,theta2,TRUE,FALSE,FALSE)
	dC<-contrast(mle2,theta2,FALSE,TRUE,FALSE)
	d2C<-contrast(mle2,theta2,FALSE,FALSE,TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE)[c(1,2,3,5)],equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE)[c(1,2,3,5),c(1,2,3,5)],equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE)[c(1,2,4)],equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE)[c(1,2,4),c(1,2,4)],equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM GQR_ARAm-m=1-fun=log vs GQR_ARA1 + PM ARA1 + mutlisystems",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,14)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0),row.names=1:22)
    mle <- mle.vam(System & Time & Type ~ (GQR_ARAm(0.7,1.3|1,log) | Weibull(0.001,2.5)) & (ARA1(0.6)),data=simData)
    mle2 <- mle.vam(System & Time & Type ~ (GQR_ARA1(0.7,1.3|log) | Weibull(0.001,2.5)) & (ARAm(0.6|1)),data=simData)
    theta<-c(0.03,2.4,0.7,1.3,0.8)
    theta2<-c(0.03,2.4,0.7,1.3,0.8)
	lnL<- logLik(mle2,theta2,TRUE,FALSE,FALSE)
	dlnL<- logLik(mle2,theta2,FALSE,TRUE,FALSE)
	d2lnL<- logLik(mle2,theta2,FALSE,FALSE,TRUE)
	C<-contrast(mle2,theta2,TRUE,FALSE,FALSE)
	dC<-contrast(mle2,theta2,FALSE,TRUE,FALSE)
	d2C<-contrast(mle2,theta2,FALSE,FALSE,TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("DataList + Weibull + MC ARAInf",{
	simData<-data.frame(Time=c(3.36,4.04,4.97,5.16),Type=c(-1,-1,-1,-1),row.names=1:4)
	simData<-list(simData)
	names(simData)<-"System1"
	mle <- mle.vam(Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)),data=simData)
	theta<-c(0.3,0.8,0.6)
	lnL<- -7.37830963135462
	dlnL<- c(9.33348796771948,5.77076155284033,1.19923836457015)
	d2lnL<- matrix(c(-44.4444444444444,-5.26230480903023,-0.723044448779292,-5.26230480903023,-7.7471885008781,-6.30726171420585,-0.723044448779292,-6.30726171420585,0.435684125802398),nrow=3,byrow=TRUE)
	C<- -5.36231016699152
	dC<-c(2.08694474533596,0.693079297460398)
	d2C<-matrix(c(-4.31732300351524,-3.55104259186756,-3.55104259186756,-0.351517905180765),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("OthersNames + Weibull + MC ARAInf",{
	simData<-data.frame(T=c(3.36,4.04,4.97,5.16),C=c(-1,-1,-1,-1),row.names=1:4)
	#simData<-list(simData)
	#names(simData)<-"System1"
	mle <- mle.vam(T & C ~ (ARAInf(0.4) | Weibull(0.001,2.5)),data=simData)
	theta<-c(0.3,0.8,0.6)
	lnL<- -7.37830963135462
	dlnL<- c(9.33348796771948,5.77076155284033,1.19923836457015)
	d2lnL<- matrix(c(-44.4444444444444,-5.26230480903023,-0.723044448779292,-5.26230480903023,-7.7471885008781,-6.30726171420585,-0.723044448779292,-6.30726171420585,0.435684125802398),nrow=3,byrow=TRUE)
	C<- -5.36231016699152
	dC<-c(2.08694474533596,0.693079297460398)
	d2C<-matrix(c(-4.31732300351524,-3.55104259186756,-3.55104259186756,-0.351517905180765),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("OthersNames2 + Weibull + MC ARAInf",{
	simData<-data.frame(C=c(-1,-1,-1,-1),T=c(3.36,4.04,4.97,5.16),row.names=1:4)
	#simData<-list(simData)
	#names(simData)<-"System1"
	mle <- mle.vam(T & C ~ (ARAInf(0.4) | Weibull(0.001,2.5)),data=simData)
	theta<-c(0.3,0.8,0.6)
	lnL<- -7.37830963135462
	dlnL<- c(9.33348796771948,5.77076155284033,1.19923836457015)
	d2lnL<- matrix(c(-44.4444444444444,-5.26230480903023,-0.723044448779292,-5.26230480903023,-7.7471885008781,-6.30726171420585,-0.723044448779292,-6.30726171420585,0.435684125802398),nrow=3,byrow=TRUE)
	C<- -5.36231016699152
	dC<-c(2.08694474533596,0.693079297460398)
	d2C<-matrix(c(-4.31732300351524,-3.55104259186756,-3.55104259186756,-0.351517905180765),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("OthersNames + DataList + Weibull + MC ARAInf",{
	simData<-data.frame(T=c(3.36,4.04,4.97,5.16),C=c(-1,-1,-1,-1),row.names=1:4)
	simData<-list(simData)
	names(simData)<-"System1"
	mle <- mle.vam(T & C ~ (ARAInf(0.4) | Weibull(0.001,2.5)),data=simData)
	theta<-c(0.3,0.8,0.6)
	lnL<- -7.37830963135462
	dlnL<- c(9.33348796771948,5.77076155284033,1.19923836457015)
	d2lnL<- matrix(c(-44.4444444444444,-5.26230480903023,-0.723044448779292,-5.26230480903023,-7.7471885008781,-6.30726171420585,-0.723044448779292,-6.30726171420585,0.435684125802398),nrow=3,byrow=TRUE)
	C<- -5.36231016699152
	dC<-c(2.08694474533596,0.693079297460398)
	d2C<-matrix(c(-4.31732300351524,-3.55104259186756,-3.55104259186756,-0.351517905180765),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("OthersNames + DataList2 + Weibull + MC ARAInf",{
	simData<-data.frame(T=c(3.36,4.04,4.97,5.16),C=c(-1,-1,-1,-1),row.names=1:4)
	simData<-list(simData)
	names(simData)<-"System1"
	mle <- mle.vam(System & T & C ~ (ARAInf(0.4) | Weibull(0.001,2.5)),data=simData)
	theta<-c(0.3,0.8,0.6)
	lnL<- -7.37830963135462
	dlnL<- c(9.33348796771948,5.77076155284033,1.19923836457015)
	d2lnL<- matrix(c(-44.4444444444444,-5.26230480903023,-0.723044448779292,-5.26230480903023,-7.7471885008781,-6.30726171420585,-0.723044448779292,-6.30726171420585,0.435684125802398),nrow=3,byrow=TRUE)
	C<- -5.36231016699152
	dC<-c(2.08694474533596,0.693079297460398)
	d2C<-matrix(c(-4.31732300351524,-3.55104259186756,-3.55104259186756,-0.351517905180765),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("OthersNames2 + DataList2 + Weibull + MC ARAInf",{
	simData<-data.frame(T=c(3.36,4.04,4.97,5.16),C=c(-1,-1,-1,-1),row.names=1:4)
	simData<-list(simData)
	names(simData)<-"S1"
	mle <- mle.vam(S & T & C ~ (ARAInf(0.4) | Weibull(0.001,2.5)),data=simData)
	theta<-c(0.3,0.8,0.6)
	lnL<- -7.37830963135462
	dlnL<- c(9.33348796771948,5.77076155284033,1.19923836457015)
	d2lnL<- matrix(c(-44.4444444444444,-5.26230480903023,-0.723044448779292,-5.26230480903023,-7.7471885008781,-6.30726171420585,-0.723044448779292,-6.30726171420585,0.435684125802398),nrow=3,byrow=TRUE)
	C<- -5.36231016699152
	dC<-c(2.08694474533596,0.693079297460398)
	d2C<-matrix(c(-4.31732300351524,-3.55104259186756,-3.55104259186756,-0.351517905180765),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("DataList + Weibull + CM GQR_ARAm-m=1-fun=log vs GQR_ARA1 + PM ARA1 + mutlisystems",{
	simData1<-data.frame(Time=c(3.36,4.04,4.97,5.16),Type=c(1,1,-1,1),row.names=1:4)
	simData2<-data.frame(Time=c(0.78,2.36,4.05,4.97),Type=c(-1,1,1,0),row.names=1:4)
	simData3<-data.frame(Time=c(2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0),row.names=1:14)
	simData<-list(simData1,simData2,simData3)
	names(simData)<-c("System1","System2","System3")
    mle <- mle.vam(System &Time & Type ~ (GQR_ARAm(0.7,1.3|1,log) | Weibull(0.001,2.5)) & (ARA1(0.6)),data=simData)
    mle2 <- mle.vam(System & Time & Type ~ (GQR_ARA1(0.7,1.3|log) | Weibull(0.001,2.5)) & (ARAm(0.6|1)),data=simData)
    theta<-c(0.03,2.4,0.7,1.3,0.8)
    theta2<-c(0.03,2.4,0.7,1.3,0.8)
	lnL<- logLik(mle2,theta2,TRUE,FALSE,FALSE)
	dlnL<- logLik(mle2,theta2,FALSE,TRUE,FALSE)
	d2lnL<- logLik(mle2,theta2,FALSE,FALSE,TRUE)
	C<-contrast(mle2,theta2,TRUE,FALSE,FALSE)
	dC<-contrast(mle2,theta2,FALSE,TRUE,FALSE)
	d2C<-contrast(mle2,theta2,FALSE,FALSE,TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("DataList2 + Weibull + CM GQR_ARAm-m=1-fun=log vs GQR_ARA1 + PM ARA1 + mutlisystems",{
	simData1<-data.frame(Time=c(3.36,4.04,4.97,5.16),Type=c(1,1,-1,1),row.names=1:4)
	simData2<-data.frame(Time=c(0.78,2.36,4.05,4.97),Type=c(-1,1,1,0),row.names=1:4)
	simData3<-data.frame(Time=c(2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0),row.names=1:14)
	simData<-list(simData1,simData2,simData3)
	names(simData)<-c("System1","System2","System3")
    mle <- mle.vam(Time & Type ~ (GQR_ARAm(0.7,1.3|1,log) | Weibull(0.001,2.5)) & (ARA1(0.6)),data=simData)
    mle2 <- mle.vam( Time & Type ~ (GQR_ARA1(0.7,1.3|log) | Weibull(0.001,2.5)) & (ARAm(0.6|1)),data=simData)
    theta<-c(0.03,2.4,0.7,1.3,0.8)
    theta2<-c(0.03,2.4,0.7,1.3,0.8)
	lnL<- logLik(mle2,theta2,TRUE,FALSE,FALSE)
	dlnL<- logLik(mle2,theta2,FALSE,TRUE,FALSE)
	d2lnL<- logLik(mle2,theta2,FALSE,FALSE,TRUE)
	C<-contrast(mle2,theta2,TRUE,FALSE,FALSE)
	dC<-contrast(mle2,theta2,FALSE,TRUE,FALSE)
	d2C<-contrast(mle2,theta2,FALSE,FALSE,TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("OthersNames + DataList + Weibull + CM GQR_ARAm-m=1-fun=log vs GQR_ARA1 + PM ARA1 + mutlisystems",{
	simData1<-data.frame(T=c(3.36,4.04,4.97,5.16),C=c(1,1,-1,1),row.names=1:4)
	simData2<-data.frame(C=c(-1,1,1,0),T=c(0.78,2.36,4.05,4.97),row.names=1:4)
	simData3<-data.frame(T=c(2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),C=c(1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0),row.names=1:14)
	simData<-list(simData1,simData2,simData3)
	names(simData)<-c("System1","System2","System3")
    mle <- mle.vam(System &T & C ~ (GQR_ARAm(0.7,1.3|1,log) | Weibull(0.001,2.5)) & (ARA1(0.6)),data=simData)
    mle2 <- mle.vam(System & T & C ~ (GQR_ARA1(0.7,1.3|log) | Weibull(0.001,2.5)) & (ARAm(0.6|1)),data=simData)
    theta<-c(0.03,2.4,0.7,1.3,0.8)
    theta2<-c(0.03,2.4,0.7,1.3,0.8)
	lnL<- logLik(mle2,theta2,TRUE,FALSE,FALSE)
	dlnL<- logLik(mle2,theta2,FALSE,TRUE,FALSE)
	d2lnL<- logLik(mle2,theta2,FALSE,FALSE,TRUE)
	C<-contrast(mle2,theta2,TRUE,FALSE,FALSE)
	dC<-contrast(mle2,theta2,FALSE,TRUE,FALSE)
	d2C<-contrast(mle2,theta2,FALSE,FALSE,TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("OthersNames + DataList2 + Weibull + CM GQR_ARAm-m=1-fun=log vs GQR_ARA1 + PM ARA1 + mutlisystems",{
	simData1<-data.frame(T=c(3.36,4.04,4.97,5.16),C=c(1,1,-1,1),row.names=1:4)
	simData2<-data.frame(C=c(-1,1,1,0),T=c(0.78,2.36,4.05,4.97),row.names=1:4)
	simData3<-data.frame(T=c(2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),C=c(1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0),row.names=1:14)
	simData<-list(simData1,simData2,simData3)
	names(simData)<-c("System1","System2","System3")
    mle <- mle.vam(T & C ~ (GQR_ARAm(0.7,1.3|1,log) | Weibull(0.001,2.5)) & (ARA1(0.6)),data=simData)
    mle2 <- mle.vam( T & C ~ (GQR_ARA1(0.7,1.3|log) | Weibull(0.001,2.5)) & (ARAm(0.6|1)),data=simData)
    theta<-c(0.03,2.4,0.7,1.3,0.8)
    theta2<-c(0.03,2.4,0.7,1.3,0.8)
	lnL<- logLik(mle2,theta2,TRUE,FALSE,FALSE)
	dlnL<- logLik(mle2,theta2,FALSE,TRUE,FALSE)
	d2lnL<- logLik(mle2,theta2,FALSE,FALSE,TRUE)
	C<-contrast(mle2,theta2,TRUE,FALSE,FALSE)
	dC<-contrast(mle2,theta2,FALSE,TRUE,FALSE)
	d2C<-contrast(mle2,theta2,FALSE,FALSE,TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("OthersNames2 + DataList + Weibull + CM GQR_ARAm-m=1-fun=log vs GQR_ARA1 + PM ARA1 + mutlisystems",{
	simData1<-data.frame(T=c(3.36,4.04,4.97,5.16),C=c(1,1,-1,1),row.names=1:4)
	simData2<-data.frame(C=c(-1,1,1,0),T=c(0.78,2.36,4.05,4.97),row.names=1:4)
	simData3<-data.frame(T=c(2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),C=c(1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0),row.names=1:14)
	simData<-list(simData1,simData2,simData3)
	names(simData)<-c("S1","S2","S3")
    mle <- mle.vam(S &T & C ~ (GQR_ARAm(0.7,1.3|1,log) | Weibull(0.001,2.5)) & (ARA1(0.6)),data=simData)
    mle2 <- mle.vam( T & C ~ (GQR_ARA1(0.7,1.3|log) | Weibull(0.001,2.5)) & (ARAm(0.6|1)),data=simData)
    theta<-c(0.03,2.4,0.7,1.3,0.8)
    theta2<-c(0.03,2.4,0.7,1.3,0.8)
	lnL<- logLik(mle2,theta2,TRUE,FALSE,FALSE)
	dlnL<- logLik(mle2,theta2,FALSE,TRUE,FALSE)
	d2lnL<- logLik(mle2,theta2,FALSE,FALSE,TRUE)
	C<-contrast(mle2,theta2,TRUE,FALSE,FALSE)
	dC<-contrast(mle2,theta2,FALSE,TRUE,FALSE)
	d2C<-contrast(mle2,theta2,FALSE,FALSE,TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("OthersNames + Weibull + CM GQR_ARAm-m=1-fun=log vs GQR_ARA1 + PM ARA1 + mutlisystems",{
	simData<-data.frame(S=c(rep(1,4),rep(2,4),rep(3,14)),T=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),C=c(1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0),row.names=1:22)
    mle <- mle.vam(S & T & C ~ (GQR_ARAm(0.7,1.3|1,log) | Weibull(0.001,2.5)) & (ARA1(0.6)),data=simData)
    mle2 <- mle.vam(S & T & C ~ (GQR_ARA1(0.7,1.3|log) | Weibull(0.001,2.5)) & (ARAm(0.6|1)),data=simData)
    theta<-c(0.03,2.4,0.7,1.3,0.8)
    theta2<-c(0.03,2.4,0.7,1.3,0.8)
	lnL<- logLik(mle2,theta2,TRUE,FALSE,FALSE)
	dlnL<- logLik(mle2,theta2,FALSE,TRUE,FALSE)
	d2lnL<- logLik(mle2,theta2,FALSE,FALSE,TRUE)
	C<-contrast(mle2,theta2,TRUE,FALSE,FALSE)
	dC<-contrast(mle2,theta2,FALSE,TRUE,FALSE)
	d2C<-contrast(mle2,theta2,FALSE,FALSE,TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM GQR-sqrt+PM GQR-ARA4-log",{
	simData<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,15)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78,15.2),Type=c(1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,2,1,-1,1,-1,0),row.names=1:23)
    theta<-c(0.3,2.3,0.9,1.3,0.7)
    mle <- mle.vam(System & Time & Type ~ (GQR(0.8|sqrt) | Weibull(0.001,2.5)) & (GQR_ARAm(1.2,0.9|log,2)+AGAN()),data=simData)
	lnL<- -19.7677532956276
	dlnL<- c(-43.6893072126687,-10.1401676628137,-20.0344727068032,-9.6021635142892,8.1377928954614)
	d2lnL<- matrix(c(-88.8888888888889,-52.735805905841,-110.620124025887,-73.2064165088806,57.1613591346954,-52.735805905841,-15.4626620386603,-26.0911341778741,-16.3124727895093,15.1653437915908,-110.620124025887,-26.0911341778741,-94.8861688132064,-65.3159909454354,28.0973156129578,-73.2064165088806,-16.3124727895093,-65.3159909454354,-36.5511188533409,23.2799421642019,57.1613591346954,15.1653437915908,28.0973156129578,23.2799421642019,-29.2947194485407),nrow=5,byrow=TRUE)
	C<- -14.4221879358273
	dC<-c(-0.315881117028209,0.573229735462233,4.03564481004899,-2.51094114471752)
	d2C<-matrix(c(-2.30514244145999,3.08976784976177,3.17804451035529,-3.42796096717849,3.08976784976177,-25.2616239090908,-11.6684468915397,0.430198562578086,3.17804451035529,-11.6684468915397,-10.894975844016,3.73273315713559,-3.42796096717849,0.430198562578086,3.73273315713559,-8.79105571582819),nrow=4,byrow=TRUE)
	expect_that(logLik(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ABAO + LeftCens",{
	simData<-data.frame(Time=c(2,3.36,4.04,4.97,5.16),Type=c(0,-1,-1,-1,-1),row.names=1:5)
    mle <- mle.vam(Time & Type ~ (ABAO() | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.3,2.5)
	lnL<- -8.81937993286691
	dlnL<- c(-41.4915492379502,-21.1451564444401)
	d2lnL<- matrix(c(-44.4444444444444,-95.3256617666342,-95.3256617666342,-48.681903539105),nrow=2,byrow=TRUE)
	C<- -2.02742310937457
	dC<-c(0.497622545672179)
	d2C<- -0.230996044748901
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
	simData2<-data.frame(Time=c(3.36,4.04,4.97,5.16),Type=c(-1,-1,-1,-1),row.names=1:4)
    mle2 <- mle.vam(Time & Type ~ (ABAO() | Weibull(0.001,2.5)),data=simData2)
    update(mle2,simData)
    expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(logLik(mle2,theta,TRUE,FALSE,FALSE),tolerance=0.00000000000001))
    mle3 <- mle.vam(Time & Type ~ (ABAO() | Weibull(0.001,2.5)),data=simData2)
    mle4 <- mle.vam(Time & Type ~ (ABAO() | Weibull(0.001,2.5)),data=simData)
    update(mle4,simData2)
    expect_that(logLik(mle3,theta,TRUE,FALSE,FALSE),equals(logLik(mle4,theta,TRUE,FALSE,FALSE),tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ABAO + LeftCens + RightCens",{
	simData<-data.frame(Time=c(2,3.36,4.04,4.97,5.16),Type=c(0,-1,-1,-1,0),row.names=1:5)
    mle <- mle.vam(Time & Type ~ (ABAO() | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.3,2.5)
	lnL<- -10.9931027296553
	dlnL<- c(-44.8248825712835,-23.186093023933)
	d2lnL<- matrix(c(-33.3333333333333,-95.3256617666342,-95.3256617666342,-48.521903539105),nrow=2,byrow=TRUE)
	C<- -2.65031513654516
	dC<-c(0.195415851148208)
	d2C<- -0.173247033561675
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
	simData2<-data.frame(Time=c(3.36,4.04,4.97,5.16),Type=c(-1,-1,-1,0),row.names=1:4)
    mle2 <- mle.vam(Time & Type ~ (ABAO() | Weibull(0.001,2.5)),data=simData2)
    update(mle2,simData)
    expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(logLik(mle2,theta,TRUE,FALSE,FALSE),tolerance=0.00000000000001))
    mle3 <- mle.vam(Time & Type ~ (ABAO() | Weibull(0.001,2.5)),data=simData2)
    mle4 <- mle.vam(Time & Type ~ (ABAO() | Weibull(0.001,2.5)),data=simData)
    update(mle4,simData2)
    expect_that(logLik(mle3,theta,TRUE,FALSE,FALSE),equals(logLik(mle4,theta,TRUE,FALSE,FALSE),tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ABAO + PM ARAInf + mutlisystems + LeftCens",{
	simData<-data.frame(System=c(rep(1,5),rep(2,4),rep(3,5),rep(4,4),rep(5,8)),Time=c(2.28,3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,1.57,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 1.95,2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,0,1, -1,-1,-1,0, 1,0,-1,-1,1, -1,1,1,0, 1,0,-1,1,-1,-1,1,0),row.names=1:26)
    mle <- mle.vam(System & Time & Type ~ (ABAO() | Weibull(0.001,2.5)) & (ARAInf(0.5)),data=simData)
    theta<-c(0.3,1.8,0.6)
	lnL<- -18.3495730327467
	dlnL<- c(-30.5997998851721,-14.4371097829312,5.2616319836793)
	d2lnL<- matrix(c(-100,-87.1905417248509,32.4548475547877,-87.1905417248509,-39.5521455750697,18.2435580647309,32.4548475547877,18.2435580647309,-7.82756357635268),nrow=3,byrow=TRUE)
	C<- -15.4974809498401
	dC<-c(-1.22908038160137,0.345220498680633)
	d2C<- matrix(c(-2.35185397287521,-0.727993035616652,-0.727993035616652,-2.76405321485956),nrow=2,byrow=TRUE)
	expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,TRUE,FALSE,FALSE),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,TRUE,FALSE),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,FALSE,FALSE,TRUE),equals(d2C,tolerance=0.00000000000001))
	simData2<-data.frame(System=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,8)),Time=c(2.28,3.36,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 1.95,2.45,2.78,3.56,4.23,5.32,6.43,6.98),Type=c(1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,0,-1,1,-1,-1,1,0),row.names=1:24)
	mle2 <- mle.vam(System & Time & Type ~ (ABAO() | Weibull(0.001,2.5)) & (ARAInf(0.5)),data=simData2)
    update(mle2,simData)
    expect_that(logLik(mle,theta,TRUE,FALSE,FALSE),equals(logLik(mle2,theta,TRUE,FALSE,FALSE),tolerance=0.00000000000001))
    mle3 <- mle.vam(System & Time & Type ~ (ABAO() | Weibull(0.001,2.5)) & (ARAInf(0.5)),data=simData2)
	mle4 <- mle.vam(System & Time & Type ~ (ABAO() | Weibull(0.001,2.5)) & (ARAInf(0.5)),data=simData)
    update(mle4,simData2)
    expect_that(logLik(mle3,theta,TRUE,FALSE,FALSE),equals(logLik(mle4,theta,TRUE,FALSE,FALSE),tolerance=0.00000000000001))
}
)

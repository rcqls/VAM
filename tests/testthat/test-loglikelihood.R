
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
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
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
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
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
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
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
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
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
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
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
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + multiSystems",{
	simData<-data.frame(System=c(1,2),Time=c(3.36,2.34),Type=c(-1,-1),row.names=1:2)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.3,0.8,0.6)
    logLikelihood(mle,theta,c(TRUE,FALSE,FALSE))
    logLikelihood(mle,theta,c(FALSE,TRUE,FALSE))
    logLikelihood(mle,theta,c(FALSE,FALSE,TRUE))
    contrast(mle,theta,c(TRUE,FALSE,FALSE))
    contrast(mle,theta,c(FALSE,TRUE,FALSE))
    contrast(mle,theta,c(FALSE,FALSE,TRUE))
    theta<-c(0.3,1.8,0.6)
	lnL<- -3.62638734326211
	dlnL<- c(-6.81229529030084,-1.22612955802,0)
	d2lnL<- matrix(c(-22.2222222222222,-14.6644419082528,0,-14.6644419082528,-5.52276726509835,0,0,0,0),nrow=3,byrow=TRUE)
	C<- -2.9907189791188
	dC<-c(0.99730223588401,0)
	d2C<-matrix(c(-0.67625373876488,0,0,0),nrow=2,byrow=TRUE)
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + CM ARAInf + Censorship + multiSystems",{
	simData<-data.frame(System=c(1,1,1,1,2,2,2,3),Time=c(3.36,4.04,4.97,5.16,2.34,3.46,5.02,4),Type=c(-1,-1,-1,0,-1,-1,-1,0),row.names=1:8)
    mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.3,0.8,0.6)
    logLikelihood(mle,theta,c(TRUE,FALSE,FALSE))
    logLikelihood(mle,theta,c(FALSE,TRUE,FALSE))
    logLikelihood(mle,theta,c(FALSE,FALSE,TRUE))
    contrast(mle,theta,c(TRUE,FALSE,FALSE))
    contrast(mle,theta,c(FALSE,TRUE,FALSE))
    contrast(mle,theta,c(FALSE,FALSE,TRUE))
    theta<-c(0.3,1.8,0.6)
	lnL<- -10.7518360550876
	dlnL<- c(-16.6082469606206,-4.41424143581748,0.565291823677669)
	d2lnL<- matrix(c(-66.6666666666667,-42.1908945206754,17.8580795012604,-42.1908945206754,-16.0970944508072,5.49620639347231,17.8580795012604,5.49620639347231,-6.72216164111632),nrow=3,byrow=TRUE)
	C<- -9.39660957648529
	dC<-c(1.32804496573541,-1.86523751785326)
	d2C<-matrix(c(-1.66488790656225,-3.08811044563958,-3.08811044563958,-3.53294273656389),nrow=2,byrow=TRUE)
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
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
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
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
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
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
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
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
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
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
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("Weibull + AGAN",{
	simData<-data.frame(Time=c(3.36,4.04,4.97,5.16),Type=c(-1,-1,-1,-1),row.names=1:4)
    mle <- mle.vam(Time & Type ~ (AGAN() | Weibull(0.001,2.5)),data=simData)
    theta<-c(0.3,0.8,1)
	lnL<- -6.90098251895707
	dlnL<- c(8.75359407434948,3.3717767832448)
	d2lnL<- matrix(c(-44.4444444444444,-2.40399936753948,-2.40399936753948,-7.6652713149061),nrow=2,byrow=TRUE)
	C<- -5.25256034408912
	dC<-1.99329429144225
	d2C<--9.26821752088532
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
}
)

test_that("LogLinear + CM ABAO",{
    simData<-data.frame(Time=c(3.36,4.04,4.97,5.16),Type=c(-1,-1,-1,-1),row.names=1:4)
    mle <- mle.vam(Time & Type ~ (ABAO() | LogLinear(0.001,2.5)),data=simData)
    theta<-c(0.3,0.8,1)
	lnL<- -13.6870254780068
	dlnL<- c(-62.9837808690102,-73.9249749593491)
	d2lnL<- matrix(c(-44.4444444444444,-304.849916531164,-304.849916531164,-390.943849373404),nrow=2,byrow=TRUE)
	C<- -1.77041141396439
	dC<-1.55193690275558
	d2C<--4.47702275378921
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
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
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
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
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
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
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
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
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
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
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
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
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
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
	expect_that(logLikelihood(mle,theta,c(TRUE,FALSE,FALSE)),equals(lnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,TRUE,FALSE)),equals(dlnL,tolerance=0.00000000000001))
	expect_that(logLikelihood(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2lnL,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(TRUE,FALSE,FALSE)),equals(C,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,TRUE,FALSE)),equals(dC,tolerance=0.00000000000001))
	expect_that(contrast(mle,theta,c(FALSE,FALSE,TRUE)),equals(d2C,tolerance=0.00000000000001))
}
)
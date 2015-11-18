#mleMultiCppTest.R
require(VAM)

nExp<-rep(10,1000)
#nExp <-10000
cat("Simulating...\n")
simulate(simCpp,nExp,as.list=length(nExp)>1) -> simDf
cat("Table of number of system:\n")
print(table(nExp))
update(mleCppMulti,data=simDf)
print(coef(mleCppMulti,par0)->res)
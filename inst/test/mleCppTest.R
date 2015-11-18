#mleCppTest.R
require(VAM)

nExp <-10000
cat("Simulating...\n")
simulate(simCpp,nExp,as.list=length(nExp)>1) -> simDf
cat("Number of events:",nExp,"\n")
update(mleCpp,data=simDf)

print(coef(mleCpp,par0)->res)
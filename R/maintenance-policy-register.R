# Register maintenance policies (but can be located in any other R file)
Periodic.maintenance.policy <-function(by,from=0,prob=1) {}
AtIntensity.maintenance.policy <- function(level=1,model=NULL) {}
AtVirtualAge.maintenance.policy <- function(level=1,,model=NULL) {} 
AtFailureProbability.maintenance.policy <- function(level=0.5,model=NULL) {}


########################### Do not consider this. It is just to remember!
## This detects automatically (like what proposed inside convert.mp inside vam.R) the args of a call Periodic.maintenance.policy(1,prob=c(.5,.5)) returns list(by=1,prob=c(.5,.5))
# pars <- as.list(match.call(
# 	switch(deparse(mp[[1]]),
# 		Periodic.maintenance.policy=function(by,from=0,prob=1) {},
# 		AtIntensity.maintenance.policy=, 
# 		AtVirtualAge.maintenance.policy=, 
# 		AtFailureProbability.maintenance.policy=function(level=1) {}
# 	),mp)
# )[-1]
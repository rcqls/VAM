#include "rcpp_maintenance_policy.h"

// Constructor from R to build different kind of 
MaintenancePolicy* newMaintenancePolicy(List policy) {
	std::string name=policy["name"];
	MaintenancePolicy*  mp=NULL; //default "None"
	//DEBUG:printf("name=%s\n",name.c_str());
	if(name.compare("Periodic.maintenance.policy") == 0) {
		//DEBUG:printf("Params:alpha=%lf,beta=%lf\n",alpha,beta);
		List pars=policy["params"];
		mp=new PeriodicMaintenancePolicy(pars);
	} //else if(name.compare("None") == 0) {
	//	mp=NULL;
	//}	
	return mp;
}
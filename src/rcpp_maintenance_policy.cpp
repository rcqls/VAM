#include "rcpp_maintenance_policy.h"
#include "rcpp_maintenance_model.h"

// Constructor from R to build different kind of 
MaintenancePolicy* newMaintenancePolicy(List policy) {
	std::string name=policy["name"];
	MaintenancePolicy*  mp=NULL; //default "None"
	//DEBUG:printf("name=%s\n",name.c_str());
	if(name.compare("Periodic.maintenance.policy") == 0) {
		//DEBUG:printf("Params:alpha=%lf,beta=%lf\n",alpha,beta);
		List pars=policy["params"];
		mp=new PeriodicMaintenancePolicy(pars);
	} else if(name.compare("AtIntensity.maintenance.policy") == 0) {
		//DEBUG:printf("Params:alpha=%lf,beta=%lf\n",alpha,beta);
		List pars=policy["params"];
		mp=new AtIntensityMaintenancePolicy(pars);
	} else if(name.compare("AtVirtualAge.maintenance.policy") == 0) {
		//DEBUG:printf("Params:alpha=%lf,beta=%lf\n",alpha,beta);
		List pars=policy["params"];
		mp=new AtVirtualAgeMaintenancePolicy(pars);
	} else if(name.compare("AtFailureProbability.maintenance.policy") == 0) {
		//DEBUG:printf("Params:alpha=%lf,beta=%lf\n",alpha,beta);
		List pars=policy["params"];
		mp=new AtFailureProbabilityMaintenancePolicy(pars);
	}
	return mp;
}

 
List PeriodicMaintenancePolicy::update(VamModel* model) {
    double current=model->time[model->k];
    // List update(double current) {   
	Function sample_int = Environment::base_env()["sample.int"];
	List res;
	res["time"] = from + (floor((current - from)/by) + 1) * by;
	//First argument not automatically wrapped in RcppWin64bits ??? => weird!
	res["type"]=sample_int(NumericVector::create(prob.size()),1,true,prob);
	return res;
 
};

List AtIntensityMaintenancePolicy::update(VamModel* model) {
    Function sample_int = Environment::base_env()["sample.int"];
    List res;
    
    res["time"] = model->models->at(model->idMod)->virtual_age_inverse(model->family->inverse_density(level[0]));
    //First argument not automatically wrapped in RcppWin64bits 
    res["type"]= 1; //sample_int(NumericVector::create(prob.size()),1,true,prob);
    return res;
};

List AtVirtualAgeMaintenancePolicy::update(VamModel* model) {
    Function sample_int = Environment::base_env()["sample.int"];
    List res;
    
    res["time"] = model->models->at(model->idMod)->virtual_age_inverse(level[0]);
    //First argument not automatically wrapped in RcppWin64bits 
    res["type"]= 1; //sample_int(NumericVector::create(prob.size()),1,true,prob);
    return res;
};

List AtFailureProbabilityMaintenancePolicy::update(VamModel* model) {
    Function sample_int = Environment::base_env()["sample.int"];
    List res;
    
    res["time"] = model->models->at(model->idMod)->virtual_age_inverse(model->family->inverse_cumulative_density(model->family->cumulative_density(model->models->at(model->idMod)->virtual_age(model->time[model->k]))-log(1-level)[0]));
    //First argument not automatically wrapped in RcppWin64bits 
    res["type"]= 1; //sample_int(NumericVector::create(prob.size()),1,true,prob);
    return res;
};



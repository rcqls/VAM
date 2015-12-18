#include "rcpp_maintenance_policy.h"
#include "rcpp_maintenance_model.h"
#include "rcpp_vam_module.h"

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
        if(as<bool>(policy["with.model"])) {
            VamModel* vmod_= as<VamModel*>(policy["model"]);
            vmod_->idMod=0;
            mp->set_vmod(vmod_);
        }
	} else if(name.compare("AtVirtualAge.maintenance.policy") == 0) {
		//DEBUG:printf("Params:alpha=%lf,beta=%lf\n",alpha,beta);
		List pars=policy["params"];
		mp=new AtVirtualAgeMaintenancePolicy(pars);
	} else if(name.compare("AtFailureProbability.maintenance.policy") == 0) {
		//DEBUG:printf("Params:alpha=%lf,beta=%lf\n",alpha,beta);
		List pars=policy["params"];
		mp=new AtFailureProbabilityMaintenancePolicy(pars);
	} else if(name.compare("MaintenancePolicyList") == 0) {
		List policies=policy["policies"];
		mp=new MaintenancePolicyList(policies);
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
	int t=as<int>(sample_int(NumericVector::create(prob.size()),1,true,prob))+get_from_type();
	res["type"]=t;
	return res;
 
};

List AtIntensityMaintenancePolicy::update(VamModel* model) {
    Function sample_int = Environment::base_env()["sample.int"];
    List res;
    VamModel* mod;

    if(get_vmod() != NULL) {
        mod = get_vmod();
        //update k time and everything needed to compute the next
        mod->k = model->k;
        mod->time = model->time;
        mod->update_Vleft(false);
        //TODO: like Simulation but to check HERE or 1 line before!
        mod->idMod = model->idMod; 
        mod->models->at(mod->idMod)->update(false);
    } else mod=model;
    
    //printf("at=%d\n",model->idMod);
    res["time"] = mod->models->at(mod->idMod)->virtual_age_inverse(mod->family->inverse_density(level[0]));
    //First argument not automatically wrapped in RcppWin64bits 
    res["type"]= 1+get_from_type(); //sample_int(NumericVector::create(prob.size()),1,true,prob);
    return res;
};

List AtVirtualAgeMaintenancePolicy::update(VamModel* model) {
    Function sample_int = Environment::base_env()["sample.int"];
    List res;
    
    res["time"] = model->models->at(model->idMod)->virtual_age_inverse(level[0]);
    //First argument not automatically wrapped in RcppWin64bits 
    res["type"]= 1+get_from_type(); //sample_int(NumericVector::create(prob.size()),1,true,prob);
    return res;
};

List AtFailureProbabilityMaintenancePolicy::update(VamModel* model) {
    Function sample_int = Environment::base_env()["sample.int"];
    List res;
    
    res["time"] = model->models->at(model->idMod)->virtual_age_inverse(model->family->inverse_cumulative_density(model->family->cumulative_density(model->models->at(model->idMod)->virtual_age(model->time[model->k]))-log(1-level)[0]));
    //First argument not automatically wrapped in RcppWin64bits 
    res["type"]= 1+get_from_type(); //sample_int(NumericVector::create(prob.size()),1,true,prob);
    return res;
};

MaintenancePolicyList::MaintenancePolicyList(List policies) {
	int ft=0;
	//int i=0;
    for(
        List::iterator lit=policies.begin();
        lit != policies.end();
        ++lit //,i++
    ) {
    	//printf("i=%d\n",i);
    	List policy=*lit;
    	MaintenancePolicy*  mp=newMaintenancePolicy(policy);
        if(!(mp == NULL)) {
        	//printf("i=%d\n",i);
        	mp->set_from_type(ft);
            policy_list.push_back(mp);
        	ft += mp -> type_size(); //update ft for next mp
        }
    }
};

MaintenancePolicyList::~MaintenancePolicyList() {
	for(
		std::vector<MaintenancePolicy*>::iterator vit=policy_list.begin();
		vit != policy_list.end();
        ++vit
    ) {
		delete *vit;
	}

};


List MaintenancePolicyList::update(VamModel* model) {
    List res,res2;
    res=at(0)->update(model);
    for(int i=1;i<size();i++) {
    	res2 = at(i)->update(model);
    	if(as<float>(res2["time"]) < as<float>(res["time"])) res=res2;
    }
    return res;
};





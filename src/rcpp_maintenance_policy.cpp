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
        VamModel* ext_mod_= as<VamModel*>(policy["model"]);
        ext_mod_->idMod=0;
        mp->set_external_model(ext_mod_);
    }
	} else if(name.compare("AtVirtualAge.maintenance.policy") == 0) {
		//DEBUG:printf("Params:alpha=%lf,beta=%lf\n",alpha,beta);
		List pars=policy["params"];
		mp=new AtVirtualAgeMaintenancePolicy(pars);
		if(as<bool>(policy["with.model"])) {
				VamModel* ext_mod_= as<VamModel*>(policy["model"]);
				ext_mod_->idMod=0;
				mp->set_external_model(ext_mod_);
		}
	} else if(name.compare("AtFailureProbability.maintenance.policy") == 0) {
		//DEBUG:printf("Params:alpha=%lf,beta=%lf\n",alpha,beta);
		List pars=policy["params"];
		mp=new AtFailureProbabilityMaintenancePolicy(pars);
		if(as<bool>(policy["with.model"])) {
				VamModel* ext_mod_= as<VamModel*>(policy["model"]);
				ext_mod_->idMod=0;
				mp->set_external_model(ext_mod_);
		}
	} else if(name.compare("MaintenancePolicyList") == 0) {
		List policies=policy["policies"];
		mp=new MaintenancePolicyList(policies);
	}
	return mp;
}

VamModel* MaintenancePolicy::update_external_model(VamModel* model) {
    VamModel* mod;

    if(get_external_model() != NULL) {
        mod = get_external_model();
        //update everything needed to compute the update step (see end of simulation step)
        mod->k = model->k;
        mod->time = model->time;//Normally, it is a pointer: it would be better to do it once but it is not such cost.
        mod->type = model->type;//Idem
        if(mod->k > 0) {//Not at the init step!
            //VERY VERY IMPORTANT: Since we want to update at previous step
            mod->k -= 1;

            //And then as what it is applied at the end of the simulation step but for mod instead of model!
            mod->update_Vleft(false,false); //mod->idMod is not yet updated!
            mod->idMod = model->idMod; //mod->idMod is then updated for the next task!
            mod->models->at(mod->idMod)->update(false,false);
        }
    } else mod=model;
    return mod;
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
		//printf("ici\n");
    VamModel* mod=update_external_model(model);
		//ToRemove: double next_time2=model->virtual_age_inverse(model->family->inverse_hazardRate(level[0]));
		double next_time=mod->virtual_age_inverse(mod->family->inverse_hazardRate(level[0]));
		//ToRemove: printf("next_time=(%lf,%lf,%lf,%lf,%lf,%lf)\n",next_time,next_time2,mod->Vright,model->Vright,mod->C,model->C );
		//printf("at=%d\n",model->idMod);
    res["time"] = next_time ;//mod->virtual_age_inverse(mod->family->inverse_hazardRate(level[0]));

    res["type"]= 1+get_from_type(); //sample_int(NumericVector::create(prob.size()),1,true,prob);
    return res;
};

void AtIntensityMaintenancePolicy::first() {
	//Needs init for external model for new simulation
	// printf("mp: %p\n",model->maintenance_policy);
	// printf("mp: %p\n",(model->maintenance_policy->get_external_model()));
	if( get_external_model() != NULL) get_external_model()->init_computation_values();
}

List AtVirtualAgeMaintenancePolicy::update(VamModel* model) {
    Function sample_int = Environment::base_env()["sample.int"];
    List res;
    VamModel* mod=update_external_model(model);

    res["time"] = mod->virtual_age_inverse(level[0]);

    res["type"]= 1+get_from_type(); //sample_int(NumericVector::create(prob.size()),1,true,prob);
    return res;
};

void AtVirtualAgeMaintenancePolicy::first() {
	//Needs init for external model for new simulation
	// printf("mp: %p\n",model->maintenance_policy);
	// printf("mp: %p\n",(model->maintenance_policy->get_external_model()));
	if(get_external_model() != NULL) get_external_model()->init_computation_values();
}

List AtFailureProbabilityMaintenancePolicy::update(VamModel* model) {
    Function sample_int = Environment::base_env()["sample.int"];
    List res;
    VamModel* mod=update_external_model(model);

    res["time"] = mod->virtual_age_inverse(mod->family->inverse_cumulative_hazardRate(mod->family->cumulative_hazardRate(mod->virtual_age(mod->time[mod->k]))-log(1-level)[0]));
    //First argument not automatically wrapped in RcppWin64bits
    res["type"]= 1+get_from_type(); //sample_int(NumericVector::create(prob.size()),1,true,prob);
    return res;
};

void AtFailureProbabilityMaintenancePolicy::first() {
	//Needs init for external model for new simulation
	// printf("mp: %p\n",model->maintenance_policy);
	// printf("mp: %p\n",(model->maintenance_policy->get_external_model()));
	if( get_external_model() != NULL) get_external_model()->init_computation_values();
}

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

void MaintenancePolicyList::first() {

	for(
		std::vector<MaintenancePolicy*>::iterator vit=policy_list.begin();
		vit != policy_list.end();
        ++vit
    ) {
			(*vit)->first();
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

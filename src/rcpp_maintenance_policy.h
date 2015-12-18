#ifndef RCPP_MAINTENANCE_POLICY_H
#define RCPP_MAINTENANCE_POLICY_H
#include <Rcpp.h>
#include "rcpp_maintenance_policy.h"
#include "rcpp_vam_model.h"

using namespace Rcpp ;
class VamModel;

class MaintenancePolicy {
public:

	MaintenancePolicy() {};

	virtual ~MaintenancePolicy() {};

	virtual List update(VamModel* model) = 0;
    //virtual List update(double current) = 0;

	virtual List get_params() = 0;

	virtual void set_params(List params) = 0;

    void set_vmod(VamModel* vmod_) {
        vmod=vmod_;
    }

    VamModel* get_vmod() {
        return vmod;
    }

    void set_from_type(int from_type_) {
        from_type=from_type_;
    };

    int get_from_type() {
        return from_type;
    };

    virtual int type_size() = 0;


private:

    int from_type;

    VamModel* vmod;

};

//IMPORTANT: Names of params of R has to be registered in maintenance-policy-register.R or in any other R file (if you think about plugin)  
// Ex: AtIntensity(level=1.2)
// Since in vam.R (inside convert.mp),  "level" is defined by default at 0.5, 
// this allow us to magically call AtIntensity(1.2) instead of AtIntensity(level=1.2). 
//THEN, if you want to provide this feature for other parameters you could set the default
// inside convert.mp of file vam.R

class PeriodicMaintenancePolicy : public MaintenancePolicy {
public:
    PeriodicMaintenancePolicy(List params) {
    	set_params(params);
        set_from_type(0);
    }

    ~PeriodicMaintenancePolicy() {};

    NumericVector from, by, prob;

    List get_params() {
    	List out;
    	out["from"]=from;out["by"]=by;out["prob"]=prob;
    	return out;
    }

    void set_params(List params) {
    	from=params["from"];by=params["by"];prob=params["prob"];
    }

    int type_size() {
        return prob.size();
    }

    List update(VamModel* model);
 
};

class AtIntensityMaintenancePolicy : public MaintenancePolicy {
public:
    AtIntensityMaintenancePolicy(List params) {
        set_params(params);
        set_from_type(0);
        set_vmod(NULL);
    }

    ~AtIntensityMaintenancePolicy() {};

    NumericVector level;

    List get_params() {
        List out;
        out["level"]=level;
        return out;
    }

    void set_params(List params) {
        level=params["level"];
    }

    int type_size() {
        return 1;
    }

    List update(VamModel* model);
 
};

class AtVirtualAgeMaintenancePolicy : public MaintenancePolicy {
public:
    AtVirtualAgeMaintenancePolicy(List params) {
        set_params(params);
        set_from_type(0);
    }

    ~AtVirtualAgeMaintenancePolicy() {};

    NumericVector level;

    List get_params() {
        List out;
        out["level"]=level;
        return out;
    }

    void set_params(List params) {
        level=params["level"];
    }

    int type_size() {
        return 1;
    }

    List update(VamModel* model);
 
};

class AtFailureProbabilityMaintenancePolicy : public MaintenancePolicy {
public:
    AtFailureProbabilityMaintenancePolicy(List params) {
        set_params(params);
        set_from_type(0);
    }

    ~AtFailureProbabilityMaintenancePolicy() {};

    NumericVector level;

    List get_params() {
        List out;
        out["level"]=level;
        return out;
    }

    void set_params(List params) {
        level=params["level"];
    }

    int type_size() {
        return 1;
    }

    List update(VamModel* model);
 
};

class MaintenancePolicyList : public MaintenancePolicy {//List of MaintenancePolicy (heterogeneous terms) 
public:
    MaintenancePolicyList(List policies);

    ~MaintenancePolicyList();

    List get_params() {
        List out;
        //run over all policies and get_params
        for(int i=0;i<size();i++) {
            out[i]=at(i)->get_params();
        }
        return out;
    }

    void set_params(List params) {//params is here a List of List
        //run over all policies and set_params
        for(int i=0;i<size();i++) {
            at(i)->set_params(params[i]);
        }
    }

    int type_size() {
        int s=0;
        for(int i=0;i<size();i++) {
            s += at(i)->type_size();
        }
        return s;
    }

    List update(VamModel* model);

    //Further ones because of list of policy

    MaintenancePolicy* at(int i) {
        return policy_list[i];
    }

    int size() {
        return policy_list.size();
    }


protected:

    std::vector<MaintenancePolicy*> policy_list; //policy list
     
};

MaintenancePolicy* newMaintenancePolicy(List policy);

#endif //RCPP_FAMILY_MODEL_H
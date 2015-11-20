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

    List update(VamModel* model);
 
};

class AtIntensityMaintenancePolicy : public MaintenancePolicy {
public:
    AtIntensityMaintenancePolicy(List params) {
        set_params(params);
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

    List update(VamModel* model);
 
};

class AtVirtualAgeMaintenancePolicy : public MaintenancePolicy {
public:
    AtVirtualAgeMaintenancePolicy(List params) {
        set_params(params);
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

    List update(VamModel* model);
 
};

class AtFailureProbabilityMaintenancePolicy : public MaintenancePolicy {
public:
    AtFailureProbabilityMaintenancePolicy(List params) {
        set_params(params);
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

    List update(VamModel* model);
 
};

MaintenancePolicy* newMaintenancePolicy(List policy);

#endif //RCPP_FAMILY_MODEL_H
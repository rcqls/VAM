#ifndef RCPP_MAINTENANCE_POLICY_H
#define RCPP_MAINTENANCE_POLICY_H
#include <Rcpp.h>

using namespace Rcpp ;

class MaintenancePolicy {
public:

	MaintenancePolicy() {};

	virtual ~MaintenancePolicy() {};

	virtual List update(double current) = 0;

	virtual List get_params() = 0;

	virtual void set_params(List params) = 0;

};

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

    List update(double current) {
  		Function sample_int = Environment::base_env()["sample.int"];
  		List res;
  		res["time"] = from + (floor((current - from)/by) + 1) * by;
        //First argument not automatically wrapped in RcppWin64bits 
		res["type"]=sample_int(NumericVector::create(prob.size()),1,true,prob);
		return res;
    }
 
};

MaintenancePolicy* newMaintenancePolicy(List policy);

#endif //RCPP_FAMILY_MODEL_H
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

    void set_external_model(VamModel* external_model_) {
        external_model=external_model_;
    }

    VamModel* get_external_model() {
        return external_model;
    }

    VamModel* update_external_model(VamModel* model);


    void set_from_type(int from_type_) {
        from_type=from_type_;
    };

    int get_from_type() {
        return from_type;
    };

    virtual int type_size() = 0;

		virtual void first() = 0;


private:

    int from_type;

    VamModel* external_model;

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

		void first() {};

};

class AtTimesMaintenancePolicy : public MaintenancePolicy {
public:
    AtTimesMaintenancePolicy(List params) {
        set_params(params);
        set_from_type(0);
    }

    ~AtTimesMaintenancePolicy() {};

    NumericVector times;
    int i,k;
    bool cycle;

    List get_params() {
        List out;
        out["times"]=NumericVector(times);out["cycle"]=cycle;
        return out;
    }

    void set_params(List params) {
        cycle=params["cycle"];
        times=params["times"];
    }

    int type_size() {
        return 1;
    }

    List update(VamModel* model);

        void first();

};

class AtIntensityMaintenancePolicy : public MaintenancePolicy {
public:
    AtIntensityMaintenancePolicy(List params) {
        set_params(params);
        set_from_type(0);
        set_external_model(NULL);
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

		void first();

};

class AtVirtualAgeMaintenancePolicy : public MaintenancePolicy {
public:
    AtVirtualAgeMaintenancePolicy(List params) {
        set_params(params);
        set_from_type(0);
        set_external_model(NULL);
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

		void first();

};

class AtFailureProbabilityMaintenancePolicy : public MaintenancePolicy {
public:
    AtFailureProbabilityMaintenancePolicy(List params) {
        set_params(params);
        set_from_type(0);
        set_external_model(NULL);
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

		void first();

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

		void first();


protected:

    std::vector<MaintenancePolicy*> policy_list; //policy list

};

MaintenancePolicy* newMaintenancePolicy(List policy);

#endif //RCPP_MAINTENANCE_POLICY_H

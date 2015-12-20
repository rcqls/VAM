#ifndef RCPP_SIM_VAM_H
#define RCPP_SIM_VAM_H
#include <Rcpp.h>
#include "rcpp_maintenance_model.h"
//#include "rcpp_stop_policy.h"

using namespace Rcpp ;

class StopPolicy;

class SimVam { 

public:

    SimVam(List model_) {
        model=new VamModel(model_);
    };

    ~SimVam() {
        delete model;
    };

    DataFrame get_data();
    
    DataFrame simulate(int nbsim);

    VamModel* get_model() {
        return model;
    }

    NumericVector get_params() {
        return model->get_params();
    }

    void set_params(NumericVector pars) {
        model->set_params(pars);
    }

    //delegate from model cache!
    List get_virtual_age_infos(double by) {
        return model->get_virtual_age_infos(by);
    }

    void add_stop_policy(List policy);

    int cache_size,size;
private:
    VamModel* model;

    StopPolicy* stop_policy;

    void init(int cache_size_);

    void resize();

};

#endif //RCPP_SIM_VAM_H
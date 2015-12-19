#ifndef RCPP_SIM_VAM_H
#define RCPP_SIM_VAM_H
#include <Rcpp.h>
#include "rcpp_maintenance_model.h"

using namespace Rcpp ;

class SimVam { 

public:

    SimVam(List model_) {
        model=new VamModel(model_);
    };

    ~SimVam() {
        delete model;
    };

    DataFrame get_data() {
        return DataFrame::create(_["Time"]=model->time,_["Type"]=model->type);
    }

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


private:
    VamModel* model;

    void init(int nbsim);

    void resize();

};

#endif //RCPP_SIM_VAM_H
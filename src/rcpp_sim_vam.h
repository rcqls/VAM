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

    //This is used inside simulate to fetch last generated data
    DataFrame get_last_data();

    //The 2 following functions are now similarly used in plot.R than mle.vam and model.vam
    //This is just a delagation to model 
    void set_data(List data_) {
        model->set_data(data_);
    }

    DataFrame get_selected_data(int i) {
        return model->get_selected_data(i);
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
    List get_virtual_age_infos(double by,double from, double to) {
        return model->get_virtual_age_infos(by,from,to);
    }

    void add_stop_policy(List policy);

    int cache_size,size;

    //Covariates related

    void select_current_system(int i) {
        model->select_current_system(i,false);
    }

	double compute_covariates() {
        return exp(-model->compute_covariates());
    }; 

    void set_covariates(List covariates_) {
        model->set_covariates(covariates_);
    }

private:
    VamModel* model;

    StopPolicy* stop_policy;

    void init(int cache_size_);

    void resize();

};

#endif //RCPP_SIM_VAM_H

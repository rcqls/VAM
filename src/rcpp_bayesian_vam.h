#ifndef RCPP_BAYESIAN_VAM_H
#define RCPP_BAYESIAN_VAM_H
#include <Rcpp.h>
#include "rcpp_mle_vam.h"

using namespace Rcpp ;

class BayesianVam {

public:

    BayesianVam(List model_,List data_) {
        mle=new MLEVam(model_,data_);
        model=mle->get_model();
    }

    ~BayesianVam() {
        //DEBUG: printf("BayesianVAM: %p, %p, %p\n",model,dS1,dS2);
        model=NULL;
        delete mle;
    };

    void set_data(List data_) {
        model->set_data(data_);
    }

    //NumericVector get_params() {
        //return model->get_params();
    //}

    void set_params(NumericVector pars) {
        //model->set_params(pars);
    }

    //List MCMC(int nb, int burn) {
    //
    //}

private:

    MLEVam *mle;

    VamModel *model;
};

#endif //RCPP_BAYESIAN_VAM_H

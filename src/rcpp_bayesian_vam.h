#ifndef RCPP_BAYESIAN_VAM_H
#define RCPP_BAYESIAN_VAM_H
#include <Rcpp.h>
#include "rcpp_mle_vam.h"
#include "rcpp_bayesian_prior.h"

using namespace Rcpp ;

class BayesianVam {

public:

    BayesianVam(List model_,List data_, List priors_) {
        mle=new MLEVam(model_,data_);
        model=mle->get_model();
        set_priors(priors_);
    }

    ~BayesianVam() {
        //DEBUG: printf("BayesianVAM: %p, %p, %p\n",model,dS1,dS2);
        model=NULL;
        delete mle;
        delete priors;
    };

    void set_priors(List priors_) {// Notice that priors are in provided in the same order as params (to fulfill set_params of VamModel)
        priors=new BayesianPriorList(priors_);
    }

    void set_data(List data_) {
        model->set_data(data_);
    }

    //NumericVector get_params() {
        //return model->get_params();
    //}

    void set_params(NumericVector pars) {
        //model->set_params(pars);
    }

    List mcmc(int nb, int burn) {
      return List::create();
    }

private:

    MLEVam *mle;

    VamModel *model;

    BayesianPriorList* priors;
};

#endif //RCPP_BAYESIAN_VAM_H

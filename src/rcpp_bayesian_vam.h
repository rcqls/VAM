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
      //initial values of sigma is 1
      for(int j=0;j<priors->size();j++) priors->at(j)->set_sigma(1);
    }

    //facility method to initialize directly in R sigma of the proposal
    void set_sigma(double sigma_, int j) {
      priors->at(j)->set_sigma(sigma_);
    }

    void set_data(List data_) {
      model->set_data(data_);
    }

    NumericVector get_params() {
        return model->get_params();
    }

    void set_params(NumericVector pars) {
        model->set_params(pars);
    }

    //Is it make sense to consider multisystem for Bayesian framework? I guess so since you can compute
    // the contrast with these data.

    List mcmc(NumericVector pars,int nb, int burn,bool alpha_fixed=false) {
      int nbParams = pars.size();
      BayesianPrior* curPrior;
      NumericVector curPars(pars);
      int i,j,jStart=0;
      double oldL,L,r,paramCurVal;
      oldL=mle->contrast(curPars,alpha_fixed)[0];
      if(alpha_fixed) jStart++;
      List res(nbParams-jStart);
      for(j=jStart;j<nbParams;j++) res[j]=List::create();
      for(i=0;i<burn;i++) {
        for(j=jStart;j<nbParams;j++) {
          curPrior=priors->at(j);
          //propose a new value!
          paramCurVal=curPrior->new_proposal(curPars[j]);
          // compute the ratio
          L=mle->contrast(curPars,alpha_fixed)[0];
          r=L/oldL*(curPrior->density(curPars[j]))/curPrior->density(paramCurVal);
          if(r > R::runif(0,1)) {
            curPrior->push_back(paramCurVal);
            curPars[j]=paramCurVal;
            oldL=L;
          }
        }
      }
      for(j=jStart;j<nbParams;j++) {
        std::vector<double> res2=priors->at(j)->get_result();
        NumericVector res3(res2.begin(),res2.end());
        res[j-jStart]=res2;
      }
      return res;
    }

private:

    MLEVam *mle;

    VamModel *model;

    BayesianPriorList* priors;
};

#endif //RCPP_BAYESIAN_VAM_H

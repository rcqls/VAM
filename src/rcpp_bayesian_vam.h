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
      for(int j=0;j<priors->size();j++) set_sigma_at(j,1);
    }

    //facility method to initialize directly in R sigma of the proposal
    void set_sigma_at( int j, double sigma_) {
      priors->at(j)->set_sigma(sigma_);
    }

    double get_sigma(int j) {
      return priors->at(j)->get_sigma();
    }

    void set_data(List data_) {
      model->set_data(data_);
    }

    NumericVector get_params() {
        return model->get_params();
    }

    NumericVector get_prior_params(int j) {
      if(j<priors->size()) {
        return priors->at(j)->get_params();
      } else {
        return NumericVector::create();
      }
    }

    void set_params(NumericVector pars) {
        model->set_params(pars);
    }

    //Is it make sense to consider multisystem for Bayesian framework? I guess so since you can compute
    // the contrast with these data.

    List mcmc(NumericVector pars,int nb, int burn, bool alpha_fixed=false) {
      int nbParams = pars.size();
      BayesianPrior* curPrior;
      NumericVector curPars=clone(pars),oldPars=clone(pars);
      int i,j,jStart=1;//jStart=1 here because alpha is not considered yet!
      double oldL,L,r,r0,paramCurVal;
      oldL=mle->contrast(oldPars,alpha_fixed)[0];
      //if(alpha_fixed) jStart++;
      List res(nbParams-jStart);
      for(j=jStart;j<nbParams;j++) priors->at(j)->clear();
      for(i=0;i<nb;i++) {
        for(j=jStart;j<nbParams;j++) {
          curPrior=priors->at(j);
          //propose a new value!
          paramCurVal=curPrior->new_proposal(oldPars[j]);
          curPars[j]=paramCurVal;
          // compute the ratio
          // for(int jj=0;jj<nbParams;jj++) {
          //   printf("[%d](%lf=%lf),",jj,oldPars[jj],curPars[jj]);
          // }
          // printf("\n");
          L=mle->contrast(curPars,alpha_fixed)[0];
          r=exp(L-oldL)*(curPrior->density(curPars[j]))/curPrior->density(oldPars[j]);
          r0=R::runif(0,1);
          //printf("r,r0=%lf,%lf\n",r,r0);
          if(r > r0) {
            if(i>burn) curPrior->push_back(curPars[j]);
            oldPars[j]=curPars[j];
            oldL=L;
          } else {
            //Put the old value;
            curPars[j]=oldPars[j];
          }
        }
      }
      // for(int jj=0;jj<nbParams;jj++) {
      //   printf("[%d](%lf=%lf),",jj,oldPars[jj],curPars[jj]);
      // }
      // L=mle->contrast(curPars,alpha_fixed)[0];
      // printf("here:%lf==%lf\n",L,oldL);
      // oldL=mle->contrast(oldPars,alpha_fixed)[0];
      // printf("here:%lf==%lf\n",L,oldL);
      // L=mle->contrast(curPars,alpha_fixed)[0];
      // printf("here:%lf==%lf\n",L,oldL);
      // oldL=mle->contrast(oldPars,alpha_fixed)[0];
      // printf("here:%lf==%lf\n",L,oldL);
      for(j=jStart;j<nbParams;j++) {
        std::vector<double> res2=priors->at(j)->get_result();
        //NumericVector res3(res2.begin(),res2.end());
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

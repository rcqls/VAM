#include "rcpp_bayesian_prior.h"

using namespace Rcpp ;


BayesianPriorList::BayesianPriorList(List priors_) {
    //int i=0;
    //int j=0;
    for(
        List::iterator lit=priors_.begin();
        lit != priors_.end();
        ++lit
    ) {
    	List prior=*lit;
    	BayesianPrior*  bp=newBayesianPrior(prior);
        if(!(bp == NULL)) {
            // bp->set_id(i++);
            // bp->set_id_params(j);
            // j+=bp->nb_params();
            prior_list.push_back(bp);
        }
    }
}

BayesianPriorList::~BayesianPriorList() {
	for(
		std::vector<BayesianPrior*>::iterator pit=prior_list.begin();
		pit != prior_list.end();
        ++pit
    ) {
		delete *pit;
	}

}


BayesianPrior* newBayesianPrior(List prior) {
	std::string name=prior["name"];
	NumericVector params=prior["params"];
	BayesianPrior*  bp=NULL;

	if(name.compare("Unif.prior") == 0) {
    bp=new UnifPrior(params);
	} else if(name.compare("Beta.prior") == 0) {
    bp=new BetaPrior(params);
  }

  return bp;
}

#ifndef RCPP_BAYESIAN_PRIOR_H
#define RCPP_BAYESIAN_PRIOR_H
#include <Rcpp.h>

using namespace Rcpp ;

//Forward declarations
class BayesianPrior;

//Effective declarations
class BayesianPriorList {//List of BayesianPrior (heterogeneous terms)
public:
    BayesianPriorList(List priors_);

    ~BayesianPriorList();

    BayesianPrior* at(int i) {
    	return prior_list[i];
    }

    int size() {
    	return prior_list.size();
    }


protected:

    std::vector<BayesianPrior*> prior_list; //prior list

};


class BayesianPrior {
public:
    BayesianPrior() {};

    virtual ~BayesianPrior() {};

    virtual double get() = 0;

};

class UnifPrior : public BayesianPrior {

public:

    UnifPrior(NumericVector params_) : BayesianPrior() {
    	a = params_[0];b=params_[1];
    }

    double get() {
      return R::runif(a,b);
    };

private:
    double a,b;

};

class BetaPrior : public BayesianPrior {

public:

    BetaPrior(NumericVector params_) : BayesianPrior() {
    	a = params_[0];b=params_[1];
    }

    double get() {
      return R::beta(a,b);
    };

private:
    double a,b;

};

BayesianPrior* newBayesianPrior(List prior);

#endif //RCPP_BAYESIAN_PRIOR_H

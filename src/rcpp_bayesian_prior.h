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

    virtual double density(double x) = 0;

    virtual NumericVector get_params() = 0;

    double new_proposal(double par) {
      return R::rnorm(par,sigma);
    }


    void set_sigma(double sigma_) {sigma=sigma_;}; //initialization could be managed directly in R

    double get_sigma() {return sigma;}; //initialization could be managed directly in R

    void clear() {result.clear();};

    void push_back(double res) {result.push_back(res);};

    std::vector<double> get_result() {return result;};

private:

    double sigma;

    std::vector<double> result;

};

class UnifPrior : public BayesianPrior {

public:

    UnifPrior(NumericVector params_) : BayesianPrior() {
    	a = params_[0];b=params_[1];
    }

    NumericVector get_params() {
      //printf("Unif:a,b=%lf,%lf\n",a,b);
      return NumericVector::create(a,b);
    }

    double get() {
      return R::runif(a,b);
    };

    double density(double x) {
      //printf("Unif:a,b=%lf,%lf\n",a,b);
      return R::dunif(x,a,b,0);
    };

private:
    double a,b;

};

class BetaPrior : public BayesianPrior {

public:

    BetaPrior(NumericVector params_) : BayesianPrior() {
    	a = params_[0];b = params_[1];
    }

    NumericVector get_params() {
      //printf("Beta:a,b=%lf,%lf\n",a,b);
      return NumericVector::create(a,b);
    }

    double get() {
      return R::beta(a,b);
    };

    double density(double x) {
      return R::dbeta(x,a,b,0);
    };

private:
    double a,b;

};

class GammaPrior : public BayesianPrior {

public:

    GammaPrior(NumericVector params_) : BayesianPrior() {
    	a = params_[0];s=params_[1];
    }

    NumericVector get_params() {
      //printf("Gamma:a,s=%lf,%lf\n",a,s);
      return NumericVector::create(a,s);
    }

    double get() {
      return R::rgamma(a,s);
    };

    double density(double x) {
      return R::dgamma(x,a,s,0);
    };

private:
    double a,s;

};

BayesianPrior* newBayesianPrior(List prior);

#endif //RCPP_BAYESIAN_PRIOR_H

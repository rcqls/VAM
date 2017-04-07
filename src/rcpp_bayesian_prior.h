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
      switch (proposal_distribution){
      case 1:
        return exp(R::rnorm(log(par)-log(1+pow(sigma/par,2.0))/2,pow(log(1+pow(sigma/par,2.0)),0.5)));
        break;
      default:
        return R::rnorm(par,sigma);
      }
    }

    bool not_symetric_proposal() {
      switch (proposal_distribution){
      case 1:
        return true;
        break;
      default:
        return false;
      }
    }

    double density_proposal(double x,double par) {
      switch (proposal_distribution){
      case 1:
        return (R::dnorm(log(x),log(par)-log(1+pow(sigma/par,2.0))/2,pow(log(1+pow(sigma/par,2.0)),0.5),false))/x;
        break;
      default:
        return R::dnorm(x,par,sigma,false);
      }
    }

    void set_proposal(int proposal_distribution_) {proposal_distribution=proposal_distribution_;}

    void set_sigma(double sigma_) {sigma=sigma_;}; //initialization could be managed directly in R

    double get_sigma() {return sigma;}; //initialization could be managed directly in R

    double get_proposal() {return proposal_distribution;}; //initialization could be managed directly in R

    void clear() {result.clear();};

    void push_back(double res) {result.push_back(res);};

    std::vector<double> get_result() {return result;};

private:

    double sigma;

    int proposal_distribution;//instrumental proposal distribution: 0 for normal and 1 for log-normal

    std::vector<double> result;

};

class NonInformativePrior : public BayesianPrior {

public:

    NonInformativePrior(NumericVector params_) : BayesianPrior() {
    }

    NumericVector get_params() {
      //printf("Unif:a,b=%lf,%lf\n",a,b);
      return NumericVector::create();
    }

    double get() {
      return R::runif(0,2);//Non sense!
    };

    double density(double x) {
      //printf("Unif:a,b=%lf,%lf\n",a,b);
      return 1/x;
    };

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
      return R::rbeta(a,b); //R::rbeta plutôt à la place R::beta!!!
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

class NormPrior : public BayesianPrior {

public:

    NormPrior(NumericVector params_) : BayesianPrior() {
    	m = params_[0];s=params_[1];
    }

    NumericVector get_params() {
      //printf("Gamma:a,s=%lf,%lf\n",a,s);
      return NumericVector::create(m,s);
    }

    double get() {
      return R::rnorm(m,s);
    };

    double density(double x) {
      return R::dnorm(x,m,s,0);
    };

private:
    double m,s;

};

class LNormPrior : public BayesianPrior {

public:

  LNormPrior(NumericVector params_) : BayesianPrior() {
    m = params_[0];s=params_[1];
  }

  NumericVector get_params() {
    //printf("Gamma:a,s=%lf,%lf\n",a,s);
    return NumericVector::create(m,s);
  }

  double get() {
    return exp(R::rnorm(m,s));
  };

  double density(double x) {
    return (R::dnorm(log(x),m,s,0))/x;
  };

private:
  double m,s;

};


BayesianPrior* newBayesianPrior(List prior);

#endif //RCPP_BAYESIAN_PRIOR_H

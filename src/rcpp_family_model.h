#ifndef RCPP_FAMILY_MODEL_H
#define RCPP_FAMILY_MODEL_H
#include <Rcpp.h> 

using namespace Rcpp ;

class FamilyModel {
public:

	FamilyModel() {};

	virtual ~FamilyModel() {};

	virtual NumericVector get_params() = 0;

	virtual void set_params(double, double) = 0;

	virtual double density(double x) = 0; 

	virtual double cumulative_density(double x) = 0;

	virtual double inverse_cumulative_density(double x) = 0;

	virtual double density_derivative(double x) = 0;

	virtual double density_param_derivative(double x) = 0;

	virtual double cumulative_density_param_derivative(double x) = 0;
};

class WeibullFamilyModel : public FamilyModel {
public:
    WeibullFamilyModel(double alpha_, double beta_) {
    	alpha=alpha_;beta=beta_;
    }

    ~WeibullFamilyModel() {};

    double alpha, beta;

    NumericVector get_params() {
    	NumericVector out(2);
    	out[0]=alpha;out[1]=beta;
    	return out;
    }

    void set_params(double alpha_, double beta_) {
    		alpha=alpha_;beta=beta_;
    }

    double density(double x) {
    	return (x<=0 ? 0 : alpha*beta*pow(x,beta-1));
    }

	double cumulative_density(double x) {
		return alpha*pow(x,beta);
	}

	double inverse_cumulative_density(double x) {
		 return pow(x/alpha,1/beta);

	}

	double density_derivative(double x) {
		return (x<=0 ? 0 : alpha*beta*(beta-1)*pow(x,beta-2));
	}

	double density_param_derivative(double x) {
		return (x==0 ? 0 : alpha*(1+beta*log(x))*pow(x,beta-1));
	}

	double cumulative_density_param_derivative(double x) {
		return (x==0 ? 0 : alpha*log(x)*pow(x,beta));
	}
 
};

// Thanks to Cecile Chauvel
class LogLinearFamilyModel : public FamilyModel {
  public:
    LogLinearFamilyModel(double alpha_, double beta_) {
      alpha=alpha_;beta=beta_;
    }
  
  ~LogLinearFamilyModel() {};
  
  double alpha, beta;
  
  NumericVector get_params() {
    NumericVector out(2);
    out[0]=alpha;out[1]=beta;
    return out;
  }
  
  void set_params(double alpha_, double beta_) {
    alpha=alpha_;beta=beta_;
  }
  
  double density(double x) {
    return (x<=0 ? 0 : alpha*exp(beta*x));
  }
  
  double cumulative_density(double x) {
    return alpha*(exp(beta*x)-1)/beta;
  }
  
  double inverse_cumulative_density(double x) {
    return log(1+x*beta/alpha)/beta;
    
  }
  
  double density_derivative(double x) {
    return (x<=0 ? 0 : alpha*beta*exp(beta*x));
  }
  
  double density_param_derivative(double x) {
    return  alpha*x*exp(beta*x) ;
  }
  
  double cumulative_density_param_derivative(double x) {
    return alpha*(x*exp(x*beta)/beta-(exp(beta*x)-1)/pow(beta,2));
  }
  
};

FamilyModel* newFamilyModel(List family);

#endif //RCPP_FAMILY_MODEL_H
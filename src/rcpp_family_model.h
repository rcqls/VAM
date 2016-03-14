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

	virtual double density(double x) = 0;  //TODO: hazard_rate

	virtual double cumulative_density(double x) = 0;

	virtual double inverse_cumulative_density(double x) = 0;

	virtual double density_derivative(double x) = 0;

  virtual double inverse_density(double x) = 0;

	virtual double density_param_derivative(double x) = 0;

	virtual double cumulative_density_param_derivative(double x) = 0;

  virtual double density_derivative_param_derivative(double x) = 0;//LD

  virtual double density_2derivative(double x) = 0;//LD

  virtual double density_param_2derivative(double x) = 0;//LD

  virtual double cumulative_density_param_2derivative(double x) = 0;//LD
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

  double inverse_density(double x) {
    return pow(x/alpha/beta,1/(beta-1));
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

  double density_derivative_param_derivative(double x) {//LD
    return (x==0 ? 0 : alpha*(2*beta-1+beta*(beta-1)*log(x))*pow(x,beta-2));//LD
  }//LD

  double density_2derivative(double x) {//LD
    return (x<=0 ? 0 : alpha*beta*(beta-1)*(beta-2)*pow(x,beta-3));//LD
  }//LD

  double density_param_2derivative(double x) {//LD
    return (x==0 ? 0 : alpha*(2+beta*log(x))*log(x)*pow(x,beta-1));//LD
  }//LD

  double cumulative_density_param_2derivative(double x) {//LD
    return (x==0 ? 0 : alpha*pow(log(x),2)*pow(x,beta));//LD
  }//LD
 
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

  double inverse_density(double x) {
    return log(x/alpha)/beta;
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
  
  double density_derivative_param_derivative(double x) {//LD
    return  alpha*exp(beta*x)*(1+beta*x) ;//LD
  }//LD

  double density_2derivative(double x) {//LD
    return (x<=0 ? 0 : alpha*pow(beta,2)*exp(beta*x));//LD
  }//LD

  double density_param_2derivative(double x) {//LD
    return  alpha*pow(x,2)*exp(beta*x) ;//LD
  }//LD
  
  double cumulative_density_param_2derivative(double x) {//LD
    return alpha*(pow(x,2)*exp(x*beta)/beta-2*x*exp(x*beta)/pow(beta,2)+2*(exp(beta*x)-1)/pow(beta,3));//LD
  }//LD


};

FamilyModel* newFamilyModel(List family);

#endif //RCPP_FAMILY_MODEL_H
#ifndef RCPP_FAMILY_MODEL_H
#define RCPP_FAMILY_MODEL_H
#include <Rcpp.h> 

using namespace Rcpp ;

class FamilyModel {
public:

	FamilyModel() {};

	virtual ~FamilyModel() {};

	virtual NumericVector get_params() = 0;

	virtual void set_params(NumericVector par) = 0;//LD3

  virtual int nb_params() = 0;//LD3

	virtual double hazardRate(double x) = 0;  //TODO: hazard_rate

	virtual double cumulative_hazardRate(double x) = 0;

	virtual double inverse_cumulative_hazardRate(double x) = 0;

	virtual double hazardRate_derivative(double x) = 0;

  virtual double inverse_hazardRate(double x) = 0;

	virtual NumericVector hazardRate_param_derivative(double x) = 0;//LD3//Size:nb_params-1

	virtual NumericVector cumulative_hazardRate_param_derivative(double x) = 0;//LD3:nb_params-1

  virtual NumericVector hazardRate_derivative_param_derivative(double x) = 0;//LD//LD3:nb_params-1

  virtual double hazardRate_2derivative(double x) = 0;//LD

  virtual NumericVector hazardRate_param_2derivative(double x) = 0;//LD//LD3:nb_params*(:nb_params-1)/2

  virtual NumericVector cumulative_hazardRate_param_2derivative(double x) = 0;//LD//LD3::nb_params*(:nb_params-1)/2
};

class WeibullFamilyModel : public FamilyModel {
public:
  WeibullFamilyModel(NumericVector par) {//LD3
  	alpha=par[0];beta=par[1];//LD3
  }

  ~WeibullFamilyModel() {};

  double alpha, beta;

  NumericVector get_params() {
  	NumericVector out(2);
  	out[0]=alpha;out[1]=beta;
  	return out;
  }

  void set_params(NumericVector par) {//LD3
  		alpha=par[0];beta=par[1];//LD3
  }

  int nb_params() {//LD3
    return 2;//LD3
  }//LD3

  double hazardRate(double x) {
  	return (x<=0 ? 0 : alpha*beta*pow(x,beta-1));
  }

  double inverse_hazardRate(double x) {
    return pow(x/alpha/beta,1/(beta-1));
  }

  double cumulative_hazardRate(double x) {
  	return alpha*pow(x,beta);
  }

  double inverse_cumulative_hazardRate(double x) {
  	 return pow(x/alpha,1/beta);

  }

  double hazardRate_derivative(double x) {
  	return (x<=0 ? 0 : alpha*beta*(beta-1)*pow(x,beta-2));
  }

  NumericVector hazardRate_param_derivative(double x) {//LD3
    NumericVector res(1);//LD3
  	res[0]=(x==0 ? 0 : alpha*(1+beta*log(x))*pow(x,beta-1));//LD3
    return res;//LD3
  }

  NumericVector cumulative_hazardRate_param_derivative(double x) {//LD3
  	NumericVector res(1);//LD3
    res[0]= (x==0 ? 0 : alpha*log(x)*pow(x,beta));//LD3
    return res;//LD3
  }

  NumericVector hazardRate_derivative_param_derivative(double x) {//LD//LD3
    NumericVector res(1);//LD3
    res[0]= (x==0 ? 0 : alpha*(2*beta-1+beta*(beta-1)*log(x))*pow(x,beta-2));//LD//LD3
    return res;//LD3
  }//LD

  double hazardRate_2derivative(double x) {//LD
    return (x<=0 ? 0 : alpha*beta*(beta-1)*(beta-2)*pow(x,beta-3));//LD
  }//LD

  NumericVector hazardRate_param_2derivative(double x) {//LD//LD3
    NumericVector res(1);//LD3
    res[0]= (x==0 ? 0 : alpha*(2+beta*log(x))*log(x)*pow(x,beta-1));//LD//LD3
    return res;//LD3
  }//LD

  NumericVector cumulative_hazardRate_param_2derivative(double x) {//LD//LD3
    NumericVector res(1);//LD3
    res[0]= (x==0 ? 0 : alpha*pow(log(x),2)*pow(x,beta));//LD//LD3
    return res;//LD3
  }//LD
 
};

// Thanks to Cecile Chauvel
class LogLinearFamilyModel : public FamilyModel {
  public:
    LogLinearFamilyModel(NumericVector par) {//LD3
      alpha=par[0];beta=par[1];//LD3
    }
  
  ~LogLinearFamilyModel() {};
  
  double alpha, beta;
  
  NumericVector get_params() {
    NumericVector out(2);
    out[0]=alpha;out[1]=beta;
    return out;
  }
  
  void set_params(NumericVector par) {//LD3
    alpha=par[0];beta=par[1];//LD3
  }

  int nb_params() {//LD3
    return 2;//LD3
  }//LD3
  
  double hazardRate(double x) {
    return (x<=0 ? 0 : alpha*exp(beta*x));
  }

  double inverse_hazardRate(double x) {
    return log(x/alpha)/beta;
  }
  
  double cumulative_hazardRate(double x) {
    return alpha*(exp(beta*x)-1)/beta;
  }
  
  double inverse_cumulative_hazardRate(double x) {
    return log(1+x*beta/alpha)/beta;
    
  }
  
  double hazardRate_derivative(double x) {
    return (x<=0 ? 0 : alpha*beta*exp(beta*x));
  }
  
  NumericVector hazardRate_param_derivative(double x) {//LD3
    NumericVector res(0);//LD3
    res[0]= alpha*x*exp(beta*x) ;//LD3
    return res;//LD3
  }
  
  NumericVector cumulative_hazardRate_param_derivative(double x) {//LD3
    NumericVector res(0);//LD3
    res[0]=  alpha*(x*exp(x*beta)/beta-(exp(beta*x)-1)/pow(beta,2));//LD3
    return res;//LD3
  }
  
  NumericVector hazardRate_derivative_param_derivative(double x) {//LD//LD3
    NumericVector res(0);//LD3
    res[0]= alpha*exp(beta*x)*(1+beta*x) ;//LD//LD3
    return res;//LD3
  }//LD

  double hazardRate_2derivative(double x) {//LD
    return (x<=0 ? 0 : alpha*pow(beta,2)*exp(beta*x));//LD
  }//LD

  NumericVector hazardRate_param_2derivative(double x) {//LD//LD3
    NumericVector res(0);//LD3
    res[0]=  alpha*pow(x,2)*exp(beta*x) ;//LD//LD3
    return res;//LD3
  }//LD
  
  NumericVector cumulative_hazardRate_param_2derivative(double x) {//LD//LD3
    NumericVector res(0);//LD3
    res[0]= alpha*(pow(x,2)*exp(x*beta)/beta-2*x*exp(x*beta)/pow(beta,2)+2*(exp(beta*x)-1)/pow(beta,3));//LD//LD3
    return res;//LD3
  }//LD


};

FamilyModel* newFamilyModel(List family);

#endif //RCPP_FAMILY_MODEL_H
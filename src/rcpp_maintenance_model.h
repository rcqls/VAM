#ifndef RCPP_VAM_MODEL_H
#define RCPP_VAM_MODEL_H
#include <Rcpp.h>
#include "rcpp_vam_model.h"

using namespace Rcpp ;

//Forward declarations
class MaintenanceModel;

//Effective declarations
class MaintenanceModelList {//List of ModelBase (heterogeneous terms) 
public:
    MaintenanceModelList(List models_,VamModel* model);

    ~MaintenanceModelList();

    MaintenanceModel* at(int i) {
    	return model_list[i];
    }

    int size() {
    	return model_list.size();
    }


protected:

    std::vector<MaintenanceModel*> model_list; //model list
     
};


class MaintenanceModel {
public:
    MaintenanceModel(VamModel* model_) {
    	model = model_;
    }

    virtual ~MaintenanceModel() {};

    virtual NumericVector get_params() = 0;

    virtual  void set_params(NumericVector par, int ind) = 0;//ind indicates the indice of vector par at which the parameters to set begin

    virtual int nb_params()=0;
 
    virtual void update(bool with_gradient,bool with_hessian) = 0;

    virtual double virtual_age(double time) = 0;

    virtual double* virtual_age_derivative(double x) = 0;

    virtual double* virtual_age_hessian(double x) = 0;

    virtual double virtual_age_inverse(double time) = 0;



    VamModel* model;

    void set_id(int id_) {
    	id=id_;
    }

    void set_id_params(int id_params_) {
        id_params=id_params_;
    }

    int id;
    int id_params;

};

class ARA1 : public MaintenanceModel { 

public:

    ARA1(NumericVector rho_,VamModel* model_) : MaintenanceModel(model_) {
    	rho = rho_[0];
    }

    NumericVector get_params() {
    	NumericVector out(1);
    	out[0]=rho;
    	return out;
    }

    void set_params(NumericVector par, int ind) {
    	rho=par[ind];
    }

    int nb_params(){
        return 1;
    }

    void update(bool with_gradient,bool with_hessian);

    double virtual_age(double time);

    double* virtual_age_derivative(double x);

    double* virtual_age_hessian(double x);

    double virtual_age_inverse(double x);

private:
    double rho;

};

class ARAInf : public MaintenanceModel { 

public:

    ARAInf(NumericVector rho_,VamModel* model_) : MaintenanceModel(model_) {
    	rho = rho_[0];
    }

    NumericVector get_params() {
    	NumericVector out(1);
    	out[0]=rho;
    	return out;
    }

    void set_params(NumericVector par, int ind) {
    	rho=par[ind];
    }

    int nb_params(){
        return 1;
    }

    void update(bool with_gradient,bool with_hessian);

    double virtual_age(double time);

    double* virtual_age_derivative(double x);

    double* virtual_age_hessian(double x);

    double virtual_age_inverse(double x);

private:

	double rho;

};

class AGAN : public MaintenanceModel {

public:

    AGAN(VamModel* model_) : MaintenanceModel(model_) {
    }

    NumericVector get_params() {
        NumericVector out(0);
        return out;
    }

    void set_params(NumericVector par, int ind) {
    }

    int nb_params(){
        return 0;
    }

    void update(bool with_gradient,bool with_hessian);

    double virtual_age(double time);

    double* virtual_age_derivative(double x);

    double* virtual_age_hessian(double x);

    double virtual_age_inverse(double x);

};

class ABAO : public MaintenanceModel { 

public:

    ABAO(VamModel* model_) : MaintenanceModel(model_) {
    }

    NumericVector get_params() {
        NumericVector out(0);
        return out;
    }

    void set_params(NumericVector par,int ind) {
    }

    int nb_params(){
        return 0;
    }

    void update(bool with_gradient,bool with_hessian);

    double virtual_age(double time);

    double* virtual_age_derivative(double x);

    double* virtual_age_hessian(double x);

    double virtual_age_inverse(double x);

};

MaintenanceModel* newMaintenanceModel(List maintenance,VamModel* model);

#endif //RCPP_VAM_MODEL_H
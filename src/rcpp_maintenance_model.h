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
        //Rprintf("at(i)=%d,%d,%d\n",i,size(),model_list.size());
    	return model_list[i];
    }

    int size() {
    	return size_;
    }


protected:

    std::vector<MaintenanceModel*> model_list; //model list
    int size_;
};


class MaintenanceModel {
public:
    MaintenanceModel(VamModel* model_) {
    	model = model_;
    }

    virtual ~MaintenanceModel() {};

    virtual NumericVector get_params() = 0;

    virtual  void set_params(double) = 0;
 
    virtual void update(bool with_gradient) = 0;

    virtual double virtual_age(double time) = 0;

    //virtual double* 
    virtual std::vector<double> virtual_age_derivative(double x) = 0;

    virtual double virtual_age_inverse(double time) = 0;

    VamModel* model;

    void set_id(int id_) {
    	id=id_;
    }

    int id;

};

class ARA1 : public MaintenanceModel { 

public:

    ARA1(double rho_,VamModel* model_) : MaintenanceModel(model_) {
    	rho = rho_;
    }

    virtual ~ARA1() {};

    NumericVector get_params() {
    	NumericVector out(1);
    	out[0]=rho;
    	return out;
    }

    void set_params(double par) {
    	rho=par;
    }

    void update(bool with_gradient); 

    double virtual_age(double time);

    //double* 
    std::vector<double> virtual_age_derivative(double x);

    double virtual_age_inverse(double x);

private:
    double rho;

};

class ARAInf : public MaintenanceModel { 

public:

    ARAInf(double rho_,VamModel* model_) : MaintenanceModel(model_) {
    	rho = rho_;
    }

    virtual ~ARAInf() {};

    NumericVector get_params() {
    	NumericVector out(1);
    	out[0]=rho;
    	return out;
    }

    void set_params(double par) {
    	rho=par;
    }

    void update(bool with_gradient);

    double virtual_age(double time);

    //double* 
    std::vector<double> virtual_age_derivative(double x);

    double virtual_age_inverse(double x);

private:

	double rho;

};

MaintenanceModel* newMaintenanceModel(List maintenance,VamModel* model);

#endif //RCPP_VAM_MODEL_H
#ifndef RCPP_VAM_CACHE_H
#define RCPP_VAM_CACHE_H
#include <Rcpp.h>
#include "rcpp_family_model.h"
#include "rcpp_maintenance_policy.h"

using namespace Rcpp ;

class MaintenanceModelList;

class VamModel {
public:

	VamModel(List model_) {
        init(model_);
    };

	VamModel(List model_,List data_) {
		init(model_);
		set_data(data_);
	};

	~VamModel() {
		//DEBUG: printf("VamModel: %p, %p, %p, %p, %p, %p, %p\n",dVright,dVleft,dS1,dS2,models,family,maintenance_policy);
		delete[] dVright;
		delete[] dVleft;
		delete[] dS1;
		delete[] dS2;
		delete models;
		delete family;
		delete maintenance_policy;
	};

	int k,nbPM,idMod,nb_system;

	List data;

	NumericVector time, type;

	double S1, S2, S3, indType;

	double Vleft, Vright, hVleft;

	double *dVleft, *dVright, *dS1, *dS2;

	MaintenanceModelList* models;

	FamilyModel* family;

	MaintenancePolicy* maintenance_policy;

	void initMLE() {
		int i;
		k=0;
		Vleft=0;
		Vright=0;
		indType=0;hVleft=0;
		dS1[0]=0;dS2[0]=0;
		for (i=0;i<nbPM+1;i++) {
			dVright[i]=0;
			dVleft[i]=0;
			dS1[i+1]=0;
			dS2[i+1]=0;
		}
		
	};

	FamilyModel* get_family() {
		return family;
	}

	List get() {
		List ret;
		ret["S1"]=NumericVector::create(S1);ret["S2"]=NumericVector::create(S2);ret["S3"]=NumericVector::create(S3);
		ret["Vright"]=NumericVector::create(Vright);ret["Vleft"]=NumericVector::create(Vleft);
		NumericVector dS1R(nbPM+2),dS2R(nbPM+2),dVrightR(nbPM+1),dVleftR(nbPM+1);
		ret["dS1"]=dS1R;ret["dS2"]=dS2R;
		ret["dVright"]=dVrightR;ret["dVleft"]=dVleftR;
		dS1R[0]=dS1[0];dS2R[0]=dS2[0];
		for (int i=0;i<nbPM+1;i++) {
			dVrightR[i]=dVright[i];
			dVleftR[i]=dVleft[i];
			dS1R[i+1]=dS1[i+1];
			dS2R[i+1]=dS2[i+1];
		}
		return ret;
	}

	void set_data(List data_);

	void select_data(int i);

	DataFrame get_selected_data(int i);

	NumericVector get_params();

    void set_params(NumericVector pars);

    void update_Vleft(bool with_gradient);

    List get_virtual_age_infos(double by);

private:
	void set_models(List models_);

    void set_family(List family_);

    void set_maintenance_policy(List maintenance_policy_);

	void init(List model_);

	void init_virtual_age_infos();

	DataFrame get_virtual_age_info(double from,double to,double by);

};

#endif
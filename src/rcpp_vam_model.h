#ifndef RCPP_VAM_CACHE_H
#define RCPP_VAM_CACHE_H
#include <Rcpp.h>
#include "rcpp_family_model.h"
#include "rcpp_maintenance_policy.h"

using namespace Rcpp ;

class MaintenanceModelList;
class MaintenancePolicy;

class VamModel {
public:

	VamModel(List model_) {
        init(model_);
    };

	VamModel(List model_,List data_) {
		init(model_);
		set_data(data_);
	};

	~VamModel();

	/***************************
	k: current position
 	nb_system: number of system
	****************************/
	int k,nb_system;

	int nbPM,idMod;
	int nb_paramsMaintenance,nb_paramsFamily;

	List data;

	//NumericVector time, type;
	std::vector<double> time;
	std::vector<int> type;

	double S1, S2, S3, indType;

	double Vleft, Vright, hVleft;

	double *dVleft, *dVright, *dS1, *dS2;
	double *d2Vleft, *d2Vright, *d2S1, *d2S2;

	MaintenanceModelList* models;

	FamilyModel* family;

	MaintenancePolicy* maintenance_policy;

	FamilyModel* get_family() {
	 	return family;
	}

	List get() {
		int j;
		int n_params=nb_paramsMaintenance+nb_paramsFamily-1;

		List ret;
		ret["S1"]=NumericVector::create(S1);ret["S2"]=NumericVector::create(S2);ret["S3"]=NumericVector::create(S3);
		ret["Vright"]=NumericVector::create(Vright);ret["Vleft"]=NumericVector::create(Vleft);
		NumericVector dS1R(n_params),dS2R(n_params);
		NumericMatrix d2S1R(n_params,n_params),d2S2R(n_params,n_params);
		ret["dS1"]=dS1R;ret["dS2"]=dS2R;
		ret["d2S1"]=d2S1R;ret["d2S2"]=d2S2R;
		for (int i=0;i<n_params;i++) {
			dS1R[i]=dS1[i];
			dS2R[i]=dS2[i];
			d2S1R(i,i)=d2S1[i*(i+1)/2+i];
			d2S2R(i,i)=d2S2[i*(i+1)/2+i];
			for (j=0;j<i;j++) {
				//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
				d2S1R(i,j)=d2S1[i*(i+1)/2+j];
				d2S2R(i,j)=d2S2[i*(i+1)/2+j];
				d2S1R(j,i)=d2S1R(i,j);
				d2S2R(j,i)=d2S2R(i,j);
			}
		}
		NumericVector dVrightR(nb_paramsMaintenance),dVleftR(nb_paramsMaintenance);
		NumericMatrix d2VrightR(nb_paramsMaintenance,nb_paramsMaintenance),d2VleftR(nb_paramsMaintenance,nb_paramsMaintenance);
		ret["dVright"]=dVrightR;ret["dVleft"]=dVleftR;
		ret["d2Vright"]=d2VrightR;ret["d2Vleft"]=d2VleftR;
		for (int i=0;i<nb_paramsMaintenance;i++) {
			dVrightR[i]=dVright[i];
			dVleftR[i]=dVleft[i];
			d2VrightR(i,i)=d2Vright[i*(i+1)/2+i];
			d2VleftR(i,i)=d2Vleft[i*(i+1)/2+i];
			for (j=0;j<i;j++) {
				//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
				d2VrightR(i,j)=d2Vright[i*(i+1)/2+j];
				d2VleftR(i,j)=d2Vleft[i*(i+1)/2+j];
				d2VrightR(j,i)=d2VrightR(i,j);
				d2VleftR(j,i)=d2VleftR(i,j);
			}
		}
		return ret;
	}

	void set_data(List data_);

	void select_data(int i);

	DataFrame get_selected_data(int i);

	NumericVector get_params();

    void set_params(NumericVector pars);

    void update_Vleft(bool with_gradient,bool with_hessian);

    List get_virtual_age_infos(double by);

		void init_computation_values();

private:
	void set_models(List models_);

    void set_family(List family_);

    void set_maintenance_policy(List maintenance_policy_);

	void init(List model_);

	void init_virtual_age_infos();

	DataFrame get_virtual_age_info(double from,double to,double by);

};

#endif

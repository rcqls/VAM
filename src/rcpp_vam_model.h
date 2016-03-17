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
	int nb_paramsMaintenance,nb_paramsFamily;//LD3

	List data;

	//NumericVector time, type;
	std::vector<double> time;
	std::vector<int> type;

	double S1, S2, S3, indType;

	double Vleft, Vright, hVleft;

	double *dVleft, *dVright, *dS1, *dS2;
	double *d2Vleft, *d2Vright, *d2S1, *d2S2;//LD

	MaintenanceModelList* models;

	FamilyModel* family;

	MaintenancePolicy* maintenance_policy;

	// void initMLE() {
	// 	int i;
	// 	int j;//LD
	// 	k=0;
	// 	Vleft=0;
	// 	Vright=0;
	// 	indType=0;hVleft=0;
	// 	if (nb_paramsCM+nb_paramsPM>0){//LD3
	// 	  dS1[0]=0;dS2[0]=0;
	// 	  d2S1[0]=0;//LD
	//       d2S2[0]=0;//LD
	// 	  for (i=0;i<nb_paramsPM+1;i++) {
	// 		dVright[i]=0;
	// 		dVleft[i]=0;
	// 		dS1[i+1]=0;
	// 		dS2[i+1]=0;
	// 		d2S1[(i+1)*(i+2)/2]=0;//LD
	// 		d2S2[(i+1)*(i+2)/2]=0;//LD
	// 	  	for (j=0;j<=i;j++) {//LD
	// 		  	d2Vright[i*(i+1)/2+j]=0;//LD
	// 		  	d2Vleft[i*(i+1)/2+j]=0;//LD
	// 		  	d2S1[(i+1)*(i+2)/2+j+1]=0;//LD
	// 		  	d2S2[(i+1)*(i+2)/2+j+1]=0;//LD
	// 	  	}//LD
	// 	  }
	// 	}//LD3

	// };

	FamilyModel* get_family() {
	 	return family;
	}

	List get() {
		int j;//LD
		int n_params=nb_paramsMaintenance+nb_paramsFamily-1;//LD3

		List ret;
		ret["S1"]=NumericVector::create(S1);ret["S2"]=NumericVector::create(S2);ret["S3"]=NumericVector::create(S3);
		ret["Vright"]=NumericVector::create(Vright);ret["Vleft"]=NumericVector::create(Vleft);
		NumericVector dS1R(n_params),dS2R(n_params);//LD3
		NumericMatrix d2S1R(n_params,n_params),d2S2R(n_params,n_params);//LD3
		ret["dS1"]=dS1R;ret["dS2"]=dS2R;
		ret["d2S1"]=d2S1R;ret["d2S2"]=d2S2R;//LD
		for (int i=0;i<n_params;i++) {//LD3
			dS1R[i]=dS1[i];//LD3
			dS2R[i]=dS2[i];//LD3
			d2S1R(i,i)=d2S1[i*(i+1)/2+i];//LD3
			d2S2R(i,i)=d2S2[i*(i+1)/2+i];//LD3
			for (j=0;j<i;j++) {//LD3
				//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
				d2S1R(i,j)=d2S1[i*(i+1)/2+j];//LD3
				d2S2R(i,j)=d2S2[i*(i+1)/2+j];//LD3
				d2S1R(j,i)=d2S1R(i,j);//LD3
				d2S2R(j,i)=d2S2R(i,j);//LD3
			}//LD3
		}//LD3
		NumericVector dVrightR(nb_paramsMaintenance),dVleftR(nb_paramsMaintenance);//LD3
		NumericMatrix d2VrightR(nb_paramsMaintenance,nb_paramsMaintenance),d2VleftR(nb_paramsMaintenance,nb_paramsMaintenance);//LD3
		ret["dVright"]=dVrightR;ret["dVleft"]=dVleftR;
		ret["d2Vright"]=d2VrightR;ret["d2Vleft"]=d2VleftR;//LD
		for (int i=0;i<nb_paramsMaintenance;i++) {//LD3
			dVrightR[i]=dVright[i];
			dVleftR[i]=dVleft[i];
			d2VrightR(i,i)=d2Vright[i*(i+1)/2+i];//LD
			d2VleftR(i,i)=d2Vleft[i*(i+1)/2+i];//LD
			for (j=0;j<i;j++) {//LD
				//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
				d2VrightR(i,j)=d2Vright[i*(i+1)/2+j];//LD
				d2VleftR(i,j)=d2Vleft[i*(i+1)/2+j];//LD
				d2VrightR(j,i)=d2VrightR(i,j);//LD
				d2VleftR(j,i)=d2VleftR(i,j);//LD
			}//LD
		}
		return ret;
	}

	void set_data(List data_);

	void select_data(int i);

	DataFrame get_selected_data(int i);

	NumericVector get_params();

    void set_params(NumericVector pars);

    void update_Vleft(bool with_gradient,bool with_hessian);//LD

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

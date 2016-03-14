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

	void initMLE() {
		int i;
		int j;//LD
		k=0;
		Vleft=0;
		Vright=0;
		indType=0;hVleft=0;
		dS1[0]=0;dS2[0]=0;
		d2S1[0]=0;//LD
	    d2S2[0]=0;//LD
		for (i=0;i<nbPM+1;i++) {
			dVright[i]=0;
			dVleft[i]=0;
			dS1[i+1]=0;
			dS2[i+1]=0;
			d2S1[(i+1)*(i+2)/2]=0;//LD
			d2S2[(i+1)*(i+2)/2]=0;//LD
		  	for (j=0;j<=i;j++) {//LD
			  	d2Vright[i*(i+1)/2+j]=0;//LD
			  	d2Vleft[i*(i+1)/2+j]=0;//LD
			  	d2S1[(i+1)*(i+2)/2+j+1]=0;//LD
			  	d2S2[(i+1)*(i+2)/2+j+1]=0;//LD
		  	}//LD
		}

	};

	FamilyModel* get_family() {
		return family;
	}

	List get() {
		int j;//LD
		List ret;
		ret["S1"]=NumericVector::create(S1);ret["S2"]=NumericVector::create(S2);ret["S3"]=NumericVector::create(S3);
		ret["Vright"]=NumericVector::create(Vright);ret["Vleft"]=NumericVector::create(Vleft);
		NumericVector dS1R(nbPM+2),dS2R(nbPM+2),dVrightR(nbPM+1),dVleftR(nbPM+1);
		NumericMatrix d2S1R(nbPM+2,nbPM+2),d2S2R(nbPM+2,nbPM+2),d2VrightR(nbPM+1,nbPM+1),d2VleftR(nbPM+1,nbPM+1);//LD
		ret["dS1"]=dS1R;ret["dS2"]=dS2R;
		ret["dVright"]=dVrightR;ret["dVleft"]=dVleftR;
		dS1R[0]=dS1[0];dS2R[0]=dS2[0];
		d2S1R(0,0)=d2S1[0];d2S2R(0,0)=d2S2[0];//LD
		for (int i=0;i<nbPM+1;i++) {
			dVrightR[i]=dVright[i];
			dVleftR[i]=dVleft[i];
			dS1R[i+1]=dS1[i+1];
			dS2R[i+1]=dS2[i+1];
			d2S1R(0,i+1)=d2S1[(i+1)*(i+2)/2];//LD
			d2S2R(0,i+1)=d2S2[(i+1)*(i+2)/2];//LD
			d2S1R(i+1,0)=d2S1R(0,i+1);//LD
			d2S2R(i+1,0)=d2S2R(0,i+1);//LD
			d2VrightR(i,i)=d2Vright[i*(i+1)/2+i];//LD
			d2VleftR(i,i)=d2Vleft[i*(i+1)/2+i];//LD
			d2S1R(i+1,i+1)=d2S1[(i+1)*(i+2)/2+i+1];//LD
			d2S2R(i+1,i+1)=d2S2[(i+1)*(i+2)/2+i+1];//LD
			for (j=0;j<i;j++) {//LD
				//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
				d2VrightR(i,j)=d2Vright[i*(i+1)/2+j];//LD
				d2VleftR(i,j)=d2Vleft[i*(i+1)/2+j];//LD
				d2S1R(i+1,j+1)=d2S1[(i+1)*(i+2)/2+j+1];//LD
				d2S2R(i+1,j+1)=d2S2[(i+1)*(i+2)/2+j+1];//LD
				d2VrightR(j,i)=d2VrightR(i,j);//LD
				d2VleftR(j,i)=d2VleftR(i,j);//LD
				d2S1R(j+1,i+1)=d2S1R(i+1,j+1);//LD
				d2S2R(j+1,i+1)=d2S2R(i+1,j+1);//LD
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

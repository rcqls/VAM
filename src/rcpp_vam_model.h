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
	int nb_paramsMaintenance,nb_paramsFamily,nb_paramsCov;

	List data;

	//Additional Covariates stuff
	DataFrame data_cov;
	NumericVector params_cov;
	double sum_cov; //to save the computation 

	//NumericVector time, type;
	std::vector<double> time;
	std::vector<int> type;

	double S1, S2, S0, S3, indType;

	double Vleft, Vright, hVleft;

	double *dVleft, *dVright, *dS1, *dS2, *dS3;
	double *d2Vleft, *d2Vright, *d2S1, *d2S2, *d2S3;

	double A;
	double *dA;
	double *d2A;

	int mu;
	double *VR_prec;
	double *dVR_prec,*d2VR_prec;



	MaintenanceModelList* models;

	FamilyModel* family;

	MaintenancePolicy* maintenance_policy;

	FamilyModel* get_family() {
	 	return family;
	}

	//Ununsed and can create bugs
	//NumericVector get_dB(int k) {
	// 	NumericVector dBkR(nb_paramsMaintenance);
	// 	for (int i=0;i<nb_paramsMaintenance;i++){
	// 		dBkR[i]=dB[k*nb_paramsMaintenance+i];
	// 	}
	// 	return dBkR;
	// }

	// NumericMatrix get_d2B(int k){
	// 	int j;
	// 	int di=k*nb_paramsMaintenance*(nb_paramsMaintenance+1)/2;
	// 	NumericMatrix d2BkR(nb_paramsMaintenance,nb_paramsMaintenance);
	// 	for (int i=0;i<nb_paramsMaintenance;i++) {
	// 		d2BkR(i,i)=d2B[di+i*(i+1)/2+i];
	// 		for (j=0;j<i;j++) {
	// 			//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
	// 			d2BkR(i,j)=d2B[di+i*(i+1)/2+j];
	// 			d2BkR(j,i)=d2BkR(i,j);
	// 		}
	// 	}
	// 	return d2BkR;
	// }

	// List get() {
	// 	int j;
	// 	int n_params=nb_paramsMaintenance+nb_paramsFamily-1;

	// 	List ret;
	// 	ret["S1"]=NumericVector::create(S1);ret["S2"]=NumericVector::create(S2);ret["S0"]=NumericVector::create(S0);ret["S3"]=NumericVector::create(S3);
	// 	ret["Vright"]=NumericVector::create(Vright);ret["Vleft"]=NumericVector::create(Vleft);
	// 	//printf("S0=%f, S1=%f, S2=%f, S3=%f, Vright=%f, Vleft=%f\n",S0,S1,S2,S3,Vright,Vleft);
	// 	NumericVector dS1R(n_params),dS2R(n_params), dS3R(nb_paramsMaintenance);
	// 	NumericMatrix d2S1R(n_params,n_params),d2S2R(n_params,n_params), d2S3R(nb_paramsMaintenance,nb_paramsMaintenance);
	// 	ret["dS1"]=dS1R;ret["dS2"]=dS2R;ret["dS3"]=dS3R;
	// 	ret["d2S1"]=d2S1R;ret["d2S2"]=d2S2R;ret["d2S3"]=d2S3R;

	// 	for (int i=0;i<n_params;i++) {
	// 		dS1R[i]=dS1[i];
	// 		dS2R[i]=dS2[i];
	// 		//printf("dS1[%d]=%f, d2S2=%f\n",i,dS1[i],dS2[i]);
	// 		d2S1R(i,i)=d2S1[i*(i+1)/2+i];
	// 		d2S2R(i,i)=d2S2[i*(i+1)/2+i];
	// 		for (j=0;j<i;j++) {
	// 			//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
	// 			d2S1R(i,j)=d2S1[i*(i+1)/2+j];
	// 			d2S2R(i,j)=d2S2[i*(i+1)/2+j];
	// 			//printf("d2S1[%d,%d][%d]=%f, d2S2=%f\n",i,j,i*(i+1)/2+j,d2S1[i*(i+1)/2+j],d2S2[i*(i+1)/2+j]);
	// 			d2S1R(j,i)=d2S1R(i,j);
	// 			d2S2R(j,i)=d2S2R(i,j);
	// 			//printf("d2S1[%d,%d][%d]=%f, d2S2=%f\n",j,i,i*(i+1)/2+j,d2S1R(i,j),d2S2R(i,j));
	// 		}
	// 	}
	// 	NumericVector dVrightR(nb_paramsMaintenance),dVleftR(nb_paramsMaintenance);
	// 	NumericMatrix d2VrightR(nb_paramsMaintenance,nb_paramsMaintenance),d2VleftR(nb_paramsMaintenance,nb_paramsMaintenance);
	// 	NumericVector dAR(nb_paramsMaintenance);
	// 	NumericMatrix d2AR(nb_paramsMaintenance,nb_paramsMaintenance);
	// 	NumericVector dCR(nb_paramsMaintenance);
	// 	NumericMatrix d2CR(nb_paramsMaintenance,nb_paramsMaintenance);

	// 	ret["dVright"]=dVrightR;ret["dVleft"]=dVleftR;
	// 	ret["d2Vright"]=d2VrightR;ret["d2Vleft"]=d2VleftR;
	// 	ret["A"]=NumericVector::create(A);
	// 	ret["dA"]=dAR;
	// 	ret["d2A"]=d2AR;
	// 	ret["C"]=NumericVector::create(C);
	// 	ret["dC"]=dCR;
	// 	ret["d2C"]=d2CR;
	// 	for (int i=0;i<nb_paramsMaintenance;i++) {
	// 		dVrightR[i]=dVright[i];
	// 		dVleftR[i]=dVleft[i];
	// 		dS3R[i]=dS3[i];
	// 		dAR[i]=dA[i];
	// 		dCR[i]=dC[i];
	// 		d2VrightR(i,i)=d2Vright[i*(i+1)/2+i];
	// 		d2VleftR(i,i)=d2Vleft[i*(i+1)/2+i];
	// 		d2S3R(i,i)=d2S3[i*(i+1)/2+i];
	// 		d2AR(i,i)=d2A[i*(i+1)/2+i];
	// 		d2CR(i,i)=d2A[i*(i+1)/2+i];
	// 		for (j=0;j<i;j++) {
	// 			//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
	// 			d2VrightR(i,j)=d2Vright[i*(i+1)/2+j];
	// 			d2VleftR(i,j)=d2Vleft[i*(i+1)/2+j];
	// 			d2S3R(i,j)=d2S3[i*(i+1)/2+j];
	// 			d2AR(i,j)=d2A[i*(i+1)/2+j];
	// 			d2CR(i,j)=d2C[i*(i+1)/2+j];
	// 			d2VrightR(j,i)=d2VrightR(i,j);
	// 			d2VleftR(j,i)=d2VleftR(i,j);
	// 			d2S3R(j,i)=d2S3R(i,j);
	// 			d2AR(j,i)=d2AR(i,j);
	// 			d2CR(j,i)=d2CR(i,j);
	// 		}
	// 	}
	// 	List BR(max_mem);
	// 	List dBR(max_mem);
	// 	List d2BR(max_mem);
	// 	ret["B"]=BR;
	// 	ret["dB"]=dBR;
	// 	ret["d2B"]=d2BR;
	// 	for (j=0;j<max_mem;j++){
	// 		BR[j]=B[j];
	// 		dBR[j]=get_dB(j);
	// 		d2BR[j]=get_d2B(j);
	// 	}

	// 	return ret;
	// }

	void set_data(List data_);

	void select_data(int i);

	DataFrame get_selected_data(int i);

	NumericVector get_params();

    void set_params(NumericVector pars);

    double virtual_age(double x) ;

    double virtual_age_inverse(double x);

    void update_Vleft(bool with_gradient,bool with_hessian);

    List get_virtual_age_infos(double by, double from, double to);

	void init_computation_values();
	
	//Covariates related
	double compute_covariates(int i); //output maybe useful inside R

private:
	void set_models(List models_);

    void set_family(List family_);

    void set_maintenance_policy(List maintenance_policy_);

	void init(List model_);

	void init_virtual_age_infos();

	DataFrame get_virtual_age_info(double from,double to,double by);

	//Covariates related
	void set_covariates(List covariates_);

};

#endif

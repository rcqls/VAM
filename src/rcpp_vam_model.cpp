#include "rcpp_vam_model.h"
#include "rcpp_family_model.h"
#include "rcpp_maintenance_model.h"
#include "rcpp_maintenance_policy.h"

using namespace Rcpp ;

VamModel::~VamModel() {
	//DEBUG: printf("VamModel: %p, %p, %p, %p, %p, %p, %p\n",dVright,dVleft,dS1,dS2,models,family,maintenance_policy);
	delete[] dS1;
	delete[] dS2;
	delete[] dS3;
	delete[] d2S1;
	delete[] d2S2;
	delete[] d2S3;
	delete[] dVright;
	delete[] dVleft;
	delete[] d2Vright;
	delete[] d2Vleft;
	delete[] dA;
	delete[] d2A;
	if(mu>0){
		delete[] VR_prec;
		delete[] dVR_prec;
		delete[] d2VR_prec;
	}
	delete models;
	delete family;
	delete maintenance_policy;
};

NumericVector VamModel::get_params() {
	int i;
	int j;
	int k;
	NumericVector pars(nb_paramsMaintenance+nb_paramsFamily+nb_paramsCov);
	NumericVector fam=family->get_params();
	for(i=0;i<nb_paramsFamily;i++){
		pars[i]=fam[i];
	}
	j=nb_paramsFamily-1;
	for(i=0;i<nbPM + 1;i++) {
		MaintenanceModel* vam=models->at(i);
		NumericVector res=vam->get_params();
		for(k=0;k<vam->nb_params();k++){
			j++;
			pars[j]=res[k];
		}
	}
	//Covariates related
	for(k=0;k<nb_paramsCov;k++) {
		j++;
		pars[j]=params_cov[k];
	}
	return pars;
}

void VamModel::set_params(NumericVector pars) {
	int i;
	int j;
	double toto;
	if (pars.size()!=nb_paramsFamily+nb_paramsMaintenance+nb_paramsCov){
		if (pars.size()>nb_paramsFamily+nb_paramsMaintenance+nb_paramsCov){
			toto=time[1];
			printf("The length of the parameter vector is too big, some values are not considered !%f \n",toto);
		} else {
			printf("The length of the parameter vector is too small, the missing values are fixed to 0.5 !\n");
		}
		NumericVector pars2(nb_paramsMaintenance+nb_paramsFamily+nb_paramsCov);
		for(i=0;i<std::min(pars.size(),pars2.size());i++){
			pars2[i]=pars[i];
		}
		for(i=std::min(pars.size(),pars2.size());i<pars2.size();i++){
			pars2[i]=0.5;
		}
		if(nb_paramsFamily>0){
			family->set_params(pars2);
		}
		j=nb_paramsFamily;
		for(i=0;i<nbPM + 1;i++) {
			MaintenanceModel* vam=models->at(i);
			vam->set_params(pars2,j);
			j=j+vam->nb_params();
		}
		for(i=0;i<nb_paramsCov;i++) {
			params_cov[i]=pars2[j];
			j++;
		}
	} else {
		if(nb_paramsFamily>0){
			family->set_params(pars);
		}
		j=nb_paramsFamily;
		for(i=0;i<nbPM + 1;i++) {
			MaintenanceModel* vam=models->at(i);
			vam->set_params(pars,j);
			j=j+vam->nb_params();
		}
		for(i=0;i<nb_paramsCov;i++) {
			params_cov[i]=pars[j];
			j++;
		}
	}
}

double VamModel::virtual_age(double x) {
    return Vright + (x  - time[k])*A;
}

double VamModel::virtual_age_inverse(double x) {
    return (x - Vright)/A + time[k] ;
}

void VamModel::update_Vleft(bool with_gradient,bool with_hessian) {
	int i;
	int j;
	/*if(model->k < 10) printf("Vleft:%lf\n", model->Vleft);*/
	Vleft = virtual_age(time[k+1]);
	//printf("Vleft:%lf\n", model->Vleft);
	if(with_hessian) {
		for(i=0;i<nb_paramsMaintenance;i++) {
            dVleft[i]=dVright[i] + (time[k+1]  - time[k])*dA[i];
            for (j=0;j<=i;j++)
                d2Vleft[i*(i+1)/2+j]= d2Vright[i*(i+1)/2+j] + (time[k+1]  - time[k])*d2A[i*(i+1)/2+j];
        }
	}
	else if(with_gradient) {
		for(i=0;i<nb_paramsMaintenance;i++)
            dVleft[i]= dVright[i] + (time[k+1]  - time[k])*dA[i];
	}

}

void VamModel::set_data(List data_) {
	data=data_;
	nb_system=data.size();
	//printf("Number of systems: %d\n",nb_system);
	select_data(0);//default when only one system no need to
}

void VamModel::select_data(int i) {
	//In particular, if no data the following is skipped!
	if(data.size() > i) {
		List data2=data[i];
		//OLD CODE before gcc6: time = data2[0]; type = data2[1];//0 stand for Time and 1 for Type
		time = as<std::vector<double> >(data2[0]); type = as<std::vector<int> >(data2[1]); //Thanks to Lea, seems to be related to introduction of gcc Version 6 (see also https://github.com/apache/incubator-mxnet/issues/2185)
	}
}

DataFrame VamModel::get_selected_data(int i) {
	select_data(i);//Skipped if data is unset (see above)
	return DataFrame::create(_["Time"]=time,_["Type"]=type);
};


void VamModel::set_models(List models_) {
    models=new MaintenanceModelList(models_,this);
}

void VamModel::set_family(List family_) {
	family=newFamilyModel(family_);
}

void VamModel::set_maintenance_policy(List maintenance_policy_) {
	maintenance_policy=newMaintenancePolicy(maintenance_policy_);
	//if(maintenance_policy==NULL) printf("maintenance_policy is NULL\n");
};

void VamModel::init_computation_values() {
	int i;
	S1=0;S2=0;S0=0;S3=0;
	Vleft=0;Vright=0;
	hVleft=0;
	A=1;
	for(i=0;i<nbPM + 1;i++) models->at(i)->init();
}

void VamModel::init(List model_) {
	mu=model_["max_memory"];mu--;
	List models_=model_["models"];
	List family_=model_["family"];
	List maintenance_policy_=model_["pm.policy"];
  	set_models(models_);
	nbPM=models->size()-1;
	nb_paramsMaintenance=0;
	for(int i=0;i<nbPM + 1;i++) {
		nb_paramsMaintenance=nb_paramsMaintenance+models->at(i)->nb_params();
	}

	set_family(family_);
	nb_paramsFamily=family->nb_params();
	set_maintenance_policy(maintenance_policy_);

	set_covariates(model_);

	// S1=0;S2=0;S0=0;
	// Vleft=0;Vright=0;
	// hVleft=0;
	init_computation_values();
	dS1=new double[nb_paramsMaintenance+nb_paramsFamily-1];
	dS2=new double[nb_paramsMaintenance+nb_paramsFamily-1];
	dS3=new double[nb_paramsMaintenance];
	d2S1=new double[(nb_paramsMaintenance+nb_paramsFamily-1)*(nb_paramsMaintenance+nb_paramsFamily)/2];//inferior diagonal part of the hessian matrice by lines
	d2S2=new double[(nb_paramsMaintenance+nb_paramsFamily-1)*(nb_paramsMaintenance+nb_paramsFamily)/2];//inferior diagonal part of the hessian matrice by lines
	d2S3=new double[(nb_paramsMaintenance)*(nb_paramsMaintenance+1)/2];//inferior diagonal part of the hessian matrice by lines
	dVright=new double[nb_paramsMaintenance];
	dVleft=new double[nb_paramsMaintenance];
	d2Vright=new double[(nb_paramsMaintenance)*(nb_paramsMaintenance+1)/2];//inferior diagonal part of the hessian matrice by lines
	d2Vleft=new double[(nb_paramsMaintenance)*(nb_paramsMaintenance+1)/2];//inferior diagonal part of the hessian matrice by lines
	dA=new double[nb_paramsMaintenance];
	d2A=new double[(nb_paramsMaintenance)*(nb_paramsMaintenance+1)/2];//inferior diagonal part of the hessian matrice by lines
	if(mu>0){
		VR_prec=new double[mu];
		dVR_prec=new double[mu*nb_paramsMaintenance];//dVR_prec[i*nb_paramsMaintenance+j] for 0<=i<mu and 0<=j<nb_paramsMaintenance, corresponds to the j th partial derivative corresponding to the i th last inter-maintenance time effect
		d2VR_prec=new double[mu*(nb_paramsMaintenance)*(nb_paramsMaintenance+1)/2];//d2VR_prec[i*(nb_paramsMaintenance)*(nb_paramsMaintenance+1)/2+j] for 0<=i<mu and 0<=j<(nb_paramsMaintenance)*(nb_paramsMaintenance+1)/2, inferior diagonal part of the hessian matrice by lines
	}
};

void VamModel::init_virtual_age_infos() {
		int i;
    	k=0;
    	idMod=0; //id of current model
    	S1 = 0;
    	Vright=0;
    	A=1;
    	for(i=0;i<nbPM + 1;i++) models->at(i)->init();
};

DataFrame VamModel::get_virtual_age_info(double from,double to, double by) {
	double s=ceil((to-from)/by);
	int n=static_cast<int>(s);
//printf("ici=%d,%lf (%lf,%lf,%lf)\n",n,s,to,from,by);
	std::vector<double> t(n+1);
	std::vector<double> v(n+1);
	std::vector<double> h(n+1); //i as intensity
	std::vector<double> H(n+1); //I for cumulative intensity
	std::vector<double> F(n+1); //F for conditional cumulative distribution function
	std::vector<double> S(n+1); //S for conditional survival function
	std::vector<double> f(n+1); //S for conditional survival function


	t[0]=from;t[n]=to;
	v[0]=virtual_age(from);v[n]=virtual_age(to);
	h[0]=A*family->hazardRate(v[0]);h[n]=A*family->hazardRate(v[n]);
	H[0]=S1;H[n]=S1+family->cumulative_hazardRate(v[n])-family->cumulative_hazardRate(v[0]);
	F[0]=0;F[n]=1-exp(-(family->cumulative_hazardRate(v[n])-family->cumulative_hazardRate(v[0])));
	S[0]=1;S[n]=exp(-(family->cumulative_hazardRate(v[n])-family->cumulative_hazardRate(v[0])));
	f[0]=A*family->hazardRate(v[0]);f[n]=A*family->hazardRate(v[n])*exp(-(family->cumulative_hazardRate(v[n])-family->cumulative_hazardRate(v[0])));
	double by_t=(t[n]-t[0])/s;
	double by_v=(v[n]-v[0])/s;

	for(int i=1;i<n;i++) {
		t[i]=t[i-1]+by_t;//printf("t[%d]=%lf\n",i,t[i]);
		v[i]=v[i-1]+by_v;
		h[i]=A*family->hazardRate(v[i]);
		H[i]=S1+family->cumulative_hazardRate(v[i])-family->cumulative_hazardRate(v[0]);
		F[i]=1-exp(-(family->cumulative_hazardRate(v[i])-family->cumulative_hazardRate(v[0])));
		S[i]=exp(-(family->cumulative_hazardRate(v[i])-family->cumulative_hazardRate(v[0])));
		f[i]=A*family->hazardRate(v[i])*exp(-(family->cumulative_hazardRate(v[i])-family->cumulative_hazardRate(v[0])));
	}

	return DataFrame::create(
		_["t"]=NumericVector(t.begin(),t.end()),
		_["v"]=NumericVector(v.begin(),v.end()),
		_["i"]=NumericVector(h.begin(),h.end()),
		_["I"]=NumericVector(H.begin(),H.end()),
		_["F"]=NumericVector(F.begin(),F.end()),
		_["S"]=NumericVector(S.begin(),S.end()),
		_["f"]=NumericVector(f.begin(),f.end())
	);
};

List VamModel::get_virtual_age_infos(double by,double from, double to) {

	// Only one system first!
	init_virtual_age_infos();
	int n=time.size() - 1;
	List res(n);
	while(k < n) {
		//printf("k=%d/n=%d,(%lf,%lf)\n",k,n,time[k],time[k+1]);
		update_Vleft(false,false);
		if(from > time[k] || time[k+1] > to ) res[k] = R_NilValue;
		else res[k]=get_virtual_age_info(time[k],time[k+1],by);
		S1 += family->cumulative_hazardRate(Vleft) - family->cumulative_hazardRate(Vright);
		//gradient_update_for_current_system();
		int type2=type[k + 1];
		if(type2 < 0) type2=0;
		models->at(type2)->update(false,false);
	}
	return res;
};

//Covariates related
void VamModel::set_covariates(List model) {
	sum_cov=0.0;
	if(model["covariates"]==R_NilValue) {
		nb_paramsCov=0;
	} else {
		List covariates_=model["covariates"];
		data_cov=covariates_["data"];
		params_cov=covariates_["params"];
		nb_paramsCov=params_cov.size();
	}
}

double VamModel::compute_covariates(int i) {
	double sum_cov=0.0;
	for(int j=0;j<nb_paramsCov + 1;j++) {
		NumericVector var=data_cov[j];
		sum_cov += params_cov[j] * var[i-1]; //i-1 because R start from 1
	}
	return sum_cov;
}

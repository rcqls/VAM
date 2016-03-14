#include "rcpp_vam_model.h"
#include "rcpp_family_model.h"
#include "rcpp_maintenance_model.h"
#include "rcpp_maintenance_policy.h"

using namespace Rcpp ;

VamModel::~VamModel() {
	//DEBUG: printf("VamModel: %p, %p, %p, %p, %p, %p, %p\n",dVright,dVleft,dS1,dS2,models,family,maintenance_policy);
	delete[] dVright;
	delete[] dVleft;
	delete[] dS1;
	delete[] dS2;
	delete[] d2Vright;//LD
	delete[] d2Vleft;//LD
	delete[] d2S1;//LD
	delete[] d2S2;//LD
	delete models;
	delete family;
	delete maintenance_policy;
};

NumericVector VamModel::get_params() {
	NumericVector pars(nbPM+3);
	NumericVector fam=family->get_params();
	pars[0]=fam[0];pars[1]=fam[1];
	for(int i=0;i<nbPM + 1;i++) {
		MaintenanceModel* vam=models->at(i);
		NumericVector res=vam->get_params();
		pars[2+i]=vam->get_params()[0];
	}
	return pars;
}

void VamModel::set_params(NumericVector pars) {
	family->set_params(pars[0],pars[1]);
	for(int i=0;i<nbPM + 1;i++) {
		MaintenanceModel* vam=models->at(i);
		vam->set_params(pars[2+i]);
	}
}

void VamModel::update_Vleft(bool with_gradient,bool with_hessian) {//LD
	int j;//LD
	int i;//LD
	/*if(model->k < 10) printf("Vleft:%lf\n", model->Vleft);*/
	Vleft =(models->at(idMod))->virtual_age(time[k+1]);
	//printf("Vleft:%lf\n", model->Vleft);
	if(with_hessian) {//LD
		double* tmp=(models->at(idMod))->virtual_age_derivative(time[k+1]);//LD
		double* dtmp=(models->at(idMod))->virtual_age_hessian(time[k+1]);//LD
		for(i=0;i<nbPM+1;i++) {//LD
			dVleft[i]=tmp[i];//LD
			for (j=0;j<=i;j++) d2Vleft[i*(i+1)/2+j]=dtmp[i*(i+1)/2+j];//LD
		}//LD
	}//LD
	else if(with_gradient) {//LD
		double* tmp=(models->at(idMod))->virtual_age_derivative(time[k+1]);
		for(i=0;i<nbPM+1;i++) dVleft[i]=tmp[i];//LD:sortie la déclaration int i de la boucle
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
		time = data2["Time"]; type = data2["Type"];
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
	S1=0;S2=0;S3=0;
	Vleft=0;Vright=0;
	hVleft=0;
}

void VamModel::init(List model_) {
	List models_=model_["models"];
	List family_=model_["family"];
	List maintenance_policy_=model_["pm.policy"];
  set_models(models_);
	nbPM=models->size()-1;

	set_family(family_);
	set_maintenance_policy(maintenance_policy_);

	// S1=0;S2=0;S3=0;
	// Vleft=0;Vright=0;
	// hVleft=0;
	init_computation_values();

	dVright=new double[nbPM+1];
	dVleft=new double[nbPM+1];
	dS1=new double[nbPM+2];
	dS2=new double[nbPM+2];
	d2Vright=new double[(nbPM+1)*(nbPM+2)/2];//inferior diagonal part of the hessian matrice by lines//LD
	d2Vleft=new double[(nbPM+1)*(nbPM+2)/2];//inferior diagonal part of the hessian matrice by lines//LD
	d2S1=new double[(nbPM+2)*(nbPM+3)/2];//inferior diagonal part of the hessian matrice by lines//LD
	d2S2=new double[(nbPM+2)*(nbPM+3)/2];//inferior diagonal part of the hessian matrice by lines//LD
	//DEBUG: printf("dVright:%p,dVleft:%p\n",dVright,dVleft);
};

void VamModel::init_virtual_age_infos() {
    	k=0;
    	idMod=0; //id of current model
    	S1 = 0;
    	Vright=0;
};

DataFrame VamModel::get_virtual_age_info(double from,double to, double by) {
	double s=ceil((to-from)/by);
	int n=static_cast<int>(s);
//printf("ici=%d,%lf (%lf,%lf,%lf)\n",n,s,to,from,by);
	std::vector<double> t(n+1);
	std::vector<double> v(n+1);
	std::vector<double> h(n+1); //i as intensity
	std::vector<double> H(n+1); //I for cumulative intensity

	t[0]=from;t[n]=to;
	v[0]=Vright;v[n]=Vleft;
	h[0]=family->density(v[0]);h[n]=family->density(v[n]);
	H[0]=S1;H[n]=S1+family->cumulative_density(v[n])-family->cumulative_density(v[0]);
	double by_t=(t[n]-t[0])/s;
	double by_v=(v[n]-v[0])/s;

	for(int i=1;i<n;i++) {
		t[i]=t[i-1]+by_t;//printf("t[%d]=%lf\n",i,t[i]);
		v[i]=v[i-1]+by_v;
		h[i]=family->density(v[i]);
		H[i]=S1+family->cumulative_density(v[i])-family->cumulative_density(v[0]);
	}

	return DataFrame::create(
		_["t"]=NumericVector(t.begin(),t.end()),
		_["v"]=NumericVector(v.begin(),v.end()),
		_["i"]=NumericVector(h.begin(),h.end()),
		_["I"]=NumericVector(H.begin(),H.end())
	);
};

List VamModel::get_virtual_age_infos(double by) {

	// Only one system first!
	init_virtual_age_infos();
	int n=time.size() - 1;
	List res(n);
	while(k < n) {
		//printf("k=%d/n=%d\n",k,n);
		update_Vleft(false,false);
		res[k]=get_virtual_age_info(time[k],time[k+1],by);
		S1 += family->cumulative_density(Vleft) - family->cumulative_density(Vright);
		//gradient_update_for_current_system();
		int type2=type[k + 1];
		if(type2 < 0) type2=0;
		models->at(type2)->update(false,false);
	}
	return res;
};

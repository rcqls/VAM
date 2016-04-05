#include "rcpp_family_model.h"

// Constructor from R to build different kind of 
FamilyModel* newFamilyModel(List family) {
	std::string name=family["name"];
	NumericVector params=family["params"];
	FamilyModel*  fam=NULL;
	//DEBUG: printf("name=%s\n",name.c_str());

	if(name.compare("Weibull.family.cm") == 0) {
		//double alpha=params[0],beta=params[1];
		NumericVector par=NumericVector::create(1,3);
		if (params.size()==0){
            printf("WARNING: Weibull baseline hazard rate needs a parameter vector of length 2 ! It has been fixed to c(1,3). \n");
        } else if (params.size()==1){
            printf("WARNING: Weibull baseline hazard rate needs a parameter vector of length 2 ! It has been fixed to c(%f,3). \n",params[0]);
            par[1]=params[0];   
        } else if(params.size()!=2){
            printf("WARNING: Weibull baseline hazard rate needs a parameter vector of length 2 !\n");
            par[0]=params[0];
            par[1]=params[1];
        } else {
            par[0]=params[0];
            par[1]=params[1];
        }
		//DEBUG:printf("Params:alpha=%lf,beta=%lf\n",alpha,beta);
		fam=new WeibullFamilyModel(par);

	} else if(name.compare("LogLinear.family.cm") == 0) {
		//double alpha=params[0],beta=params[1];
		NumericVector par=NumericVector::create(1,3);
		if (params.size()==0){
            printf("WARNING: LogLinear baseline hazard rate needs a parameter vector of length 2 ! It has been fixed to c(1,3). \n");
        } else if (params.size()==1){
            printf("WARNING: LogLinear baseline hazard rate needs a parameter vector of length 2 ! It has been fixed to c(%f,3). \n",params[0]);
            par[1]=params[0];   
        } else if(params.size()!=2){
            printf("WARNING: LogLinear baseline hazard rate needs a parameter vector of length 2 !\n");
            par[0]=params[0];
            par[1]=params[1];
        } else {
            par[0]=params[0];
            par[1]=params[1];
        }
		fam=new LogLinearFamilyModel(par);

	} else if(name.compare("Weibull3.family.cm") == 0) {
		NumericVector par=NumericVector::create(1,3,0);
		if (params.size()==0){
            printf("WARNING: Weibull3 baseline hazard rate needs a parameter vector of length 3 ! It has been fixed to c(1,3,0). \n");
        } else if (params.size()==1){
            printf("WARNING: Weibull3 baseline hazard rate needs a parameter vector of length 3 ! It has been fixed to c(%f,3,0). \n",params[0]);
            par[1]=params[0];  
        } else if (params.size()==2){
            printf("WARNING: Weibull3 baseline hazard rate needs a parameter vector of length 3 ! It has been fixed to c(%f,%f,0). \n",params[0],params[1]);
            par[1]=params[0];
            par[1]=params[1];   
        } else if(params.size()!=3){
            printf("WARNING: Weibull3 baseline hazard rate needs a parameter vector of length 3 !\n");
            par[0]=params[0];
            par[1]=params[1];
            par[2]=params[2];
        } else {
            par[0]=params[0];
            par[1]=params[1];
            par[2]=params[2];
        }
		//double alpha=params[0],beta=params[1];
		fam=new Weibull3FamilyModel(par);

	} else {
    printf("WARNING: %s is not a proper baseline hazard rate!\n",name.c_str());
  }
	return fam;
}
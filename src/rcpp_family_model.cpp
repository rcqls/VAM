#include "rcpp_family_model.h"

// Constructor from R to build different kind of 
FamilyModel* newFamilyModel(List family) {
	std::string name=family["name"];
	NumericVector params=family["params"];
	FamilyModel*  fam=NULL;
	//DEBUG: printf("name=%s\n",name.c_str());
	if(name.compare("Weibull.family.cm") == 0) {
		//double alpha=params[0],beta=params[1];
		NumericVector par=NumericVector::create(params[0],params[1]);//LD3
		//DEBUG:printf("Params:alpha=%lf,beta=%lf\n",alpha,beta);
		fam=new WeibullFamilyModel(par);//LD3
	} else if(name.compare("LogLinear.family.cm") == 0) {
		//double alpha=params[0],beta=params[1];
		NumericVector par=NumericVector::create(params[0],params[1]);//LD3
		fam=new LogLinearFamilyModel(par);//LD3
	} else {
    printf("WARNING: %s is not a proper baseline hazard rate!\n",name.c_str());
  }
	return fam;
}
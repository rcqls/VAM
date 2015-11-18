#include "rcpp_family_model.h"

// Constructor from R to build different kind of 
FamilyModel* newFamilyModel(List family) {
	std::string name=family["name"];
	NumericVector params=family["params"];
	FamilyModel*  fam=NULL;
	//DEBUG: printf("name=%s\n",name.c_str());
	if(name.compare("Weibull.family.cm") == 0) {
		double alpha=params[0],beta=params[1];
		//DEBUG:printf("Params:alpha=%lf,beta=%lf\n",alpha,beta);
		fam=new WeibullFamilyModel(alpha,beta);
	} else if(name.compare("LogLinear.family.cm") == 0) {
		double alpha=params[0],beta=params[1];
		fam=new LogLinearFamilyModel(alpha,beta);
	}	
	return fam;
}
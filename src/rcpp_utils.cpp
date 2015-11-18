#include <Rcpp.h>
using namespace Rcpp ;

//TO detect if an external pointer is nil
SEXP is_xptr_null(SEXP objR) {
	void *objPtr;
	bool ans;
	objPtr=R_ExternalPtrAddr(objR);
  	if(objPtr==NULL) {
  		ans=1;
	} else {
		ans=0;
	}
	return Rcpp::wrap(ans);
}

Environment envFromListWithParent( List li, Environment parent) {

	Environment env=parent.new_child(true);
	CharacterVector noms=li.names();
	int l=li.size();

    for(int i = 0; i < l ; i++) {
		SEXP name = Rf_install(Rf_translateChar(STRING_ELT(noms, i)));
		Rf_defineVar(name, VECTOR_ELT(li, i), env);
    }

    return env ;
}


Environment envFromList( List li) {
    return envFromListWithParent(li,Environment::global_env());
}

//Not necessarily useful since mainly used in C++
RCPP_MODULE(tools_module) {
	function( "is_xptr_null", &is_xptr_null );
	function( "envFromListWithParent", &envFromListWithParent );
	function( "envFromList", &envFromList );
}
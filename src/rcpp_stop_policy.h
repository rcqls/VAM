#ifndef RCPP_STOP_POLICY_H
#define RCPP_STOP_POLICY_H
#include <Rcpp.h>
#include "rcpp_sim_vam.h"

using namespace Rcpp ;
class VamModel;

class StopPolicy {
public:

	StopPolicy(SimVam* sim_) {};

	virtual ~StopPolicy() {};

	virtual bool check() = 0;

	 

private:

    SimVam* sim;

};


class AfterNumberEventStopPolicy : public StopPolicy {
public:
    AfterNumberEventStopPolicy(SimVam* sim_): StopPolicy(sim_) {
        
    }

    ~AfterNumberEventStopPolicy() {};

    
 
};

#endif //RCPP_STOP_POLICY_H
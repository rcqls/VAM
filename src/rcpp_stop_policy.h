#ifndef RCPP_STOP_POLICY_H
#define RCPP_STOP_POLICY_H
#include <Rcpp.h>
#include "rcpp_sim_vam.h"

using namespace Rcpp ;
class SimVam;

class StopPolicy {
public:

	StopPolicy(SimVam* sim_) {sim=sim_;};

	virtual ~StopPolicy() {};

	virtual bool ok() = 0;

    SimVam* sim;

};


class AtRunStopPolicy : public StopPolicy {
public:
    AtRunStopPolicy(SimVam* sim_,int nb_): StopPolicy(sim_) {
        nb=nb_;
    }

    ~AtRunStopPolicy() {};

    bool ok();

    int nb;
    
 
};

class AtTimeStopPolicy : public StopPolicy {
public:
    AtTimeStopPolicy(SimVam* sim_,double time_): StopPolicy(sim_) {
        time=time_;
    }

    ~AtTimeStopPolicy() {};

    bool ok();

    double time;
    
 
};

#endif //RCPP_STOP_POLICY_H
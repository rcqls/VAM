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

	virtual void first() = 0;

    SimVam* sim;

};

StopPolicy* newStopPolicy(SimVam* sim,List policy);

class SizeGreaterThanStopPolicy : public StopPolicy {
public:
    SizeGreaterThanStopPolicy(SimVam* sim_,int size_): StopPolicy(sim_) {
        size=size_;
    }

    ~SizeGreaterThanStopPolicy() {};

    bool ok();

		void first() {};

    int size;


};

//M is for CM or PM, type=-1 or 1,2,...
class SizeOfTypeGreaterThanStopPolicy : public StopPolicy {
public:
    SizeOfTypeGreaterThanStopPolicy(SimVam* sim_,int type_,int size_): StopPolicy(sim_) {
        type=type_;
        size=size_;
        count=0;
    }

    ~SizeOfTypeGreaterThanStopPolicy() {};

    bool ok();

		void first() {};

    int size;

    int type;

    int count;

};

class TimeGreaterThanCensorshipStopPolicy : public StopPolicy {
public:
    TimeGreaterThanCensorshipStopPolicy(SimVam* sim_,double time_,Language expr_,Environment env_): StopPolicy(sim_) {
        time=time_;expr=expr_;env=env_;
				to_init=(time<0);
    }

    ~TimeGreaterThanCensorshipStopPolicy() {};

    bool ok();

		void first();

		bool to_init;
		
    double time;

		Language expr;

		Environment env;


};

class TimeGreaterThanStopPolicy : public StopPolicy {
public:
    TimeGreaterThanStopPolicy(SimVam* sim_,double time_): StopPolicy(sim_) {
        time=time_;
    }

    ~TimeGreaterThanStopPolicy() {};

    bool ok();

		void first() {};

    double time;


};

class AndStopPolicy : public StopPolicy {
public:
    AndStopPolicy(SimVam* sim_,List policies_): StopPolicy(sim_) {
        for(
            List::iterator it=policies_.begin();
            it != policies_.end();
            ++it
        ) {
            List policy=*it;
            StopPolicy*  sp=newStopPolicy(sim_,policy);
            if(!(sp == NULL)) policies.push_back(sp);
        }
    }

    ~AndStopPolicy() {
        for(
            std::vector<StopPolicy*>::iterator it=policies.begin();
            it != policies.end();
            ++it
        ) {
            delete *it;
        }
    };

    bool ok();

		void first();

    std::vector<StopPolicy*> policies;


};

class OrStopPolicy : public StopPolicy {
public:
    OrStopPolicy(SimVam* sim_,List policies_): StopPolicy(sim_) {
        for(
            List::iterator it=policies_.begin();
            it != policies_.end();
            ++it
        ) {
            List policy=*it;
            StopPolicy*  sp=newStopPolicy(sim_,policy);
            if(!(sp == NULL)) policies.push_back(sp);
        }
    }

    ~OrStopPolicy() {
        for(
            std::vector<StopPolicy*>::iterator it=policies.begin();
            it != policies.end();
            ++it
        ) {
            delete *it;
        }
    };

    bool ok();

		void first();

    std::vector<StopPolicy*> policies;


};

#endif //RCPP_STOP_POLICY_H

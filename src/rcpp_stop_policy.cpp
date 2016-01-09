#include "rcpp_stop_policy.h"
#include "rcpp_sim_vam.h"

StopPolicy* newStopPolicy(SimVam* sim,List policy) {
    StopPolicy* sp=NULL;
    std::string name=policy["name"];
    //DEBUG:printf("name=%s\n",name.c_str());
    if(name.compare("SizeGreaterThan.stop.policy") == 0) {
        //DEBUG:printf("Params:alpha=%lf,beta=%lf\n",alpha,beta);
        int size_=policy["size"];
        sp=new SizeGreaterThanStopPolicy(sim,size_);
    } else if(name.compare("SizeOfTypeGreaterThan.stop.policy") == 0) {
        //DEBUG:printf("Params:alpha=%lf,beta=%lf\n",alpha,beta);
        int size_=policy["size"];
        int type_=policy["type"];
        sp=new SizeOfTypeGreaterThanStopPolicy(sim,type_,size_);
    } else if(name.compare("TimeGreaterThanCensorship.stop.policy") == 0) {
        double time_=policy["time"];
        List expr_=policy["time.expr"];
        Language l_=expr_["expr"];
        Environment env_=expr_["env"];
        sp=new TimeGreaterThanCensorshipStopPolicy(sim,time_,l_,env_);
    } else if(name.compare("TimeGreaterThan.stop.policy") == 0) {
        double time_=policy["time"];
        sp=new TimeGreaterThanStopPolicy(sim,time_);
    } else if(name.compare("And.stop.policy") == 0) {
        List policies_=policy["policies"];
        sp=new AndStopPolicy(sim,policies_);
    } else if(name.compare("Or.stop.policy") == 0) {
        List policies_=policy["policies"];
        sp=new OrStopPolicy(sim,policies_);
    }
    return sp;
}

bool SizeGreaterThanStopPolicy::ok() {
    return (sim->get_model()->k < size);
}

bool SizeOfTypeGreaterThanStopPolicy::ok() {
    VamModel* mod=sim->get_model();
    //incr counter
    //printf("k=%d,t=%lf,ty=%d ->",mod->k,mod->time[mod->k],mod->type[mod->k]);
    if(mod->k>0 && mod->type[mod->k] == type) {
        count++;
        //printf("type=%d,count=%d",type,count);
    }
    //printf("\n");
    return (count < size);
}

bool TimeGreaterThanCensorshipStopPolicy::ok() {
    VamModel* mod=sim->get_model();
    bool ok=(mod->time[mod->k] < time);
    if(!ok) {
        //update result of sim!
        mod->time[mod->k]=time;
        mod->type[mod->k]=0;
    }
    return ok;
}

void TimeGreaterThanCensorshipStopPolicy::first() {
    //printf("time=%lf\n",time);
     if(to_init) {
       time=as<double>(Rf_eval(expr,env));
     };
}

//Same as AtTime but the next time is not missing...
bool TimeGreaterThanStopPolicy::ok() {
    VamModel* mod=sim->get_model();
    bool ok=(mod->time[mod->k] < time);
    return ok;
}


bool AndStopPolicy::ok() {
    bool ans=false;
    for(
        std::vector<StopPolicy*>::iterator it=policies.begin();
        it != policies.end();
        ++it
    ) {
        if((*it)->ok()) ans |= true; //every cond is tested because some init is done there!
    }
    return ans;
}

bool OrStopPolicy::ok() {
    bool ans=true;
    for(
        std::vector<StopPolicy*>::iterator it=policies.begin();
        it != policies.end();
        ++it
    ) {
        if(!(*it)->ok()) ans &= false; //every cond is tested because some init is done there!
    }
    return ans;
}


// Exactly the same first method for the following classes
void AndStopPolicy::first() {
    for(
        std::vector<StopPolicy*>::iterator it=policies.begin();
        it != policies.end();
        ++it
    ) {
        (*it)->first();
    }
}

void OrStopPolicy::first() {
    for(
        std::vector<StopPolicy*>::iterator it=policies.begin();
        it != policies.end();
        ++it
    ) {
        (*it)->first();
    }
}

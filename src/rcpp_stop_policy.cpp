#include "rcpp_stop_policy.h"
#include "rcpp_sim_vam.h"

StopPolicy* newStopPolicy(SimVam* sim,List policy) {
    StopPolicy* sp=NULL;
    std::string name=policy["name"];
    //DEBUG:printf("name=%s\n",name.c_str());
    if(name.compare("AtSize.stop.policy") == 0) {
        //DEBUG:printf("Params:alpha=%lf,beta=%lf\n",alpha,beta);
        int size_=policy["size"];
        sp=new AtSizeStopPolicy(sim,size_);
    } else if(name.compare("AtMSize.stop.policy") == 0) {
        //DEBUG:printf("Params:alpha=%lf,beta=%lf\n",alpha,beta);
        int size_=policy["size"];
        int type_=policy["type"];
        sp=new AtMSizeStopPolicy(sim,type_,size_);
    } else if(name.compare("AtTime.stop.policy") == 0) {
        double time_=policy["time"];
        sp=new AtTimeStopPolicy(sim,time_);
    } else if(name.compare("And.stop.policy") == 0) {
        List policies_=policy["policies"];
        sp=new AndStopPolicy(sim,policies_);
    } else if(name.compare("Or.stop.policy") == 0) {
        List policies_=policy["policies"];
        sp=new OrStopPolicy(sim,policies_);
    }
    return sp;
}

bool AtSizeStopPolicy::ok() {
    return (sim->get_model()->k < size);
}

bool AtMSizeStopPolicy::ok() {
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

bool AtTimeStopPolicy::ok() {
    VamModel* mod=sim->get_model();
    bool ok=(mod->time[mod->k] < time);
    if(!ok) {
        //update result of sim!
        mod->time[mod->k]=time;
        mod->type[mod->k]=0;
    }
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

#include "rcpp_sim_vam.h"
#include "rcpp_stop_policy.h"

DataFrame SimVam::simulate(int nbsim) {
    init(nbsim);

    while(stop_policy->ok()) {//model->k < nbsim) {
        //### modAV <- if(Type[k]<0) obj$vam.CM[[1]]$model else obj$vam.PM$models[[obj$data$Type[k]]]
        //# Here, obj$model$k means k-1
        //#print(c(obj$model$Vleft,obj$model$Vright))
        double timePM= 0.0, timeCM = model->models->at(model->idMod)->virtual_age_inverse(model->family->inverse_cumulative_density(model->family->cumulative_density(model->models->at(model->idMod)->virtual_age(model->time[model->k]))-log(runif(1))[0]));
        //TODO: submodels
        int idMod;
        List timeAndTypePM;
        if(model->maintenance_policy != NULL) {
            timeAndTypePM = model->maintenance_policy->update(model); //# Peut-être ajout Vright comme argument de update
            //timeAndTypePM = model->maintenance_policy->update(model->time[model->k]); //# Peut-être ajout Vright comme argument de update

            NumericVector tmp=timeAndTypePM["time"];
            timePM=tmp[0];
        }
        if(model->maintenance_policy == NULL || timeCM < timePM) {
            model->time[model->k + 1]=timeCM;
            model->type[model->k + 1]=-1;
            idMod=0;
        } else {
            model->time[model->k + 1]=timePM;
            NumericVector tmp2=timeAndTypePM["type"];
            int typePM=tmp2[0];
            model->type[model->k + 1]=typePM;
            idMod=timeAndTypePM["type"];
        }
        //printf("k=%d: cm=%lf,pm=%lf\n",model->k,timeCM,timePM);
        //# used in the next update
        model->update_Vleft(false);
        
        //# update the next k, and save model in model too!
        model->models->at(idMod)->update(false);
        

        //To dynamically increase the size of simulation
        resize();
    }

    return get_data();
}

void SimVam::add_stop_policy(List policy) {
    std::string name=policy["name"];
    //DEBUG:printf("name=%s\n",name.c_str());
    if(name.compare("AtRun.stop.policy") == 0) {
        //DEBUG:printf("Params:alpha=%lf,beta=%lf\n",alpha,beta);
        int nb=policy["nb"];
        stop_policy=new AtRunStopPolicy(this,nb);
    } else if(name.compare("AtTime.stop.policy") == 0) {
        double time=policy["time"];
        stop_policy=new AtTimeStopPolicy(this,time);
    }

}

void SimVam::init(int cache_size_) {
    model->Vright=0;
    model->k=0;
    cache_size=cache_size_;
    model->idMod=0; // Since no maintenance is possible!
    model->time=rep(0,cache_size+1);
    model->type= rep(1,cache_size+1);
}

void SimVam::resize() {
    int inf=std::numeric_limits<int>::max();
}

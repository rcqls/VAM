#include "rcpp_sim_vam.h"
#include "rcpp_stop_policy.h"


DataFrame SimVam::get_last_data() {
    //printf("size:%d,%d\n",(model->time).size(),model->k+1);
    if((model->time).size() > model->k+1) {
        size = model->k+1;
        //printf("data size:%d\n",size);
        (model->time).resize(size);
        (model->type).resize(size);
    }
    return DataFrame::create(_["Time"]=model->time,_["Type"]=model->type);
}

DataFrame SimVam::simulate(int nbsim) {
    init(nbsim);

    stop_policy->first();

    if(model->maintenance_policy != NULL) model->maintenance_policy->first();

    while(stop_policy->ok()) {//model->k < nbsim) {
        //printf("k=%d\n",model->k);
        //To dynamically increase the size of simulation
        resize();
        //printf("k2=%d\n",model->k);

        //### modAV <- if(Type[k]<0) obj$vam.CM[[1]]$model else obj$vam.PM$models[[obj$data$Type[k]]]
        //# Here, obj$model$k means k-1
        //#print(c(obj$model$Vleft,obj$model$Vright))

        RNGScope rngScope;

        double timePM= 0.0, timeCM = model->virtual_age_inverse(model->family->inverse_cumulative_hazardRate(model->family->cumulative_hazardRate(model->virtual_age(model->time[model->k]))-log(runif(1))[0]));
        //TODO: submodels
        int idMod;
        List timeAndTypePM;
        if(model->maintenance_policy != NULL) {
            timeAndTypePM = model->maintenance_policy->update(model); //# Peut-être ajout Vright comme argument de update
            //timeAndTypePM = model->maintenance_policy->update(model->time[model->k]); //# Peut-être ajout Vright comme argument de update

            NumericVector tmp=timeAndTypePM["time"];
            timePM=tmp[0];
            //DEBUG: printf("sim: timePM:%lf, timeCM=%lf\n",timePM,timeCM);
            if(timePM<timeCM && timePM<model->time[model->k]) {
          		printf("Warning: PM ignored since next_time(=%lf)<current_time(=%lf) at rank %d.\n",timePM,model->time[model->k],model->k);
            }
        }
        if(model->maintenance_policy == NULL || timeCM < timePM || timePM<model->time[model->k]) {
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
        model->update_Vleft(false,false);


        //# update the next k, and save model in model too!
        model->models->at(idMod)->update(false,false);

    }

    return get_last_data();
}

void SimVam::add_stop_policy(List policy) {
    stop_policy=newStopPolicy(this,policy);
}

void SimVam::init(int cache_size_) {
    // Almost everything in the 5 following lines are defined in model->init_computation_values() (but this last one initializes more than this 5 lines)
    model->Vright=0;
    model->A=1;
    model->C=0;
    model->k=0;
    for(int i=0;i<model->nbPM + 1;i++) model->models->at(i)->init();

    size=cache_size_+1;cache_size=cache_size_;
    model->idMod=0; // Since no maintenance is possible!
    //model->time=rep(0,size);
    //model->type= rep(1,size);
    (model->time).clear();
    (model->type).clear();
    (model->time).resize(size,0);
    (model->type).resize(size,1);
}

#define print_vector(x)                                                                     \
    for (std::vector<double>::const_iterator i = x.begin(); i != x.end(); ++i)   \
    std::cout << *i << ' '; \
    std::cout << "\n";

void SimVam::resize() {
    if(model->k > size-2) {
        //printf("RESIZE!\n");
        //print_vector((model->time))
        //printf("SIZE=%d",size);
        size += cache_size;//printf("->%d\n",size);
        (model->time).resize(size);
        //print_vector((model->time))
        //printf("model->SIZE=%d\n",(model->time).size());
        (model->type).resize(size);
    }
}

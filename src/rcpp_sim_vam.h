#ifndef RCPP_SIM_VAM_H
#define RCPP_SIM_VAM_H
#include <Rcpp.h>
#include "rcpp_maintenance_model.h"
#include "rcpp_maintenance_policy.h"

using namespace Rcpp ;

class SimVam { 

public:

    SimVam(List model_) {
        model=new VamModel(model_);
    };

    ~SimVam() {
        //printf("simvam final\n");
        delete model;
    };

    DataFrame get_data() {
        return DataFrame::create(_["Time"]=model->time,_["Type"]=model->type);
    }

    DataFrame simulate(int nbsim) {
        init(nbsim);

        while(model->k < nbsim) {
             
            //Rprintf("",model->models.size())

            /*** This is a version of the synthetic expression of timeCM below
            double timePM;
            double tmp0=model->time[model->k];
            double tmp1=model->models->at(model->idMod)->virtual_age(tmp0);
            double tmp2=model->family->cumulative_density(tmp1)-log(runif(1))[0];
            double tmp3=model->family->inverse_cumulative_density(tmp2);
            double tmp4 = model->models->at(model->idMod)->virtual_age_inverse(tmp3);
            double timeCM = tmp4; ***/

            double timePM, timeCM = model->models->at(model->idMod)->virtual_age_inverse(model->family->inverse_cumulative_density(model->family->cumulative_density(model->models->at(model->idMod)->virtual_age(model->time[model->k]))-log(runif(1))[0]));
            int idMod;
             
            std::pair<double,int> timeAndTypePM;
            if(model->maintenance_policy != NULL) {
                timeAndTypePM = model->maintenance_policy->update(model);
                 
                timePM=timeAndTypePM.first; //time
            }

            if(model->maintenance_policy == NULL || timeCM < timePM) {
                model->time[model->k + 1]=timeCM;
                model->type[model->k + 1]=-1;
                idMod=0;
            } else {
                model->time[model->k + 1]=timePM;
                int typePM=timeAndTypePM.second; //type
                model->type[model->k + 1]=typePM;
                idMod=typePM;
            }
            //# used in the next update
            model->update_Vleft(false);
            //# update the next k, and save model in model too!
            model->models->at(idMod)->update(false);
        }

        return get_data();
    }

    VamModel* get_model() {
        return model;
    }

    NumericVector get_params() {
        return model->get_params();
    }

    void set_params(NumericVector pars) {
        model->set_params(pars);
    }

    //delegate from model cache!
    List get_virtual_age_infos(double by) {
        return model->get_virtual_age_infos(by);
    }


private:
    VamModel* model;

    void init(int nbsim) {
        model->Vright=0;
        model->k=0;
        model->idMod=0; // Since no maintenance is possible!
        (model->time).resize(nbsim+1,0); //equivalent R: rep(0,nbsim+1)
        (model->type).resize(nbsim+1,0); //equivalent R: rep(1,nbsim+1)
    } 

};

#endif //RCPP_SIM_VAM_H
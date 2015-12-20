#include "rcpp_stop_policy.h"
#include "rcpp_sim_vam.h"

bool AtRunStopPolicy::ok() {
    return (sim->get_model()->k < nb);
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

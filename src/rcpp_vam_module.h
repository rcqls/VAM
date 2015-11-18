#ifndef RCPP_VAM_MODULE_H
#define RCPP_VAM_MODULE_H
#include "rcpp_sim_vam.h"
#include "rcpp_mle_vam.h"

RCPP_EXPOSED_AS(SimVam);
RCPP_EXPOSED_WRAP(SimVam);

RCPP_EXPOSED_AS(MLEVam);
RCPP_EXPOSED_WRAP(MLEVam);

#endif //RCPP_VAM_MODULE_H
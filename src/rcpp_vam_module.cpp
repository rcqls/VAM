#include "rcpp_vam_module.h"

// static void finalizer_of_sim_vam( SimVam* ptr ){
//       printf("finalizer sim_vam has been called\n");
// }

// static void finalizer_of_mle_vam( MLEVam* ptr ){
//       printf("finalizer mle_vam %p has been called\n",ptr);
// }

RCPP_MODULE(vam_module) {
	class_<VamModel>("ModelVam")
	.constructor<List>()
	.constructor<List,List>()
	.field( "time", &VamModel::time, "time" )
	.field( "type", &VamModel::type, "type" )
	.method( "get", &VamModel::get, "get many informations" )
	.method( "family", &VamModel::get_family, "get family" )
	.method("get_params",&VamModel::get_params,"get params")
    .method("set_params",&VamModel::set_params,"set params")
    .method("get_virtual_age_infos",&VamModel::get_virtual_age_infos,"get infos related to virtual ages")
	.method("set_data",&VamModel::set_data,"set data")
    .method("get_data",&VamModel::get_selected_data,"get (selected) data")
    ;

    class_<SimVam>( "SimVam" )
    .constructor<List>()
    //.finalizer( &finalizer_of_sim_vam)
    .method("model",&SimVam::get_model,"model accessor")
    .method("simulate",&SimVam::simulate,"simulate")
    .method("get_params",&SimVam::get_params,"get params")
    .method("set_params",&SimVam::set_params,"set params")
    .method("get_data",&SimVam::get_data,"get data")
    .method("get_virtual_age_infos",&SimVam::get_virtual_age_infos,"get infos related to virtual ages")
    ;

    class_<MLEVam>( "MLEVam" )
    .constructor<List,List>()
    //.finalizer( &finalizer_of_mle_vam)
    .method("model",&MLEVam::get_model,"model accessor")
    .method("set_data",&MLEVam::set_data,"set data")
    //.method("get_selected_data",&MLEVam::get_selected_data,"get selected data")
    .method("contrast",&MLEVam::contrast,"compute contrast")
    .method("gradient",&MLEVam::gradient,"compute gradient")
    .method("alpha_est",&MLEVam::get_alpha_est,"get alpha estimation")
    .method("get_params",&MLEVam::get_params,"get params")
    .method("set_params",&MLEVam::set_params,"set params")
    .method("get_virtual_age_infos",&MLEVam::get_virtual_age_infos,"get infos related to virtual ages")
    .method("get_data",&MLEVam::get_selected_data,"get (selected) data")
    ;

    class_<FamilyModel>("FamilyModel")
    .method("density",&FamilyModel::density,"density")
    .method("cumulative_density",&FamilyModel::cumulative_density,"cumulative density")
    .method("density_derivative",&FamilyModel::density_derivative,"density_derivative")
    .method("inverse_cumulative_density",&FamilyModel::inverse_cumulative_density,"inverse cumulative density")
    .method("density_param_derivative",&FamilyModel::density_param_derivative,"density derivative with respect to beta")
    .method("cumulative_density_param_derivative",&FamilyModel::cumulative_density_param_derivative,"cumulative density derivative with respect to beta")
    ;

}
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
	//.method( "get", &VamModel::get, "get many informations" )
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
	//.method("model",&SimVam::get_model,"model accessor")
	.method("simulate",&SimVam::simulate,"simulate")
	.method("get_params",&SimVam::get_params,"get params")
	.method("set_params",&SimVam::set_params,"set params")
	.method("get_data",&SimVam::get_selected_data,"get (selected) data")
	.method("set_data",&SimVam::set_data,"set data")
	.method("get_virtual_age_infos",&SimVam::get_virtual_age_infos,"get infos related to virtual ages")
	.method("add_stop_policy",&SimVam::add_stop_policy,"add stop policy")
	;

	class_<MLEVam>( "MLEVam" )
	.constructor<List,List>()
	//.finalizer( &finalizer_of_mle_vam)
	//.method("model",&MLEVam::get_model,"model accessor")
	.method("set_data",&MLEVam::set_data,"set data")
	.method("set_leftCensors",&MLEVam::set_leftCensors,"set leftCensors")
	//.method("get_selected_data",&MLEVam::get_selected_data,"get selected data")
	.method("contrast",&MLEVam::contrast,"compute contrast")
	.method("gradient",&MLEVam::gradient,"compute gradient")
	.method("hessian",&MLEVam::hessian,"compute hessian")
	.method("alpha_est",&MLEVam::get_alpha_est,"get alpha estimation")
	.method("get_params",&MLEVam::get_params,"get params")
	.method("set_params",&MLEVam::set_params,"set params")
	.method("get_virtual_age_infos",&MLEVam::get_virtual_age_infos,"get infos related to virtual ages")
	.method("get_data",&MLEVam::get_selected_data,"get (selected) data")
	;

	class_<BayesianVam>( "BayesianVam" )
	.constructor<List,List,List>()
	//.finalizer( &finalizer_of_bayesian_vam)
	.method("get_params",&BayesianVam::get_params,"get params")
	.method("set_params",&BayesianVam::set_params,"set params")
	.method("set_sigma",&BayesianVam::set_sigma_at,"set sigma")
	.method("get_sigma",&BayesianVam::get_sigma,"get sigma")
	.method("get_prior_params",&BayesianVam::get_prior_params,"get prior params")
	.method("get_virtual_age_infos",&BayesianVam::get_virtual_age_infos,"get infos related to virtual ages")
	.method("get_data",&BayesianVam::get_selected_data,"get (selected) data")
	.method("mcmc",&BayesianVam::mcmc,"mcmc method")
	.method("mcmc_history",&BayesianVam::mcmc_history,"mcmc method with result as data.frame")
	;

	class_<FamilyModel>("FamilyModel")
	.method("hazardRate",&FamilyModel::hazardRate,"hazard rate")
	.method("cumulative_hazardRate",&FamilyModel::cumulative_hazardRate,"cumulative hazardrd rate")
	.method("hazardRate_derivative",&FamilyModel::hazardRate_derivative,"hazard rate derivative")
	.method("hazardRate_2derivative",&FamilyModel::hazardRate_2derivative,"hazard rate second order derivative")
	.method("inverse_cumulative_hazardRate",&FamilyModel::inverse_cumulative_hazardRate,"inverse cumulative hazard rate")
	.method("hazardRate_param_derivative",&FamilyModel::hazardRate_param_derivative,"hazard rate derivative with respect to beta")
	.method("cumulative_hazardRate_param_derivative",&FamilyModel::cumulative_hazardRate_param_derivative,"cumulative hazard rate derivative with respect to beta")
	.method("hazardRate_param_2derivative",&FamilyModel::hazardRate_param_2derivative,"hazard rate second order derivative with respect to beta")
	.method("cumulative_hazardRate_param_2derivative",&FamilyModel::cumulative_hazardRate_param_2derivative,"cumulative hazard rate second order derivative with respect to beta")
	;

	function( "newMaintenancePolicy", &newMaintenancePolicy );

	class_<MaintenancePolicy>("MaintenancePolicy")
	.method("update",&MaintenancePolicy::update,"update")
	.method("get_params",&MaintenancePolicy::get_params,"get params")
	.method("set_params",&MaintenancePolicy::set_params,"set params")
	;

}

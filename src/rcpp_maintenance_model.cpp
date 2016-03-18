#include "rcpp_maintenance_model.h"

using namespace Rcpp ;

//Forward declarations
//class VamModel;
//class MaintenanceModel;
//MaintenanceModel* newMaintenanceModel(List model,VamModel* model);

//Effective declarations
MaintenanceModelList::MaintenanceModelList(List models_,VamModel* model) {
    int i=0;
    int j=0;
    for(
        List::iterator lit=models_.begin();
        lit != models_.end();
        ++lit
    ) {
    	List maintenance=*lit;
    	MaintenanceModel*  vam=newMaintenanceModel(maintenance,model);
        if(!(vam == NULL)) {
            vam->set_id(i++);
            vam->set_id_params(j);
            j+=vam->nb_params();
            model_list.push_back(vam);
        }
    }
}

MaintenanceModelList::~MaintenanceModelList() {
	for(
		std::vector<MaintenanceModel*>::iterator vit=model_list.begin();
		vit != model_list.end();
        ++vit
    ) {
		delete *vit;
	}

}

void ARA1::update(bool with_gradient,bool with_hessian) {
    /*# next step
    obj$vam$model$k <- obj$vam$model$k + 1
    # At T(k)
    obj$vam$model$Vright <- obj$vam$model$Vright + (1-obj$rho)*(dVlr <-(obj$vam$model$Vleft-obj$vam$model$Vright))
    if(with.gradient) {
        # only the rho parameters
        #obj$vam$model$dVright <- obj$vam$model$dVright + rep(0,1+length(obj$vam$vam.PM$models))
        i <- match(obj$id,seq(obj$vam$vam.PM$models),nomatch=0)+1
        obj$vam$model$dVright[i] <- obj$vam$model$dVright[i] - dVlr
    }
    # save old model
    obj$model$mod <- obj
    */
    model->k += 1;
    double dVlr = model->Vleft- model->Vright;
    model->Vright += (1-rho) * dVlr;
    if(with_gradient) {
        model->dVright[id_params] += -dVlr;
    }
    model->idMod = id;
}

double ARA1::virtual_age(double time) {
    //max(0.0000001,obj$vam$model$Vright+time-obj$vam$data$Time[obj$vam$model$k])
    //printf("virtual_age:%lf,%lf,%lf\n",model -> Vright, time,model->time[model->k]);
    return model -> Vright + time  - model->time[model->k];
}

double* ARA1::virtual_age_derivative(double x) {
    return model->dVright;
}

double* ARA1::virtual_age_hessian(double x) {
    return model->d2Vright;
}

double ARA1::virtual_age_inverse(double time) {
    return time + model->time[model->k] - model->Vright;
}

void ARAInf::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    model->k += 1;
    model->Vright = (1-rho) * model->Vleft;
    if (with_hessian){
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2Vright[i*(i+1)/2+j] = (1-rho) * model->d2Vright[i*(i+1)/2+j];
            }
        }
        for(j=0;j<=id_params;j++) {
            //i(<=id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
            model->d2Vright[id_params*(id_params+1)/2+j] = model->d2Vright[id_params*(id_params+1)/2+j] - model->dVleft[j];
        }
        for(i=id_params;i<model->nb_paramsMaintenance;i++) {
             //id and i(>=id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            model->d2Vright[i*(i+1)/2+id_params] = model->d2Vright[i*(i+1)/2+id_params] - model->dVleft[i];
        }
    }
    if(with_gradient||with_hessian) {
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dVright[i] = (1-rho) * model->dVright[i];
        }
        model->dVright[id_params] = model->dVright[id_params] - model->Vleft;
    }
    // save old model
    model->idMod = id;
}

double ARAInf::virtual_age(double time) {
    //max(0.0000001,obj$vam$model$Vright+time-obj$vam$data$Time[obj$vam$model$k])
    return model -> Vright + time  - model->time[model->k];
}

double* ARAInf::virtual_age_derivative(double x) {
    return model->dVright;
}

double* ARAInf::virtual_age_hessian(double x) {
    return model->d2Vright;
}

double ARAInf::virtual_age_inverse(double time) {
    return time + model->time[model->k] - model->Vright;
}

void AGAN::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    model->k += 1;
    model->Vright = 0;
    if (with_hessian){
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2Vright[i*(i+1)/2+j] = 0;
            }
        }
    }
    if(with_gradient||with_hessian) {
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dVright[i] = 0;
        }
    }
    // save old model
    model->idMod = id;
}

double AGAN::virtual_age(double time) {
    //max(0.0000001,obj$vam$model$Vright+time-obj$vam$data$Time[obj$vam$model$k])
    return model -> Vright + time  - model->time[model->k];
}

double* AGAN::virtual_age_derivative(double x) {
    return model->dVright;
}

double* AGAN::virtual_age_hessian(double x) {
    return model->d2Vright;
}

double AGAN::virtual_age_inverse(double time) {
    return time + model->time[model->k] - model->Vright;
}

void ABAO::update(bool with_gradient,bool with_hessian) {
    model->k += 1;
    model->Vright += model->Vleft- model->Vright;

    // save old model
    model->idMod = id;
}

double ABAO::virtual_age(double time) {
    //max(0.0000001,obj$vam$model$Vright+time-obj$vam$data$Time[obj$vam$model$k])
    return model -> Vright + time  - model->time[model->k];
}

double* ABAO::virtual_age_derivative(double x) {
    return model->dVright;
}

double* ABAO::virtual_age_hessian(double x) {
    return model->d2Vright;
}

double ABAO::virtual_age_inverse(double time) {
    return time + model->time[model->k] - model->Vright;
}

MaintenanceModel* newMaintenanceModel(List maintenance,VamModel* model) {
	std::string name=maintenance["name"];
	NumericVector params=maintenance["params"];
	MaintenanceModel*  mm=NULL;
	if(name.compare("ARA1.va.model") == 0) {
		//double rho=params[0];
        NumericVector rho=NumericVector::create(params[0]);
		mm=new ARA1(rho,model);
	} else if(name.compare("ARAInf.va.model") == 0) {
		//double rho=params[0];
        NumericVector rho=NumericVector::create(params[0]);
		mm=new ARAInf(rho,model);
	} else if(name.compare("AGAN.va.model") == 0) {
    //double rho=1.0;
    //mm=new ARAInf(rho,model);
        mm=new AGAN(model);
  } else if(name.compare("ABAO.va.model") == 0) {
    //double rho=0.0;
    //mm=new ARAInf(rho,model);
        mm=new ABAO(model);
  } else {
    printf("WARNING: %s is not a proper maintenance model!\n",name.c_str());
  }
	return mm;
}

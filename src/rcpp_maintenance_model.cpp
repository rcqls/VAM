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
    int i;
    int j;
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

    if (with_hessian){
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2Vright[i*(i+1)/2+j] = (1-rho) * model->d2Vleft[i*(i+1)/2+j] + rho * model->d2Vright[i*(i+1)/2+j];
            }
        }
        for(j=0;j<=id_params;j++) {
            //i(<=id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
            model->d2Vright[id_params*(id_params+1)/2+j] = model->d2Vright[id_params*(id_params+1)/2+j] - model->dVleft[j]+model->dVright[j];
        }
        for(i=id_params;i<model->nb_paramsMaintenance;i++) {
             //id and i(>=id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            model->d2Vright[i*(i+1)/2+id_params] = model->d2Vright[i*(i+1)/2+id_params] - model->dVleft[i]+model->dVright[i];
        }
    }
    if(with_gradient||with_hessian) {
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dVright[i] = (1-rho) * model->dVleft[i]+ rho * model->dVright[i];
        }
        model->dVright[id_params] = model->dVright[id_params] - model->Vleft + model->Vright;
    }

    model->Vright= (1-rho)*model->Vleft + rho*model->Vright;

    model->idMod = id;
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
                model->d2Vright[i*(i+1)/2+j] = (1-rho) * model->d2Vleft[i*(i+1)/2+j];
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
            model->dVright[i] = (1-rho) * model->dVleft[i];
        }
        model->dVright[id_params] = model->dVright[id_params] - model->Vleft;
    }
    // save old model
    model->idMod = id;
}

void AGAN::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    model->k += 1;
    model->Vright = 0;
    model->A=1;
    if (with_hessian){
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dVright[i] = 0;
            model->dA[i] = 0;
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2Vright[i*(i+1)/2+j] = 0;
                model->d2A[i*(i+1)/2+j] = 0;
            }
        }
    }
    if(with_gradient) {
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dVright[i] = 0;
            model->dA[i] = 0;
        }
    }
    // save old model
    model->idMod = id;
}

void ABAO::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    model->k += 1;
    if (with_hessian){
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2Vright[i*(i+1)/2+j] = model->d2Vleft[i*(i+1)/2+j];
            }
        }
    }
    if(with_gradient||with_hessian) {
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dVright[i] =  model->dVleft[i];
        }
    }

    model->Vright = model->Vleft;

    // save old model
    model->idMod = id;
}

void AGAP::update(bool with_gradient,bool with_hessian) {
    model->k += 1;

    // save old model
    model->idMod = id;
}

void QAGAN::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    model->k += 1;
    model->Vright = 0;

    if (with_hessian){
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dVright[i] = 0;
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2Vright[i*(i+1)/2+j] = 0;
            }
        }
    }
    if(with_gradient) {
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dVright[i] = 0;
        }
    }
    // save old model
    model->idMod = id;
}

void QR::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    model->k += 1;
    model->Vright = 0;
    if (with_hessian){
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2A[i*(i+1)/2+j] = rho* model->d2A[i*(i+1)/2+j];
                model->d2Vright[i*(i+1)/2+j] = 0;
            }
        }
        for(j=0;j<=id_params;j++) {
            //i(<=id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
            model->d2A[id_params*(id_params+1)/2+j] = model->d2A[id_params*(id_params+1)/2+j] + model->dA[j];
        }
        for(i=id_params;i<model->nb_paramsMaintenance;i++) {
             //id and i(>=id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            model->d2A[i*(i+1)/2+id_params] = model->d2A[i*(i+1)/2+id_params] + model->dA[i];
        }
    }
    if(with_gradient||with_hessian) {
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dA[i] = rho *  model->dA[i];
            model->dVright[i] = 0;
        }
        model->dA[id_params] = model->dA[id_params] +  model->A;
    }
    model->A=rho*model->A;
    // save old model
    model->idMod = id;
}

void GQR::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    model->k += 1;
    model->Vright = 0;
    K++;
    if (with_hessian){
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2A[i*(i+1)/2+j] = pow(rho,f->eval(K)-f->eval(K-1))* model->d2A[i*(i+1)/2+j];
                model->d2Vright[i*(i+1)/2+j] = 0;
            }
        }
        for(j=0;j<id_params;j++) {
            //i(<=id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
            model->d2A[id_params*(id_params+1)/2+j] = model->d2A[id_params*(id_params+1)/2+j] + pow(rho,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho*model->dA[j];
        }
        model->d2A[id_params*(id_params+1)/2+id_params] = model->d2A[id_params*(id_params+1)/2+id_params] + pow(rho,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho*(2*model->dA[id_params]+(f->eval(K)-f->eval(K-1)-1)/rho*model->A);
        for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
             //id and i(>=id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            model->d2A[i*(i+1)/2+id_params] = model->d2A[i*(i+1)/2+id_params] + pow(rho,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho*model->dA[i];
        }
    }
    if(with_gradient||with_hessian) {
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dA[i] = pow(rho,f->eval(K)-f->eval(K-1)) *  model->dA[i];
            model->dVright[i] = 0;
        }
        model->dA[id_params] = model->dA[id_params] + pow(rho,f->eval(K)-f->eval(K-1)) * (f->eval(K)-f->eval(K-1))/rho * model->A;
    }
    model->A=pow(rho,f->eval(K)-f->eval(K-1))*model->A;
    // save old model
    model->idMod = id;
}

void GQR_ARA1::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    model->k += 1;
    K++;
    if (with_hessian){
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2A[i*(i+1)/2+j] = pow(rho_QR,f->eval(K)-f->eval(K-1))* model->d2A[i*(i+1)/2+j];
                model->d2Vright[i*(i+1)/2+j] = (1-rho_ARA) * model->d2Vleft[i*(i+1)/2+j] + rho_ARA * model->d2Vright[i*(i+1)/2+j];
            }
        }
        for(j=0;j<id_params;j++) {
            //i(<=id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
            model->d2A[id_params*(id_params+1)/2+j] = model->d2A[id_params*(id_params+1)/2+j] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*model->dA[j];
            model->d2Vright[(id_params+1)*(id_params+2)/2+j] = model->d2Vright[(id_params+1)*(id_params+2)/2+j] - model->dVleft[j]+model->dVright[j];
        }
        model->d2A[id_params*(id_params+1)/2+id_params] = model->d2A[id_params*(id_params+1)/2+id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*(2*model->dA[id_params]+(f->eval(K)-f->eval(K-1)-1)/rho_QR*model->A);
        model->d2Vright[(id_params+1)*(id_params+2)/2+id_params] = model->d2Vright[(id_params+1)*(id_params+2)/2+id_params] - model->dVleft[id_params]+model->dVright[id_params];
        model->d2A[(id_params+1)*(id_params+2)/2+id_params] = model->d2A[(id_params+1)*(id_params+2)/2+id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*model->dA[id_params+1];
        model->d2Vright[(id_params+1)*(id_params+2)/2+id_params+1] = model->d2Vright[(id_params+1)*(id_params+2)/2+id_params+1] - 2*model->dVleft[id_params+1]+2*model->dVright[id_params+1];
        for(i=id_params+2;i<model->nb_paramsMaintenance;i++) {
             //id and i(>=id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            model->d2A[i*(i+1)/2+id_params] = model->d2A[i*(i+1)/2+id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*model->dA[i];
            model->d2Vright[i*(i+1)/2+id_params+1] = model->d2Vright[i*(i+1)/2+id_params+1] - model->dVleft[i]+model->dVright[i];
        }
    }
    if(with_gradient||with_hessian) {
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dA[i] = pow(rho_QR,f->eval(K)-f->eval(K-1)) *  model->dA[i];
            model->dVright[i] = (1-rho_ARA) * model->dVleft[i]+ rho_ARA * model->dVright[i];
        }
        model->dA[id_params] = model->dA[id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1)) * (f->eval(K)-f->eval(K-1))/rho_QR * model->A;
        model->dVright[id_params+1] = model->dVright[id_params+1] - model->Vleft + model->Vright;
    }
    model->A=pow(rho_QR,f->eval(K)-f->eval(K-1))*model->A;
    model->Vright= (1-rho_ARA)*model->Vleft + rho_ARA*model->Vright;
    // save old model
    model->idMod = id;
}

void GQR_ARAInf::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    model->k += 1;
    K++;
    if (with_hessian){
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2A[i*(i+1)/2+j] = pow(rho_QR,f->eval(K)-f->eval(K-1))* model->d2A[i*(i+1)/2+j];
                model->d2Vright[i*(i+1)/2+j] = (1-rho_ARA) * model->d2Vleft[i*(i+1)/2+j];
            }
        }
        for(j=0;j<id_params;j++) {
            //i(<=id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
            model->d2A[id_params*(id_params+1)/2+j] = model->d2A[id_params*(id_params+1)/2+j] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*model->dA[j];
            model->d2Vright[(id_params+1)*(id_params+2)/2+j] = model->d2Vright[(id_params+1)*(id_params+2)/2+j] - model->dVleft[j];
        }
        model->d2A[id_params*(id_params+1)/2+id_params] = model->d2A[id_params*(id_params+1)/2+id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*(2*model->dA[id_params]+(f->eval(K)-f->eval(K-1)-1)/rho_QR*model->A);
        model->d2Vright[(id_params+1)*(id_params+2)/2+id_params] = model->d2Vright[(id_params+1)*(id_params+2)/2+id_params] - model->dVleft[id_params];
        model->d2A[(id_params+1)*(id_params+2)/2+id_params] = model->d2A[(id_params+1)*(id_params+2)/2+id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*model->dA[id_params+1];
        model->d2Vright[(id_params+1)*(id_params+2)/2+id_params+1] = model->d2Vright[(id_params+1)*(id_params+2)/2+id_params+1] - 2*model->dVleft[id_params+1];
        for(i=id_params+2;i<model->nb_paramsMaintenance;i++) {
             //id and i(>=id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            model->d2A[i*(i+1)/2+id_params] = model->d2A[i*(i+1)/2+id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*model->dA[i];
            model->d2Vright[i*(i+1)/2+id_params+1] = model->d2Vright[i*(i+1)/2+id_params+1] - model->dVleft[i];
        }
    }
    if(with_gradient||with_hessian) {
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dA[i] = pow(rho_QR,f->eval(K)-f->eval(K-1)) *  model->dA[i];
            model->dVright[i] = (1-rho_ARA) * model->dVleft[i];
        }
        model->dA[id_params] = model->dA[id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1)) * (f->eval(K)-f->eval(K-1))/rho_QR * model->A;
        model->dVright[id_params+1] = model->dVright[id_params+1] - model->Vleft;
    }
    model->A=pow(rho_QR,f->eval(K)-f->eval(K-1))*model->A;
    model->Vright= (1-rho_ARA)*model->Vleft;
    // save old model
    model->idMod = id;
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
    } else if(name.compare("AGAP.va.model") == 0) {
    //double rho=0.0;
    //mm=new ARAInf(rho,model);
        mm=new AGAP(model);
    } else if(name.compare("QAGAN.va.model") == 0) {
    //double rho=0.0;
    //mm=new ARAInf(rho,model);
        mm=new QAGAN(model);
    } else if(name.compare("QR.va.model") == 0) {
        NumericVector rho=NumericVector::create(params[0]);
        mm=new QR(rho,model);
    } else if(name.compare("GQR.va.model") == 0) {
        NumericVector rho=NumericVector::create(params[0]);
        if (maintenance.containsElementNamed("extra")) {
            std::string fun=maintenance["extra"];
            mm=new GQR(rho,fun,model);
        } else {
            std::string fun="identity";
            mm=new GQR(rho,fun,model);
        }
    } else if(name.compare("GQR_ARA1.va.model") == 0) {
        NumericVector rho=NumericVector::create(params[0],params[1]);
        if (maintenance.containsElementNamed("extra")) {
            std::string fun=maintenance["extra"];
            mm=new GQR_ARA1(rho,fun,model);
        } else {
            std::string fun="identity";
            mm=new GQR_ARA1(rho,fun,model);
        }
    } else if(name.compare("GQR_ARAInf.va.model") == 0) {
        NumericVector rho=NumericVector::create(params[0],params[1]);
        if (maintenance.containsElementNamed("extra")) {
            std::string fun=maintenance["extra"];
            mm=new GQR_ARAInf(rho,fun,model);
        } else {
            std::string fun="identity";
            mm=new GQR_ARAInf(rho,fun,model);
        }
    } else {
        printf("WARNING: %s is not a proper maintenance model!\n",name.c_str());
    }
	return mm;
}

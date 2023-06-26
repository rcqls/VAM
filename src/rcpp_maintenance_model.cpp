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
    int k;
    double prov;
    model->k += 1;

    int nk=model->k;
    if(nk>model->mu){
        nk=model->mu;
    }
    //printf("ARA1 k=%d,max_mem=%d, nk=%d\n",model->k,model->mu,nk);
    if (with_hessian){
        for(k=nk-1;k>0;k--){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
            }
        }
        if(nk>0) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    prov= (1-rho)*model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
                    model->d2VR_prec[i*(i+1)/2+j] = prov;
                    model->d2Vright[i*(i+1)/2+j]+=prov;
                }
            }
            for(j=0;j<=id_params;j++) {
                prov=model->dA[j]*(model->time[model->k]-model->time[model->k - 1]);
                model->d2VR_prec[id_params*(id_params+1)/2+j] -= prov;
                model->d2Vright[id_params*(id_params+1)/2+j] -= prov;
            }
            for(i=id_params;i<model->nb_paramsMaintenance;i++) {
                prov=model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
                model->d2VR_prec[i*(i+1)/2+id_params] -= prov;
                model->d2Vright[i*(i+1)/2+id_params] -= prov;
            }
        } else {
           for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2Vright[i*(i+1)/2+j]+=(1-rho)*model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
                }
            }
            for(j=0;j<=id_params;j++) {
                model->d2Vright[id_params*(id_params+1)/2+j] -= model->dA[j]*(model->time[model->k]-model->time[model->k - 1]);
            }
            for(i=id_params;i<model->nb_paramsMaintenance;i++) {
                model->d2Vright[i*(i+1)/2+id_params] -= model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
            } 
        }
    }
    if(with_gradient||with_hessian) {
        for(k=nk-1;k>0;k--){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVR_prec[k*model->nb_paramsMaintenance+i]=model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
            }
        }
        if(nk>0) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                prov=(1-rho)*model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
                model->dVR_prec[i]=prov;
                model->dVright[i]+=prov;
                // printf("ARA1: model->dVright[i=%d]=%lf prov=%lf\n",i, model->dVright[i], prov);
            }
            prov=model->A*(model->time[model->k]-model->time[model->k - 1]);
            model->dVR_prec[id_params]-= prov;
            model->dVright[id_params]-=prov;
        } else {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVright[i]+=(1-rho)*model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
                // if (i == 1 && model->dA[i] > 0) printf("ARA1: model->dVright[i=%d]=%lf dA[i]=%lf\n",i, model->dVright[i], model->dA[i]);
            }
            model->dVright[id_params]-=model->A*(model->time[model->k]-model->time[model->k - 1]);
        }
    }
    for(k=nk-1;k>0;k--) model->VR_prec[k]=model->VR_prec[k-1];
    prov=(1-rho)*model->A*(model->time[model->k]-model->time[model->k - 1]);
    if(nk>0) model->VR_prec[0]=prov;
    model->Vright+=prov;

    // save old model
    model->idMod = id;
}

void ARAInf::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    model->k += 1;


    int nk=model->k;
    if(nk>model->mu){
        nk=model->mu;
    }
    //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model->k,model->mu,nk);
    if (with_hessian){
        for(k=nk-1;k>0;k--){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=(1-rho)*model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
            }
            for(j=0;j<=id_params;j++) {
                model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+id_params*(id_params+1)/2+j]-= model->dVR_prec[(k-1)*model->nb_paramsMaintenance+j];
            }
            for(i=id_params;i<model->nb_paramsMaintenance;i++) {
                model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+id_params] -= model->dVR_prec[(k-1)*(model->nb_paramsMaintenance)+i];
            }
        }
        if(nk>0){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2VR_prec[i*(i+1)/2+j] = (1-rho)*model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
                    model->d2Vright[i*(i+1)/2+j] = (1-rho)*model->d2Vleft[i*(i+1)/2+j];
                }
            }
            for(j=0;j<=id_params;j++) {
                model->d2VR_prec[id_params*(id_params+1)/2+j] -= model->dA[j]*(model->time[model->k]-model->time[model->k - 1]);
                model->d2Vright[id_params*(id_params+1)/2+j] -= model->dVleft[j];
            }
            for(i=id_params;i<model->nb_paramsMaintenance;i++) {
                model->d2VR_prec[i*(i+1)/2+id_params] -= model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
                model->d2Vright[i*(i+1)/2+id_params] -= model->dVleft[i];
            }
        } else {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2Vright[i*(i+1)/2+j] = (1-rho)*model->d2Vleft[i*(i+1)/2+j];
                }
            }
            for(j=0;j<=id_params;j++) {
                model->d2Vright[id_params*(id_params+1)/2+j] -= model->dVleft[j];
            }
            for(i=id_params;i<model->nb_paramsMaintenance;i++) {
                model->d2Vright[i*(i+1)/2+id_params] -= model->dVleft[i];
            }
        }
    }
    if(with_gradient||with_hessian) {
        for(k=nk-1;k>0;k--){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVR_prec[k*model->nb_paramsMaintenance+i]=(1-rho)*model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
            }
            model->dVR_prec[k*model->nb_paramsMaintenance+id_params]-=model->VR_prec[k-1];
        }
        if(nk>0) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVR_prec[i]=(1-rho)*model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
                model->dVright[i]=(1-rho)*model->dVleft[i];
                // if (i == 1 && model->dVleft[i] > 0) printf("ARAInf: model->dVright[i=%d]=%lf model->dVleft[i]=%lf\n",i, model->dVright[i], model->dVleft[i]);
            }
            model->dVR_prec[id_params]-= model->A*(model->time[model->k]-model->time[model->k - 1]);
            model->dVright[id_params]-= model->Vleft;
            // if(id_params == 1 && model->Vleft != 0) printf("ARAInf (nk=%d): model->dVright[i=%d]=%lf model->Vleft=%lf\n",nk, id_params, model->dVright[id_params], model->Vleft);
        } else {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVright[i]=(1-rho)*model->dVleft[i];
                // if (i == 1 && model->dVleft[i] > 0) printf("ARAInf: model->dVright[i=%d]=%lf model->dVleft[i]=%lf\n",i, model->dVright[i], model->dVleft[i]);
            }
            model->dVright[id_params]-= model->Vleft;
            // if(id_params == 1 && model->Vleft != 0) printf("ARAInf (nk=%d): model->dVright[i=%d]=%lf model->Vleft=%lf\n",nk, id_params, model->dVright[id_params], model->Vleft);
        }
    }
    model->Vright=(1-rho)*model->Vleft;
    for(k=nk-1;k>0;k--) model->VR_prec[k]=(1-rho)*model->VR_prec[k-1];
    if(nk>0) model->VR_prec[0]=(1-rho)*model->A*(model->time[model->k]-model->time[model->k - 1]);  
    // save old model
    model->idMod = id;
}

void AGAN::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    model->k += 1;

    int nk=model->k;
    if(nk>model->mu){
        nk=model->mu;
    }
    model->A=1;
    model->Vright=0;
    for(k=0;k<nk;k++) {
        model->VR_prec[k]=0;
    }
    if (with_hessian){
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dA[i] = 0;
            model->dVright[i] = 0;
            for(k=0;k<nk;k++) {
                model->dVR_prec[k*model->nb_paramsMaintenance+i]=0;
            }
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2A[i*(i+1)/2+j] = 0;
                model->d2Vright[i*(i+1)/2+j] = 0;
                for(k=0;k<nk;k++) {
                    model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=0;
                }
            }
        }
    }
    if(with_gradient) {
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dA[i] = 0;
            model->dVright[i] = 0;
            for(k=0;k<nk;k++) {
                model->dVR_prec[k*model->nb_paramsMaintenance+i]=0;
            }
        }
    }

    // save old model
    model->idMod = id;

    //init QR and GQR type models
    for(i=0;i<model->nbPM + 1;i++) model->models->at(i)->init();
}

void ABAO::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    double prov;
    model->k += 1;

    int nk=model->k;
    if(nk>model->mu){
        nk=model->mu;
    }
    if (with_hessian){
        for(k=nk-1;k>0;k--){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
            }
        }
        if(nk>0) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    prov= model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
                    model->d2VR_prec[i*(i+1)/2+j] = prov;
                    model->d2Vright[i*(i+1)/2+j]+=prov;
                }
            }
        } else {
           for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2Vright[i*(i+1)/2+j]+=model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
                }
            }
        }
    }
    if(with_gradient||with_hessian) {
        for(k=nk-1;k>0;k--){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVR_prec[k*model->nb_paramsMaintenance+i]=model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
            }
        }
        if(nk>0) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                prov=model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
                model->dVR_prec[i]=prov;
                model->dVright[i]+=prov;
            }
        } else {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVright[i]+=model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
            }
        }
    }
    for(k=nk-1;k>0;k--) model->VR_prec[k]=model->VR_prec[k-1];
    prov=model->A*(model->time[model->k]-model->time[model->k - 1]);
    if(nk>0) model->VR_prec[0]=prov;
    model->Vright+=prov;  

    printf("ABAO Vright = %lf\n", model->Vright);  

    // save old model
    model->idMod = id;
}

void AGAP::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    model->k += 1;

    int nk=model->k;
    if(nk>model->mu){
        nk=model->mu;
    }
    if (with_hessian){
        for(k=nk-1;k>0;k--){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
            }
        }
        if(nk>0) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) model->d2VR_prec[i*(i+1)/2+j] = 0;
            }
        }
    }
    if(with_gradient||with_hessian) {
        for(k=nk-1;k>0;k--){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVR_prec[k*model->nb_paramsMaintenance+i]=model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
            }
        }
        if(nk>0) {
            for(i=0;i<model->nb_paramsMaintenance;i++) model->dVR_prec[i]=0;
        }
    }
    for(k=nk-1;k>0;k--) model->VR_prec[k]=model->VR_prec[k-1];
    if(nk>0) model->VR_prec[0]=0;

    // save old model
    model->idMod = id;
}

void QAGAN::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    model->k += 1;

    int nk=model->k;
    if(nk>model->mu){
        nk=model->mu;
    }
    model->Vright=0;
    for(k=0;k<nk;k++) {
        model->VR_prec[k]=0;
    }
    if (with_hessian){
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dVright[i] = 0;
            for(k=0;k<nk;k++) {
                model->dVR_prec[k*model->nb_paramsMaintenance+i]=0;
            }
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2Vright[i*(i+1)/2+j] = 0;
                for(k=0;k<nk;k++) {
                    model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=0;
                }
            }
        }
    }
    if(with_gradient) {
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dVright[i] = 0;
            for(k=0;k<nk;k++) {
                model->dVR_prec[k*model->nb_paramsMaintenance+i]=0;
            }
        }
    }

    // save old model
    model->idMod = id;
}

void QR::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    model->k += 1;

    int nk=model->k;
    if(nk>model->mu){
        nk=model->mu;
    }
    if (with_hessian){
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2A[i*(i+1)/2+j] = rho* model->d2A[i*(i+1)/2+j];
                model->d2Vright[i*(i+1)/2+j] = 0;
                for(k=0;k<nk;k++) {
                    model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=0;
                }
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
            for(k=0;k<nk;k++) {
                model->dVR_prec[k*model->nb_paramsMaintenance+i]=0;
            }
        }
        model->dA[id_params] = model->dA[id_params] +  model->A;
    }
    model->A=rho*model->A;
    model->Vright=0;
    for(k=0;k<nk;k++) {
        model->VR_prec[k]=0;
    }
    model->idMod = id;
}

void GQR::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    model->k += 1;
    int nk=model->k;
    if(nk>model->mu){
        nk=model->mu;
    }
    K++;
    // printf("K=%lf delta=%lf\n",K,f->eval(K)-f->eval(K-1));
    if (with_hessian){
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2A[i*(i+1)/2+j] = pow(rho,f->eval(K)-f->eval(K-1))* model->d2A[i*(i+1)/2+j];
                model->d2Vright[i*(i+1)/2+j] = 0;
                for(k=0;k<nk;k++) {
                    model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=0;
                }
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
            for(k=0;k<nk;k++) {
                model->dVR_prec[k*model->nb_paramsMaintenance+i]=0;
            }
        }
        model->dA[id_params] = model->dA[id_params] + pow(rho,f->eval(K)-f->eval(K-1)) * (f->eval(K)-f->eval(K-1))/rho * model->A;
    }
    model->A=pow(rho,f->eval(K)-f->eval(K-1))*model->A;
    model->Vright=0;
    for(k=0;k<nk;k++) {
        model->VR_prec[k]=0;
    }
    model->idMod = id;
}

void GQR_ARA1::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    double prov;
    K++;
    model->k += 1;

    int nk=model->k;
    if(nk>model->mu){
        nk=model->mu;
    }
    //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model->k,model->max_mem,nk);
    if (with_hessian){
        for(k=nk-1;k>0;k--){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
            }
        }
        if(nk>0) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    prov= (1-rho_ARA)*model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
                    model->d2VR_prec[i*(i+1)/2+j] = prov;
                    model->d2Vright[i*(i+1)/2+j]+=prov;
                }
            }
            for(j=0;j<=id_params+1;j++) {
                prov=model->dA[j]*(model->time[model->k]-model->time[model->k - 1]);
                model->d2VR_prec[(id_params+1)*(id_params+2)/2+j] -= prov;
                model->d2Vright[(id_params+1)*(id_params+2)/2+j] -= prov;
            }
            for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
                prov=model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
                model->d2VR_prec[i*(i+1)/2+id_params+1] -= prov;
                model->d2Vright[i*(i+1)/2+id_params+1] -= prov;
            }
        } else {
           for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2Vright[i*(i+1)/2+j]+=(1-rho_ARA)*model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
                }
            }
            for(j=0;j<=id_params+1;j++) {
                model->d2Vright[(id_params+1)*(id_params+2)/2+j] -= model->dA[j]*(model->time[model->k]-model->time[model->k - 1]);
            }
            for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
                model->d2Vright[i*(i+1)/2+id_params+1] -= model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
            } 
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2A[i*(i+1)/2+j] = pow(rho_QR,f->eval(K)-f->eval(K-1))* model->d2A[i*(i+1)/2+j];
            }
        }
        for(j=0;j<id_params;j++) {
            //i(<=id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
            model->d2A[id_params*(id_params+1)/2+j] = model->d2A[id_params*(id_params+1)/2+j] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*model->dA[j];
        }
        model->d2A[id_params*(id_params+1)/2+id_params] = model->d2A[id_params*(id_params+1)/2+id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*(2*model->dA[id_params]+(f->eval(K)-f->eval(K-1)-1)/rho_QR*model->A);
        for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
             //id and i(>=id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            model->d2A[i*(i+1)/2+id_params] = model->d2A[i*(i+1)/2+id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*model->dA[i];
        }
    }
    if(with_gradient||with_hessian) {
        for(k=nk-1;k>0;k--){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVR_prec[k*model->nb_paramsMaintenance+i]=model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
            }
        }
        if(nk>0) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                prov=(1-rho_ARA)*model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
                model->dVR_prec[i]=prov;
                model->dVright[i]+=prov;
            }
            prov=model->A*(model->time[model->k]-model->time[model->k - 1]);
            model->dVR_prec[id_params+1]-= prov;
            model->dVright[id_params+1]-=prov;
        } else {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVright[i]+=(1-rho_ARA)*model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
            }
            model->dVright[id_params+1]-=model->A*(model->time[model->k]-model->time[model->k - 1]);
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dA[i] = pow(rho_QR,f->eval(K)-f->eval(K-1)) *  model->dA[i];
        }
        model->dA[id_params] = model->dA[id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1)) * (f->eval(K)-f->eval(K-1))/rho_QR * model->A;

    }
    for(k=nk-1;k>0;k--) model->VR_prec[k]=model->VR_prec[k-1];
    prov=(1-rho_ARA)*model->A*(model->time[model->k]-model->time[model->k - 1]);
    if(nk>0) model->VR_prec[0]=prov;
    model->Vright+=prov;
    model->A=pow(rho_QR,f->eval(K)-f->eval(K-1))*model->A;

    // save old model
    model->idMod = id;
}

void GQR_ARAInf::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    K++;
    model->k += 1;


    int nk=model->k;
    if(nk>model->mu){
        nk=model->mu;
    }
    //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model->k,model->max_mem,nk);
    if (with_hessian){
        for(k=nk-1;k>0;k--){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=(1-rho_ARA)*model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
            }
            for(j=0;j<=id_params+1;j++) {
                model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+(id_params+1)*(id_params+2)/2+j]-= model->dVR_prec[(k-1)*model->nb_paramsMaintenance+j];
            }
            for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
                model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+id_params+1] -= model->dVR_prec[(k-1)*(model->nb_paramsMaintenance)+i];
            }
        }
        if(nk>0){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2VR_prec[i*(i+1)/2+j] = (1-rho_ARA)*model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
                    model->d2Vright[i*(i+1)/2+j] = (1-rho_ARA)*model->d2Vleft[i*(i+1)/2+j];
                }
            }
            for(j=0;j<=id_params+1;j++) {
                model->d2VR_prec[(id_params+1)*(id_params+2)/2+j] -= model->dA[j]*(model->time[model->k]-model->time[model->k - 1]);
                model->d2Vright[(id_params+1)*(id_params+2)/2+j] -= model->dVleft[j];
            }
            for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
                model->d2VR_prec[i*(i+1)/2+id_params+1] -= model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
                model->d2Vright[i*(i+1)/2+id_params+1] -= model->dVleft[i];
            }
        } else {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2Vright[i*(i+1)/2+j] = (1-rho_ARA)*model->d2Vleft[i*(i+1)/2+j];
                }
            }
            for(j=0;j<=id_params+1;j++) {
                model->d2Vright[(id_params+1)*(id_params+2)/2+j] -= model->dVleft[j];
            }
            for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
                model->d2Vright[i*(i+1)/2+id_params+1] -= model->dVleft[i];
            }
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2A[i*(i+1)/2+j] = pow(rho_QR,f->eval(K)-f->eval(K-1))* model->d2A[i*(i+1)/2+j];
            }
        }
        for(j=0;j<id_params;j++) {
            //i(<=id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
            model->d2A[id_params*(id_params+1)/2+j] = model->d2A[id_params*(id_params+1)/2+j] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*model->dA[j];
        }
        model->d2A[id_params*(id_params+1)/2+id_params] = model->d2A[id_params*(id_params+1)/2+id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*(2*model->dA[id_params]+(f->eval(K)-f->eval(K-1)-1)/rho_QR*model->A);
        for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
             //id and i(>=id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            model->d2A[i*(i+1)/2+id_params] = model->d2A[i*(i+1)/2+id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*model->dA[i];
        }
    }
    if(with_gradient||with_hessian) {
        for(k=nk-1;k>0;k--){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVR_prec[k*model->nb_paramsMaintenance+i]=(1-rho_ARA)*model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
            }
            model->dVR_prec[k*model->nb_paramsMaintenance+id_params+1]-=model->VR_prec[k-1];
        }
        if(nk>0) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVR_prec[i]=(1-rho_ARA)*model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
                model->dVright[i]=(1-rho_ARA)*model->dVleft[i];
            }
            model->dVR_prec[id_params+1]-= model->A*(model->time[model->k]-model->time[model->k - 1]);
            model->dVright[id_params+1]-= model->Vleft;
        } else {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVright[i]=(1-rho_ARA)*model->dVleft[i];
            }
            model->dVright[id_params+1]-= model->Vleft;
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dA[i] = pow(rho_QR,f->eval(K)-f->eval(K-1)) *  model->dA[i];
        }
        model->dA[id_params] = model->dA[id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1)) * (f->eval(K)-f->eval(K-1))/rho_QR * model->A;
    }
    model->Vright=(1-rho_ARA)*model->Vleft;
    for(k=nk-1;k>0;k--) model->VR_prec[k]=(1-rho_ARA)*model->VR_prec[k-1];
    if(nk>0) model->VR_prec[0]=(1-rho_ARA)*model->A*(model->time[model->k]-model->time[model->k - 1]);
    model->A=pow(rho_QR,f->eval(K)-f->eval(K-1))*model->A;

    // save old model
    model->idMod = id;
}

void ARAm::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    double prov;
    model->k += 1;


    int nk=model->k;
    if(nk>model->mu){
        nk=model->mu;
    }
    int nk2=nk;
    if(nk>m-1){
        nk2=m-1;
    }
    
    // printf("ARAm=%d, k=%d,max_mem=%d, nk=%d, nk2=%d\n",m,model->k,model->mu,nk, nk2);
    if (with_hessian){
        for(k=nk-1;k>nk2;k--){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
            }
        }

        if ((model->k>=m)&&(nk2>0)) {
            if(nk>nk2) {
                for(i=0;i<model->nb_paramsMaintenance;i++) {
                    for(j=0;j<=i;j++) {
                        model->d2Vright[i*(i+1)/2+j]-=rho*model->d2VR_prec[(nk2-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                        model->d2VR_prec[nk2*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=(1-rho)*model->d2VR_prec[(nk2-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                    }
                }
                for(j=0;j<=id_params;j++) {
                    model->d2Vright[id_params*(id_params+1)/2+j]-= model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+j];
                    model->d2VR_prec[nk2*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+id_params*(id_params+1)/2+j]-= model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+j];
                }
                for(i=id_params;i<model->nb_paramsMaintenance;i++) {
                    model->d2Vright[i*(i+1)/2+id_params] -= model->dVR_prec[(nk2-1)*(model->nb_paramsMaintenance)+i];
                    model->d2VR_prec[nk2*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+id_params] -= model->dVR_prec[(nk2-1)*(model->nb_paramsMaintenance)+i];
                }
            } else {
                for(i=0;i<model->nb_paramsMaintenance;i++) {
                    for(j=0;j<=i;j++) model->d2Vright[i*(i+1)/2+j]-=rho*model->d2VR_prec[(nk2-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
                for(j=0;j<=id_params;j++) model->d2Vright[id_params*(id_params+1)/2+j]-= model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+j];
                for(i=id_params;i<model->nb_paramsMaintenance;i++) model->d2Vright[i*(i+1)/2+id_params] -= model->dVR_prec[(nk2-1)*(model->nb_paramsMaintenance)+i];
            }
        }

        for(k=nk2-1;k>0;k--){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2Vright[i*(i+1)/2+j]-=rho*model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                    model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=(1-rho)*model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
            }
            for(j=0;j<=id_params;j++) {
                model->d2Vright[id_params*(id_params+1)/2+j]-= model->dVR_prec[(k-1)*model->nb_paramsMaintenance+j];
                model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+id_params*(id_params+1)/2+j]-= model->dVR_prec[(k-1)*model->nb_paramsMaintenance+j];
            }
            for(i=id_params;i<model->nb_paramsMaintenance;i++) {
                model->d2Vright[i*(i+1)/2+id_params] -= model->dVR_prec[(k-1)*(model->nb_paramsMaintenance)+i];
                model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+id_params] -= model->dVR_prec[(k-1)*(model->nb_paramsMaintenance)+i];
            }
        }
        if(nk>0) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    prov= (1-rho)*model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
                    model->d2VR_prec[i*(i+1)/2+j] = prov;
                    model->d2Vright[i*(i+1)/2+j]+=prov;
                }
            }
            for(j=0;j<=id_params;j++) {
                prov=model->dA[j]*(model->time[model->k]-model->time[model->k - 1]);
                model->d2VR_prec[id_params*(id_params+1)/2+j] -= prov;
                model->d2Vright[id_params*(id_params+1)/2+j] -= prov;
            }
            for(i=id_params;i<model->nb_paramsMaintenance;i++) {
                prov=model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
                model->d2VR_prec[i*(i+1)/2+id_params] -= prov;
                model->d2Vright[i*(i+1)/2+id_params] -= prov;
            }
        } else {
           for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2Vright[i*(i+1)/2+j]+=(1-rho)*model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
                }
            }
            for(j=0;j<=id_params;j++) {
                model->d2Vright[id_params*(id_params+1)/2+j] -= model->dA[j]*(model->time[model->k]-model->time[model->k - 1]);
            }
            for(i=id_params;i<model->nb_paramsMaintenance;i++) {
                model->d2Vright[i*(i+1)/2+id_params] -= model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
            } 
        }
    }
    if(with_gradient||with_hessian) {
        for(k=nk-1;k>nk2;k--){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVR_prec[k*model->nb_paramsMaintenance+i]=model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
            }
        }
        if ((model->k>=m)&&(nk2>0)) {
            if(nk>nk2) {
                for(i=0;i<model->nb_paramsMaintenance;i++) {
                    model->dVright[i]-=rho*model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+i];
                    model->dVR_prec[nk2*model->nb_paramsMaintenance+i]=(1-rho)*model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+i];
                }
                model->dVR_prec[nk2 * (model->nb_paramsMaintenance)+id_params] -= model->VR_prec[nk2-1];
                model->dVright[id_params]-=model->VR_prec[nk2-1];
            } else {
                for(i=0;i<model->nb_paramsMaintenance;i++) model->dVright[i]-=rho*model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+i];
                model->dVright[id_params]-=model->VR_prec[nk2-1];
            }
        }
        for(k=nk2-1;k>0;k--){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVright[i]-=rho*model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
                model->dVR_prec[k*model->nb_paramsMaintenance+i]=(1-rho)*model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
            }
            model->dVR_prec[k * (model->nb_paramsMaintenance)+id_params] -= model->VR_prec[k-1];
            model->dVright[id_params]-=model->VR_prec[k-1];
        }
        if(nk>0) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                prov=(1-rho)*model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
                model->dVR_prec[i]=prov;
                model->dVright[i]+=prov;
            }
            prov=model->A*(model->time[model->k]-model->time[model->k - 1]);
            model->dVR_prec[id_params]-= prov;
            model->dVright[id_params]-=prov;
        } else {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVright[i]+=(1-rho)*model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
            }
            model->dVright[id_params]-=model->A*(model->time[model->k]-model->time[model->k - 1]);
        }
    }
    //printf("Avant nk=%i, nk2=%i, VRprec[0]=%f, VRprec[1]=%f, VRprec[2]=%f\n",nk,nk2,model->VR_prec[0],model->VR_prec[1],model->VR_prec[2]); 
    for(k=nk-1;k>nk2;k--) {
        model->VR_prec[k]=model->VR_prec[k-1];
    }
    if ((model->k>=m)&&(nk2>0)) {
        if(nk>nk2) {
            model->Vright-=rho*model->VR_prec[nk2-1];
            model->VR_prec[nk2]=(1-rho)*model->VR_prec[nk2-1];
        } else model->Vright-=rho*model->VR_prec[nk2-1];
    }
    for(k=nk2-1;k>0;k--) {
        model->Vright-=rho*model->VR_prec[k-1];
        model->VR_prec[k]=(1-rho)*model->VR_prec[k-1];
    } 
    prov=(1-rho)*model->A*(model->time[model->k]-model->time[model->k - 1]);
    //printf("Vright=%f, rho=%f, A=%f, Tk=%f, Tk-1=%f\n",model->Vright,rho,model->A,model->time[model->k],model->time[model->k - 1]);
    if(nk>0) model->VR_prec[0]=prov;
    model->Vright+=prov;
    //printf("Apres Vright=%f, nk=%i, nk2=%i, VRprec[0]=%f, VRprec[1]=%f, VRprec[2]=%f\n",model->Vright,nk,nk2,model->VR_prec[0],model->VR_prec[1],model->VR_prec[2]); 

    // save old model
    model->idMod = id;
}

void GQR_ARAm::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    double prov;
    model->k += 1;
    K++;

    int nk=model->k;
    if(nk>model->mu){
        nk=model->mu;
    }
    int nk2=nk;
    if(nk>m-1){
        nk2=m-1;
    }
    //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model->k,model->max_mem,nk);
    if (with_hessian){
        for(k=nk-1;k>nk2;k--){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
            }
        }

        if ((model->k>=m)&&(nk2>0)) {
            if(nk>nk2) {
                for(i=0;i<model->nb_paramsMaintenance;i++) {
                    for(j=0;j<=i;j++) {
                        model->d2Vright[i*(i+1)/2+j]-=rho_ARA*model->d2VR_prec[(nk2-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                        model->d2VR_prec[nk2*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=(1-rho_ARA)*model->d2VR_prec[(nk2-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                    }
                }
                for(j=0;j<=id_params+1;j++) {
                    model->d2Vright[(id_params+1)*(id_params+2)/2+j]-= model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+j];
                    model->d2VR_prec[nk2*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+(id_params+1)*(id_params+2)/2+j]-= model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+j];
                }
                for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
                    model->d2Vright[i*(i+1)/2+id_params+1] -= model->dVR_prec[(nk2-1)*(model->nb_paramsMaintenance)+i];
                    model->d2VR_prec[nk2*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+id_params+1] -= model->dVR_prec[(nk2-1)*(model->nb_paramsMaintenance)+i];
                }
            } else {
                for(i=0;i<model->nb_paramsMaintenance;i++) {
                    for(j=0;j<=i;j++) model->d2Vright[i*(i+1)/2+j]-=rho_ARA*model->d2VR_prec[(nk2-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
                for(j=0;j<=id_params+1;j++) model->d2Vright[(id_params+1)*(id_params+2)/2+j]-= model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+j];
                for(i=id_params+1;i<model->nb_paramsMaintenance;i++) model->d2Vright[i*(i+1)/2+id_params+1] -= model->dVR_prec[(nk2-1)*(model->nb_paramsMaintenance)+i];
            }
        }

        for(k=nk2-1;k>0;k--){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2Vright[i*(i+1)/2+j]-=rho_ARA*model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                    model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=(1-rho_ARA)*model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
            }
            for(j=0;j<=id_params+1;j++) {
                model->d2Vright[(id_params+1)*(id_params+2)/2+j]-= model->dVR_prec[(k-1)*model->nb_paramsMaintenance+j];
                model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+(id_params+1)*(id_params+2)/2+j]-= model->dVR_prec[(k-1)*model->nb_paramsMaintenance+j];
            }
            for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
                model->d2Vright[i*(i+1)/2+id_params+1] -= model->dVR_prec[(k-1)*(model->nb_paramsMaintenance)+i];
                model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+id_params+1] -= model->dVR_prec[(k-1)*(model->nb_paramsMaintenance)+i];
            }
        }
        if(nk>0) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    prov= (1-rho_ARA)*model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
                    model->d2VR_prec[i*(i+1)/2+j] = prov;
                    model->d2Vright[i*(i+1)/2+j]+=prov;
                }
            }
            for(j=0;j<=id_params+1;j++) {
                prov=model->dA[j]*(model->time[model->k]-model->time[model->k - 1]);
                model->d2VR_prec[(id_params+1)*(id_params+2)/2+j] -= prov;
                model->d2Vright[(id_params+1)*(id_params+2)/2+j] -= prov;
            }
            for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
                prov=model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
                model->d2VR_prec[i*(i+1)/2+id_params+1] -= prov;
                model->d2Vright[i*(i+1)/2+id_params+1] -= prov;
            }
        } else {
           for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2Vright[i*(i+1)/2+j]+=(1-rho_ARA)*model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
                }
            }
            for(j=0;j<=id_params+1;j++) {
                model->d2Vright[(id_params+1)*(id_params+2)/2+j] -= model->dA[j]*(model->time[model->k]-model->time[model->k - 1]);
            }
            for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
                model->d2Vright[i*(i+1)/2+id_params+1] -= model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
            } 
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2A[i*(i+1)/2+j] = pow(rho_QR,f->eval(K)-f->eval(K-1))* model->d2A[i*(i+1)/2+j];
            }
        }
        for(j=0;j<id_params;j++) {
            //i(<=id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
            model->d2A[id_params*(id_params+1)/2+j] = model->d2A[id_params*(id_params+1)/2+j] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*model->dA[j];
        }
        model->d2A[id_params*(id_params+1)/2+id_params] = model->d2A[id_params*(id_params+1)/2+id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*(2*model->dA[id_params]+(f->eval(K)-f->eval(K-1)-1)/rho_QR*model->A);
        for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
             //id and i(>=id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            model->d2A[i*(i+1)/2+id_params] = model->d2A[i*(i+1)/2+id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*model->dA[i];
        }
    }
    if(with_gradient||with_hessian) {
        for(k=nk-1;k>nk2;k--){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVR_prec[k*model->nb_paramsMaintenance+i]=model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
            }
        }
        if ((model->k>=m)&&(nk2>0)) {
            if(nk>nk2) {
                for(i=0;i<model->nb_paramsMaintenance;i++) {
                    model->dVright[i]-=rho_ARA*model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+i];
                    model->dVR_prec[nk2*model->nb_paramsMaintenance+i]=(1-rho_ARA)*model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+i];
                }
                model->dVR_prec[nk2 * (model->nb_paramsMaintenance)+id_params+1] -= model->VR_prec[nk2-1];
                model->dVright[id_params+1]-=model->VR_prec[nk2-1];
            } else {
                for(i=0;i<model->nb_paramsMaintenance;i++) model->dVright[i]-=rho_ARA*model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+i];
                model->dVright[id_params+1]-=model->VR_prec[nk2-1];
            }
        }
        for(k=nk2-1;k>0;k--){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVright[i]-=rho_ARA*model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
                model->dVR_prec[k*model->nb_paramsMaintenance+i]=(1-rho_ARA)*model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
            }
            model->dVR_prec[k * (model->nb_paramsMaintenance)+id_params+1] -= model->VR_prec[k-1];
            model->dVright[id_params+1]-=model->VR_prec[k-1];
        }
        if(nk>0) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                prov=(1-rho_ARA)*model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
                model->dVR_prec[i]=prov;
                model->dVright[i]+=prov;
            }
            prov=model->A*(model->time[model->k]-model->time[model->k - 1]);
            model->dVR_prec[id_params+1]-= prov;
            model->dVright[id_params+1]-=prov;
        } else {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVright[i]+=(1-rho_ARA)*model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
            }
            model->dVright[id_params+1]-=model->A*(model->time[model->k]-model->time[model->k - 1]);
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dA[i] = pow(rho_QR,f->eval(K)-f->eval(K-1)) *  model->dA[i];
        }
        model->dA[id_params] = model->dA[id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1)) * (f->eval(K)-f->eval(K-1))/rho_QR * model->A;
    }
    //printf("Avant nk=%i, nk2=%i, VRprec[0]=%f, VRprec[1]=%f, VRprec[2]=%f\n",nk,nk2,model->VR_prec[0],model->VR_prec[1],model->VR_prec[2]); 
    for(k=nk-1;k>nk2;k--) {
        model->VR_prec[k]=model->VR_prec[k-1];
    }
    if ((model->k>=m)&&(nk2>0)) {
        if(nk>nk2) {
            model->Vright-=rho_ARA*model->VR_prec[nk2-1];
            model->VR_prec[nk2]=(1-rho_ARA)*model->VR_prec[nk2-1];
        } else model->Vright-=rho_ARA*model->VR_prec[nk2-1];
    }
    for(k=nk2-1;k>0;k--) {
        model->Vright-=rho_ARA*model->VR_prec[k-1];
        model->VR_prec[k]=(1-rho_ARA)*model->VR_prec[k-1];
    } 
    prov=(1-rho_ARA)*model->A*(model->time[model->k]-model->time[model->k - 1]);
    //printf("Vright=%f, rho=%f, A=%f, Tk=%f, Tk-1=%f\n",model->Vright,rho,model->A,model->time[model->k],model->time[model->k - 1]);
    if(nk>0) model->VR_prec[0]=prov;
    model->Vright+=prov;
    //printf("Apres Vright=%f, nk=%i, nk2=%i, VRprec[0]=%f, VRprec[1]=%f, VRprec[2]=%f\n",model->Vright,nk,nk2,model->VR_prec[0],model->VR_prec[1],model->VR_prec[2]); 
    model->A=pow(rho_QR,f->eval(K)-f->eval(K-1))*model->A;

    // save old model
    model->idMod = id;
}

MaintenanceModel* newMaintenanceModel(List maintenance,VamModel* model) {
	std::string name=maintenance["name"];
	NumericVector params=maintenance["params"];
	MaintenanceModel*  mm=NULL;
    int mem;

	if(name.compare("ARA1.va.model") == 0) {
        NumericVector rho=NumericVector::create(0.5);
        if (params.size()==0){
            printf("WARNING: ARA1 model needs a parameter vector of length 1 ! It has been fixed to 0.5. \n");
        } else if(params.size()!=1){
            printf("WARNING: ARA1 model needs a parameter vector of length 1 !\n");
            rho[0]=params[0];
        } else {
            rho[0]=params[0];
        }
        if (maintenance.containsElementNamed("extra")) {
            std::string fun=maintenance["extra"];
            std::cout<<"WARNING: no function argument like "<< fun << " can be considered for ARA1 model !\n";
        }
        if (maintenance.containsElementNamed("m")) {
            mem=maintenance["m"];
            printf("WARNING: no memory argument like %d can be considered in ARA1 model : ARAm model must be used !\n",mem);
        }
		mm=new ARA1(rho,model);

	} else if(name.compare("ARAInf.va.model") == 0) {
        NumericVector rho=NumericVector::create(0.5);
        if (params.size()==0){
            printf("WARNING: ARAInf model needs a parameter vector of length 1 ! It has been fixed to 0.5. \n");
        } else if(params.size()!=1){
            printf("WARNING: ARAInf model needs a parameter vector of length 1 !\n");
            rho[0]=params[0];
        } else {
            rho[0]=params[0];
        }
        if (maintenance.containsElementNamed("extra")) {
            std::string fun=maintenance["extra"];
            std::cout<<"WARNING: no function argument like "<< fun << " can be considered for ARAInf model !\n";
        }
        if (maintenance.containsElementNamed("m")) {
            mem=maintenance["m"];
            printf("WARNING: no memory argument like %d can be considered in ARAInf model : ARAm model must be used !\n",mem);
        }
		mm=new ARAInf(rho,model);

	} else if(name.compare("AGAN.va.model") == 0) {
        if (params.size()!=0){
            printf("WARNING: AGAN model needs no parameter vector ! \n");
            NumericVector rho=NumericVector::create(0.5);
        }
        if (maintenance.containsElementNamed("extra")) {
            std::string fun=maintenance["extra"];
            std::cout<<"WARNING: no function argument like "<< fun << " can be considered for AGAN model !\n";
        }
        if (maintenance.containsElementNamed("m")) {
            mem=maintenance["m"];
            printf("WARNING: no memory argument like %d can be considered in AGAN model !\n",mem);
        }
        mm=new AGAN(model);

    } else if(name.compare("ABAO.va.model") == 0) {
        if (params.size()!=0){
            printf("WARNING: ABAO model needs no parameter vector ! \n");
            NumericVector rho=NumericVector::create(0.5);
        }
        if (maintenance.containsElementNamed("extra")) {
            std::string fun=maintenance["extra"];
            std::cout<<"WARNING: no function argument like "<< fun << " can be considered for ABAO model !\n";
        }
        if (maintenance.containsElementNamed("m")) {
            mem=maintenance["m"];
            printf("WARNING: no memory argument like %d can be considered in ABAO model !\n",mem);
        }
        mm=new ABAO(model);

    } else if(name.compare("AGAP.va.model") == 0) {
        if (params.size()!=0){
            printf("WARNING: AGAP model needs no parameter vector ! \n");
            NumericVector rho=NumericVector::create(0.5);
        }
        if (maintenance.containsElementNamed("extra")) {
            std::string fun=maintenance["extra"];
            std::cout<<"WARNING: no function argument like "<< fun << " can be considered for AGAP model !\n";
        }
        if (maintenance.containsElementNamed("m")) {
            mem=maintenance["m"];
            printf("WARNING: no memory argument like %d can be considered in AGAP model !\n",mem);
        }
        mm=new AGAP(model);

    } else if(name.compare("QAGAN.va.model") == 0) {
        if (params.size()!=0){
            printf("WARNING: QAGAN model needs no parameter vector ! \n");
            NumericVector rho=NumericVector::create(0.5);
        }
        if (maintenance.containsElementNamed("extra")) {
            std::string fun=maintenance["extra"];
            std::cout<<"WARNING: no function argument like "<< fun << " can be considered for QAGAN model !\n";
        }
        if (maintenance.containsElementNamed("m")) {
            mem=maintenance["m"];
            printf("WARNING: no memory argument like %d can be considered in QAGAN model !\n",mem);
        }
        mm=new QAGAN(model);

    } else if(name.compare("QR.va.model") == 0) {
        NumericVector rho=NumericVector::create(0.5);
        if (params.size()==0){
            printf("WARNING: QR model needs a parameter vector of length 1 ! It has been fixed to 0.5. \n");
        } else if(params.size()!=1){
            printf("WARNING: QR model needs a parameter vector of length 1 !\n");
            rho[0]=params[0];
        } else {
            rho[0]=params[0];
        }
        if (maintenance.containsElementNamed("extra")) {
            std::string fun=maintenance["extra"];
            std::cout<<"WARNING: no function argument like "<< fun << " can be considered for QR model : GQR model must be used !\n";
        }
        if (maintenance.containsElementNamed("m")) {
            mem=maintenance["m"];
            printf("WARNING: no memory argument like %d can be considered in QR model !\n",mem);
        }
        mm=new QR(rho,model);

    } else if(name.compare("GQR.va.model") == 0) {
        NumericVector rho=NumericVector::create(0.5);
        if (params.size()==0){
            printf("WARNING: GQR model needs a parameter vector of length 1 ! It has been fixed to 0.5. \n");
        } else if(params.size()!=1){
            printf("WARNING: GQR model needs a parameter vector of length 1 !\n");
            rho[0]=params[0];
        } else {
            rho[0]=params[0];
        }
        if (maintenance.containsElementNamed("m")) {
            mem=maintenance["m"];
            printf("WARNING: no memory argument like %d can be considered in GQR model !\n",mem);
        }
        if (maintenance.containsElementNamed("extra")) {
            std::string fun=maintenance["extra"];
            mm=new GQR(rho,fun,model);
        } else {
            std::string fun="identity";
            mm=new GQR(rho,fun,model);
        }

    } else if(name.compare("GQR_ARA1.va.model") == 0) {
        NumericVector rho=NumericVector::create(0.5,0.5);
        if (params.size()==0){
            printf("WARNING: GQR_ARA1 model needs a parameter vector of length 2 ! It has been fixed to c(0.5,0.5). \n");
        } else if (params.size()==1){
            printf("WARNING: GQR_ARA1 model needs a parameter vector of length 2 ! It has been fixed to c(%f,0.5). \n",params[0]);
            rho[1]=params[0];
        } else if(params.size()!=2){
            printf("WARNING: GQR_ARA1 model needs a parameter vector of length 2 !\n");
            rho[0]=params[0];
            rho[1]=params[1];
        } else {
            rho[0]=params[0];
            rho[1]=params[1];
        }
        if (maintenance.containsElementNamed("m")) {
            mem=maintenance["m"];
            printf("WARNING: no memory argument like %d can be considered in GQR_ARA1 model : GQR_ARAm model must be used !\n",mem);
        }
        if (maintenance.containsElementNamed("extra")) {
            std::string fun=maintenance["extra"];
            mm=new GQR_ARA1(rho,fun,model);
        } else {
            std::string fun="identity";
            mm=new GQR_ARA1(rho,fun,model);
        }

    } else if(name.compare("GQR_ARAInf.va.model") == 0) {
        NumericVector rho=NumericVector::create(0.5,0.5);
        if (params.size()==0){
            printf("WARNING: GQR_ARAInf model needs a parameter vector of length 2 ! It has been fixed to c(0.5,0.5). \n");
        } else if (params.size()==1){
            printf("WARNING: GQR_ARAInf model needs a parameter vector of length 2 ! It has been fixed to c(%f,0.5). \n",params[0]);
            rho[1]=params[0];
        } else if(params.size()!=2){
            printf("WARNING: GQR_ARAInf model needs a parameter vector of length 2 !\n");
            rho[0]=params[0];
            rho[1]=params[1];
        } else {
            rho[0]=params[0];
            rho[1]=params[1];
        }
        if (maintenance.containsElementNamed("m")) {
            mem=maintenance["m"];
            printf("WARNING: no memory argument like %d can be considered in GQR_ARAInf model : GQR_ARAm model must be used !\n",mem);
        }
        if (maintenance.containsElementNamed("extra")) {
            std::string fun=maintenance["extra"];
            mm=new GQR_ARAInf(rho,fun,model);
        } else {
            std::string fun="identity";
            mm=new GQR_ARAInf(rho,fun,model);
        }

    } else if(name.compare("ARAm.va.model") == 0) {
        NumericVector rho=NumericVector::create(0.5);
        if (params.size()==0){
            printf("WARNING: ARAm model needs a parameter vector of length 1 ! It has been fixed to 0.5. \n");
        } else if(params.size()!=1){
            printf("WARNING: ARAm model needs a parameter vector of length 1 !\n");
            rho[0]=params[0];;
        } else {
            rho[0]=params[0];
        }
        if (maintenance.containsElementNamed("m")) {
            mm=new ARAm(rho,maintenance["m"],model);
        } else if(maintenance.containsElementNamed("extra")){
            std::string fun=maintenance["extra"];
            std::cout<<"WARNING:"<< fun << "is not a proper memory for ARAm model, it has been replaced by 1!\n";
            mm=new ARAm(rho,1,model);
        } else {
            mm=new ARAm(rho,1,model);
        }

    } else if(name.compare("GQR_ARAm.va.model") == 0) {
        NumericVector rho=NumericVector::create(0.5,0.5);
        if (params.size()==0){
            printf("WARNING: GQR_ARAm model needs a parameter vector of length 2 ! It has been fixed to c(0.5,0.5). \n");
        } else if (params.size()==1){
            printf("WARNING: GQR_ARAm model needs a parameter vector of length 2 ! It has been fixed to c(%f,0.5). \n",params[0]);
            rho[1]=params[0];
        } else if(params.size()!=2){
            printf("WARNING: GQR_ARAm model needs a parameter vector of length 2 !\n");
            rho[0]=params[0];
            rho[1]=params[1];
        } else {
            rho[0]=params[0];
            rho[1]=params[1];
        }
        if (maintenance.containsElementNamed("m")) {
            if (maintenance.containsElementNamed("extra")) {
                std::string fun=maintenance["extra"];
                mm=new GQR_ARAm(rho,fun,maintenance["m"],model);
            } else {
                std::string fun="identity";
                mm=new GQR_ARAm(rho,fun,maintenance["m"],model);
            }
        } else {
            if (maintenance.containsElementNamed("extra")) {
                std::string fun=maintenance["extra"];
                mm=new GQR_ARAm(rho,fun,1,model);
            } else {
                std::string fun="identity";
                mm=new GQR_ARAm(rho,fun,1,model);
            }
        }

    } else {
        stop("ERROR: %s is not a proper maintenance model!\n",name.c_str());
    }
	return mm;
}

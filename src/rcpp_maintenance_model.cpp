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

void MaintenanceModel::update_Vright(bool with_gradient,bool with_hessian){
    int i;
    int j;
    int k;
    int nk=model->k;
    if(model->k>model->max_mem){
        nk=model->max_mem;
    }
    model->Vright =model->C;

    if (with_hessian) {
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dVright[i] = model->dC[i];
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2Vright[i*(i+1)/2+j]=model->d2C[i*(i+1)/2+j];
            }
        }
        for(k=0;k<nk;k++){
            model->Vright += model->B[k] *(model->time[model->k - k] - model->time[model->k - k-1]);
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVright[i] += model->dB[k*model->nb_paramsMaintenance+i] *(model->time[model->k - k] - model->time[model->k - k-1]);
                for(j=0;j<=i;j++) {
                    model->d2Vright[i*(i+1)/2+j] += model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j] *(model->time[model->k - k] - model->time[model->k - k-1]);
                }
            }
        }
    } else if (with_gradient) {
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dVright[i] = model->dC[i];
        }
        for(k=0;k<nk;k++){
            model->Vright += model->B[k] *(model->time[model->k - k] - model->time[model->k - k-1]);
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dVright[i] += model->dB[k*model->nb_paramsMaintenance+i] *(model->time[model->k - k] - model->time[model->k - k-1]);
            }
        }
    } else {
        for(k=0;k<nk;k++){
            model->Vright += model->B[k] *(model->time[model->k - k] - model->time[model->k - k-1]);
        }
    }
}

void ARA1::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    model->k += 1;


    int nk=model->k;
    if(model->k>model->max_mem){
        nk=model->max_mem;
    }
    //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model->k,model->max_mem,nk);
    if (with_hessian){
        if(model->k>model->max_mem){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2C[i*(i+1)/2+j]=(model->d2B[(model->max_mem - 1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->d2C[i*(i+1)/2+j]);
                }
            }
        }
        for(k=nk-1;k>0;k--) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=model->d2B[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
            }
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                model->d2B[i*(i+1)/2+j] = (1-rho)*model->d2A[i*(i+1)/2+j];
            }
        }
        for(j=0;j<=id_params;j++) {
            model->d2B[id_params*(id_params+1)/2+j] -= model->dA[j];
        }
        for(i=id_params;i<model->nb_paramsMaintenance;i++) {
            model->d2B[i*(i+1)/2+id_params] -= model->dA[i];
        }
    }  
    if(with_gradient||with_hessian) {
        if(model->k>model->max_mem){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dC[i]=(model->dB[(model->max_mem - 1)*model->nb_paramsMaintenance+i]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->dC[i]);
            }
        } 
        for(k=nk-1;k>0;k--) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dB[k*model->nb_paramsMaintenance+i]=model->dB[(k-1)*model->nb_paramsMaintenance+i];
            }
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dB[i]=(1-rho)*model->dA[i];
        }
        model->dB[id_params]-= model->A;
    }
    if(model->k>model->max_mem){
        model->C=((model->B[model->max_mem - 1])*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->C);
    }
    for(k=nk-1;k>0;k--) {
        model->B[k]=model->B[k-1];
    }
    model->B[0]=(1-rho)*model->A;
    update_Vright(with_gradient,with_hessian);

    // save old model
    model->idMod = id;
}

void ARAInf::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    model->k += 1;


    int nk=model->k;
    if(model->k>model->max_mem){
        nk=model->max_mem;
    }
    //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model->k,model->max_mem,nk);
    if (with_hessian){
        if(model->k>model->max_mem){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2C[i*(i+1)/2+j]=(1-rho)*(model->d2B[(model->max_mem - 1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->d2C[i*(i+1)/2+j]);
                }
            }
            for(j=0;j<=id_params;j++) {
                model->d2C[id_params*(id_params+1)/2+j]-= model->dB[(model->max_mem - 1)*(model->nb_paramsMaintenance)+j]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->dC[j];
            }
            for(i=id_params;i<model->nb_paramsMaintenance;i++) {
                model->d2C[i*(i+1)/2+id_params] -= model->dB[(model->max_mem - 1)*(model->nb_paramsMaintenance)+i]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->dC[i];
            }
        }
        for(k=nk-1;k>0;k--) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=(1-rho)*model->d2B[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
            }
            for(j=0;j<=id_params;j++) {
                model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+id_params*(id_params+1)/2+j]-= model->dB[(k-1)*model->nb_paramsMaintenance+j];
            }
            for(i=id_params;i<model->nb_paramsMaintenance;i++) {
                model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+id_params] -= model->dB[(k-1)*(model->nb_paramsMaintenance)+i];
            }
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                model->d2B[i*(i+1)/2+j] = (1-rho)*model->d2A[i*(i+1)/2+j];
            }
        }
        for(j=0;j<=id_params;j++) {
            model->d2B[id_params*(id_params+1)/2+j] -= model->dA[j];
        }
        for(i=id_params;i<model->nb_paramsMaintenance;i++) {
            model->d2B[i*(i+1)/2+id_params] -= model->dA[i];
        }
    }  
    if(with_gradient||with_hessian) {
        if(model->k>model->max_mem){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dC[i]=(1-rho)*(model->dB[(model->max_mem - 1)*model->nb_paramsMaintenance+i]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->dC[i]);
            }
            model->dC[id_params]-= (model->B[model->max_mem - 1])*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1]) +model->C ;
        } 
        for(k=nk-1;k>0;k--) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dB[k*model->nb_paramsMaintenance+i]=(1-rho)*model->dB[(k-1)*model->nb_paramsMaintenance+i];
            }
            model->dB[k * (model->nb_paramsMaintenance)+id_params] -= model->B[k-1];
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dB[i]=(1-rho)*model->dA[i];
        }
        model->dB[id_params]-= model->A;
    }
    if(model->k>model->max_mem){
        model->C=(1-rho)*((model->B[model->max_mem - 1])*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->C);
    }
    for(k=nk-1;k>0;k--) {
        model->B[k]=(1-rho)*model->B[k-1];
    }
    model->B[0]=(1-rho)*model->A;
    update_Vright(with_gradient,with_hessian);

    // save old model
    model->idMod = id;
}

void AGAN::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    model->k += 1;

    int nk=model->k;
    if(model->k>model->max_mem){
        nk=model->max_mem;
    }
    model->A=1;
    for(k=0;k<nk;k++) {
        model->B[k]=0;
    }
    model->C=0;
    if (with_hessian){
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dA[i] = 0;
            model->dC[i] = 0;
            for(k=0;k<nk;k++) {
                model->dB[k*model->nb_paramsMaintenance+i]=0;
            }
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2A[i*(i+1)/2+j] = 0;
                model->d2C[i*(i+1)/2+j] = 0;
                for(k=0;k<nk;k++) {
                    model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=0;
                }
            }
        }
    }
    if(with_gradient) {
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dA[i] = 0;
            model->dC[i] = 0;
            for(k=0;k<nk;k++) {
                model->dB[k*model->nb_paramsMaintenance+i]=0;
            }
        }
    }
    update_Vright(with_gradient,with_hessian);

    // save old model
    model->idMod = id;
}

void ABAO::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    model->k += 1;

    int nk=model->k;
    if(model->k>model->max_mem){
        nk=model->max_mem;
        model->C+=model->B[model->max_mem - 1]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1]);
    }
    model->B[0]=model->A;
    for(k=nk-1;k>0;k--) {
        model->B[k]=model->B[k-1];
    }

    if (with_hessian){
        if(model->k>model->max_mem){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dC[i]+=model->dB[(model->max_mem - 1)*model->nb_paramsMaintenance+i]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1]);
                for(j=0;j<=i;j++) {
                    model->d2C[i*(i+1)/2+j]+=model->d2B[(model->max_mem - 1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1]);
                }
            }
        }
        for(k=nk-1;k>0;k--) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dB[k*model->nb_paramsMaintenance+i]=model->dB[(k-1)*model->nb_paramsMaintenance+i];
                for(j=0;j<=i;j++) {
                    model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=model->d2B[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
            }
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dB[i]=model->dA[i];
            for(j=0;j<=i;j++) {
                model->d2B[i*(i+1)/2+j] = model->d2A[i*(i+1)/2+j];
            }
        }
    } else if(with_gradient) {
        if(model->k>model->max_mem){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dC[i]+=model->dB[(model->max_mem - 1)*model->nb_paramsMaintenance+i]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1]);
            }
        }
        for(k=nk-1;k>0;k--) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dB[k*model->nb_paramsMaintenance+i]=model->dB[(k-1)*model->nb_paramsMaintenance+i];
            }
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dB[i]=model->dA[i];
        }
    }

    update_Vright(with_gradient,with_hessian);

    // save old model
    model->idMod = id;
}

void AGAP::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    model->k += 1;

    int nk=model->k;
    if(model->k>model->max_mem){
        nk=model->max_mem;
        model->C+=model->B[model->max_mem - 1]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1]);
    }
    model->B[0]=0;
    for(k=nk-1;k>0;k--) {
        model->B[k]=model->B[k-1];
    }

    if (with_hessian){
        if(model->k>model->max_mem){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dC[i]+=model->dB[(model->max_mem - 1)*model->nb_paramsMaintenance+i]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1]);
                for(j=0;j<=i;j++) {
                    model->d2C[i*(i+1)/2+j]+=model->d2B[(model->max_mem - 1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1]);
                }
            }
        }
        for(k=nk-1;k>0;k--) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dB[k*model->nb_paramsMaintenance+i]=model->dB[(k-1)*model->nb_paramsMaintenance+i];
                for(j=0;j<=i;j++) {
                    model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=model->d2B[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
            }
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dB[i]=0;
            for(j=0;j<=i;j++) {
                model->d2B[i*(i+1)/2+j] = 0;
            }
        }
    } else if(with_gradient) {
        if(model->k>model->max_mem){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dC[i]+=model->dB[(model->max_mem - 1)*model->nb_paramsMaintenance+i]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1]);
            }
        }
        for(k=nk-1;k>0;k--) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dB[k*model->nb_paramsMaintenance+i]=model->dB[(k-1)*model->nb_paramsMaintenance+i];
            }
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dB[i]=0;
        }
    }

    update_Vright(with_gradient,with_hessian);

    // save old model
    model->idMod = id;
}

void QAGAN::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    model->k += 1;

    int nk=model->k;
    if(model->k>model->max_mem){
        nk=model->max_mem;
    }
    model->C=0;
    for(k=0;k<nk;k++) {
        model->B[k]=0;
    }
    if (with_hessian){
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dC[i] = 0;
            for(k=0;k<nk;k++) {
                model->dB[k*model->nb_paramsMaintenance+i]=0;
            }
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2C[i*(i+1)/2+j] = 0;
                for(k=0;k<nk;k++) {
                    model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=0;
                }
            }
        }
    }
    if(with_gradient) {
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dC[i] = 0;
            for(k=0;k<nk;k++) {
                model->dB[k*model->nb_paramsMaintenance+i]=0;
            }
        }
    }
    update_Vright(with_gradient,with_hessian);

    // save old model
    model->idMod = id;
}

void QR::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    model->k += 1;
    int nk=model->k;
    if(model->k>model->max_mem){
        nk=model->max_mem;
    }
    if (with_hessian){
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2A[i*(i+1)/2+j] = rho* model->d2A[i*(i+1)/2+j];
                model->d2C[i*(i+1)/2+j] = 0;
                for(k=0;k<nk;k++) {
                    model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=0;
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
            model->dC[i] = 0;
            for(k=0;k<nk;k++) {
                model->dB[k*model->nb_paramsMaintenance+i]=0;
            }
        }
        model->dA[id_params] = model->dA[id_params] +  model->A;
    }
    model->A=rho*model->A;
    for(k=0;k<nk;k++) {
        model->B[k]=0;
    }
    model->C=0;
    update_Vright(with_gradient,with_hessian);
    // save old model
    model->idMod = id;
}

void GQR::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    model->k += 1;
    int nk=model->k;
    if(model->k>model->max_mem){
        nk=model->max_mem;
    }
    K++;
    if (with_hessian){
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2A[i*(i+1)/2+j] = pow(rho,f->eval(K)-f->eval(K-1))* model->d2A[i*(i+1)/2+j];
                model->d2C[i*(i+1)/2+j] = 0;
                for(k=0;k<nk;k++) {
                    model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=0;
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
            model->dC[i] = 0;
            for(k=0;k<nk;k++) {
                model->dB[k*model->nb_paramsMaintenance+i]=0;
            }
        }
        model->dA[id_params] = model->dA[id_params] + pow(rho,f->eval(K)-f->eval(K-1)) * (f->eval(K)-f->eval(K-1))/rho * model->A;
    }
    model->A=pow(rho,f->eval(K)-f->eval(K-1))*model->A;
    for(k=0;k<nk;k++) {
        model->B[k]=0;
    }
    model->C=0;
    update_Vright(with_gradient,with_hessian);
    // save old model
    model->idMod = id;
}

void GQR_ARA1::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    K++;
    model->k += 1;


    int nk=model->k;
    if(model->k>model->max_mem){
        nk=model->max_mem;
    }
    //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model->k,model->max_mem,nk);
    if (with_hessian){
        if(model->k>model->max_mem){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2C[i*(i+1)/2+j]=(model->d2B[(model->max_mem - 1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->d2C[i*(i+1)/2+j]);
                }
            }
        }
        for(k=nk-1;k>0;k--) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=model->d2B[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
            }
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                model->d2B[i*(i+1)/2+j] = (1-rho_ARA)*model->d2A[i*(i+1)/2+j];
            }
        }
        for(j=0;j<=id_params+1;j++) {
            model->d2B[(id_params+1)*(id_params+2)/2+j] -= model->dA[j];
        }
        for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
            model->d2B[i*(i+1)/2+id_params+1] -= model->dA[i];
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
        if(model->k>model->max_mem){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dC[i]=(model->dB[(model->max_mem - 1)*model->nb_paramsMaintenance+i]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->dC[i]);
            }
        } 
        for(k=nk-1;k>0;k--) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dB[k*model->nb_paramsMaintenance+i]=model->dB[(k-1)*model->nb_paramsMaintenance+i];
            }
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dB[i]=(1-rho_ARA)*model->dA[i];
        }
        model->dB[id_params+1]-= model->A;
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dA[i] = pow(rho_QR,f->eval(K)-f->eval(K-1)) *  model->dA[i];
        }
        model->dA[id_params] = model->dA[id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1)) * (f->eval(K)-f->eval(K-1))/rho_QR * model->A;

    }
    if(model->k>model->max_mem){
        model->C=((model->B[model->max_mem - 1])*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->C);
    }
    for(k=nk-1;k>0;k--) {
        model->B[k]=model->B[k-1];
    }
    model->B[0]=(1-rho_ARA)*model->A;
    model->A=pow(rho_QR,f->eval(K)-f->eval(K-1))*model->A;
    update_Vright(with_gradient,with_hessian);

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
    if(model->k>model->max_mem){
        nk=model->max_mem;
    }
    //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model->k,model->max_mem,nk);
    if (with_hessian){
        if(model->k>model->max_mem){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2C[i*(i+1)/2+j]=(1-rho_ARA)*(model->d2B[(model->max_mem - 1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->d2C[i*(i+1)/2+j]);
                }
            }
            for(j=0;j<=id_params+1;j++) {
                model->d2C[(id_params+1)*(id_params+2)/2+j]-= model->dB[(model->max_mem - 1)*(model->nb_paramsMaintenance)+j]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->dC[j];
            }
            for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
                model->d2C[i*(i+1)/2+id_params+1] -= model->dB[(model->max_mem - 1)*(model->nb_paramsMaintenance)+i]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->dC[i];
            }
        }
        for(k=nk-1;k>0;k--) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=(1-rho_ARA)*model->d2B[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
            }
            for(j=0;j<=id_params+1;j++) {
                model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+(id_params+1)*(id_params+2)/2+j]-= model->dB[(k-1)*model->nb_paramsMaintenance+j];
            }
            for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
                model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+id_params+1] -= model->dB[(k-1)*(model->nb_paramsMaintenance)+i];
            }
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                model->d2B[i*(i+1)/2+j] = (1-rho_ARA)*model->d2A[i*(i+1)/2+j];
            }
        }
        for(j=0;j<=id_params+1;j++) {
            model->d2B[(id_params+1)*(id_params+2)/2+j] -= model->dA[j];
        }
        for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
            model->d2B[i*(i+1)/2+id_params+1] -= model->dA[i];
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
        if(model->k>model->max_mem){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dC[i]=(1-rho_ARA)*(model->dB[(model->max_mem - 1)*model->nb_paramsMaintenance+i]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->dC[i]);
            }
            model->dC[id_params+1]-= (model->B[model->max_mem - 1])*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1]) +model->C ;
        } 
        for(k=nk-1;k>0;k--) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dB[k*model->nb_paramsMaintenance+i]=(1-rho_ARA)*model->dB[(k-1)*model->nb_paramsMaintenance+i];
            }
            model->dB[k * (model->nb_paramsMaintenance)+id_params+1] -= model->B[k-1];
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dB[i]=(1-rho_ARA)*model->dA[i];
        }
        model->dB[id_params+1]-= model->A;
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dA[i] = pow(rho_QR,f->eval(K)-f->eval(K-1)) *  model->dA[i];
        }
        model->dA[id_params] = model->dA[id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1)) * (f->eval(K)-f->eval(K-1))/rho_QR * model->A;
    }
    if(model->k>model->max_mem){
        model->C=(1-rho_ARA)*((model->B[model->max_mem - 1])*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->C);
    }
    for(k=nk-1;k>0;k--) {
        model->B[k]=(1-rho_ARA)*model->B[k-1];
    }
    model->B[0]=(1-rho_ARA)*model->A;
    model->A=pow(rho_QR,f->eval(K)-f->eval(K-1))*model->A;
    update_Vright(with_gradient,with_hessian);

    // save old model
    model->idMod = id;
}

void ARAm::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    model->k += 1;


    int nk=model->k;
    if(model->k>model->max_mem){
        nk=model->max_mem;
    }
    int nk2=nk;
    if(nk>m){
        nk2=m;
    }
    //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model->k,model->max_mem,nk);
    if (with_hessian){
        if(model->k>model->max_mem){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2C[i*(i+1)/2+j]=(model->d2B[(model->max_mem - 1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->d2C[i*(i+1)/2+j]);
                }
            }
        }
        for(k=nk-1;k>nk2-1;k--) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=model->d2B[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
            }
        }
        for(k=nk2-1;k>0;k--) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=(1-rho)*model->d2B[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
            }
            for(j=0;j<=id_params;j++) {
                model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+id_params*(id_params+1)/2+j]-= model->dB[(k-1)*model->nb_paramsMaintenance+j];
            }
            for(i=id_params;i<model->nb_paramsMaintenance;i++) {
                model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+id_params] -= model->dB[(k-1)*(model->nb_paramsMaintenance)+i];
            }
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                model->d2B[i*(i+1)/2+j] = (1-rho)*model->d2A[i*(i+1)/2+j];
            }
        }
        for(j=0;j<=id_params;j++) {
            model->d2B[id_params*(id_params+1)/2+j] -= model->dA[j];
        }
        for(i=id_params;i<model->nb_paramsMaintenance;i++) {
            model->d2B[i*(i+1)/2+id_params] -= model->dA[i];
        }
    }  
    if(with_gradient||with_hessian) {
        if(model->k>model->max_mem){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dC[i]=(model->dB[(model->max_mem - 1)*model->nb_paramsMaintenance+i]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->dC[i]);
            }
        }
        for(k=nk-1;k>nk2-1;k--) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dB[k*model->nb_paramsMaintenance+i]=model->dB[(k-1)*model->nb_paramsMaintenance+i];
            }
        } 
        for(k=nk2-1;k>0;k--) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dB[k*model->nb_paramsMaintenance+i]=(1-rho)*model->dB[(k-1)*model->nb_paramsMaintenance+i];
            }
            model->dB[k * (model->nb_paramsMaintenance)+id_params] -= model->B[k-1];
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dB[i]=(1-rho)*model->dA[i];
        }
        model->dB[id_params]-= model->A;
    }
    if(model->k>model->max_mem){
        model->C=((model->B[model->max_mem - 1])*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->C);
    }
    for(k=nk-1;k>nk2-1;k--) {
        model->B[k]=model->B[k-1];
    }
    for(k=nk2-1;k>0;k--) {
        model->B[k]=(1-rho)*model->B[k-1];
    }
    model->B[0]=(1-rho)*model->A;
    update_Vright(with_gradient,with_hessian);

    // save old model
    model->idMod = id;
}

void GQR_ARAm::update(bool with_gradient,bool with_hessian) {
    int i;
    int j;
    int k;
    model->k += 1;
    K++;

    int nk=model->k;
    if(model->k>model->max_mem){
        nk=model->max_mem;
    }
    int nk2=nk;
    if(nk>m){
        nk2=m;
    }
    //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model->k,model->max_mem,nk);
    if (with_hessian){
        if(model->k>model->max_mem){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2C[i*(i+1)/2+j]=(model->d2B[(model->max_mem - 1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->d2C[i*(i+1)/2+j]);
                }
            }
        }
        for(k=nk-1;k>nk2-1;k--) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=model->d2B[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
            }
        }
        for(k=nk2-1;k>0;k--) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                for(j=0;j<=i;j++) {
                    model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=(1-rho_ARA)*model->d2B[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
                }
            }
            for(j=0;j<=id_params+1;j++) {
                model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+(id_params+1)*(id_params+2)/2+j]-= model->dB[(k-1)*model->nb_paramsMaintenance+j];
            }
            for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
                model->d2B[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+id_params+1] -= model->dB[(k-1)*(model->nb_paramsMaintenance)+i];
            }
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            for(j=0;j<=i;j++) {
                model->d2B[i*(i+1)/2+j] = (1-rho_ARA)*model->d2A[i*(i+1)/2+j];
            }
        }
        for(j=0;j<=id_params+1;j++) {
            model->d2B[(id_params+1)*(id_params+2)/2+j] -= model->dA[j];
        }
        for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
            model->d2B[i*(i+1)/2+id_params+1] -= model->dA[i];
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
        if(model->k>model->max_mem){
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dC[i]=(model->dB[(model->max_mem - 1)*model->nb_paramsMaintenance+i]*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->dC[i]);
            }
        }
        for(k=nk-1;k>nk2-1;k--) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dB[k*model->nb_paramsMaintenance+i]=model->dB[(k-1)*model->nb_paramsMaintenance+i];
            }
        } 
        for(k=nk2-1;k>0;k--) {
            for(i=0;i<model->nb_paramsMaintenance;i++) {
                model->dB[k*model->nb_paramsMaintenance+i]=(1-rho_ARA)*model->dB[(k-1)*model->nb_paramsMaintenance+i];
            }
            model->dB[k * (model->nb_paramsMaintenance)+id_params+1] -= model->B[k-1];
        }
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dB[i]=(1-rho_ARA)*model->dA[i];
        }
        model->dB[id_params+1]-= model->A;
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            model->dA[i] = pow(rho_QR,f->eval(K)-f->eval(K-1)) *  model->dA[i];
        }
        model->dA[id_params] = model->dA[id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1)) * (f->eval(K)-f->eval(K-1))/rho_QR * model->A;
    }
    if(model->k>model->max_mem){
        model->C=((model->B[model->max_mem - 1])*(model->time[model->k - model->max_mem]-model->time[model->k - model->max_mem - 1])+model->C);
    }
    for(k=nk-1;k>nk2-1;k--) {
        model->B[k]=model->B[k-1];
    }
    for(k=nk2-1;k>0;k--) {
        model->B[k]=(1-rho_ARA)*model->B[k-1];
    }
    model->B[0]=(1-rho_ARA)*model->A;
    model->A=pow(rho_QR,f->eval(K)-f->eval(K-1))*model->A;
    update_Vright(with_gradient,with_hessian);

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
        printf("WARNING: %s is not a proper maintenance model!\n",name.c_str());
    }
	return mm;
}

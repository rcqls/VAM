#ifndef RCPP_MLE_VAM_H
#define RCPP_MLE_VAM_H
#include <Rcpp.h>
#include "rcpp_vam_model.h"

using namespace Rcpp ;

class MLEVam {

public:

    MLEVam(List model_,List data_) {
        model=new VamModel(model_,data_);
        //Sum of d<?>S<??> defined in model
        dS1=new double[model->nb_paramsMaintenance+model->nb_paramsFamily-1+model->nb_paramsCov];
        dS2=new double[model->nb_paramsMaintenance+model->nb_paramsFamily-1];
        dS3=new double[model->nb_paramsMaintenance];
        if(model->nb_paramsCov>0) dS4=new double[model->nb_paramsCov];
        d2S1=new double[(model->nb_paramsMaintenance+model->nb_paramsFamily-1)*(model->nb_paramsMaintenance+model->nb_paramsFamily)/2];//inferior diagonal part of the hessian matrice by lines
        d2S2=new double[(model->nb_paramsMaintenance+model->nb_paramsFamily-1)*(model->nb_paramsMaintenance+model->nb_paramsFamily)/2];//inferior diagonal part of the hessian matrice by lines
        d2S3=new double[(model->nb_paramsMaintenance)*(model->nb_paramsMaintenance+1)/2];
        leftCensors=NULL;
        leftCensor=0;
    }

    ~MLEVam() {
        //DEBUG: printf("MLEVAM: %p, %p, %p\n",model,dS1,dS2);
        if(model->nb_paramsCov) delete[] dS4;
        delete[] dS1;
        delete[] dS2;
        delete[] dS3;
        delete[] d2S1;
        delete[] d2S2;
        delete[] d2S3;
        if(leftCensors != NULL) delete[] leftCensors;
        delete model;
    };

    void set_leftCensors(IntegerVector leftCensorsR) {
      leftCensors=new int[model->nb_system];
      for(int i=0; i<leftCensorsR.size();i++) leftCensors[i]=leftCensorsR[i];
    }

    void reset_leftCensors() {
      if(leftCensors != NULL) delete[] leftCensors;
      leftCensors=NULL;
      leftCensor=0;
    }

    void select_leftCensor(int i) {
      if(leftCensors != NULL) leftCensor=leftCensors[i];
    }

    void set_data(List data_) {
        model->set_data(data_);
    }

    NumericVector get_params() {
        return model->get_params();
    }

    void set_params(NumericVector pars) {
        model->set_params(pars);
    }

    void contrast_for_current_system() {
    	init_mle_vam_for_current_system(false,false);
  		int n=(model->time).size() - 1;
  		while(model->k < n) {
              //printf("  Time=%f, Type=%d\n",model->time[model->k+1],model->type[model->k+1]);
  			contrast_update_for_current_system(false,false);
  			// previous model for the next step
  			int type=model->type[model->k + 1 ];
  			if(type < 0) type=0;
  			//model->indMode = (type < 0 ? 0 : type);
  			model->models->at(type)->update(false,false);
  		}
        //model updated for current system: S1,S2,S0
        S1 += model->S1 *(model->nb_paramsCov > 0 ? exp(model->sum_cov) : 1.0 );S2 += model->S2; S0 += model->S0; S3 += model->S3;
        if(model->nb_paramsCov>0) S4 += model->S0 * model->sum_cov;
        //printf("Conclusion : S1=%f, S2=%f, S0=%f\n",model->S1,model->S2,model->S0);

    }

    NumericVector contrast(NumericVector param, bool alpha_fixed=false) {
        NumericVector res(1);
        double alpha=param[0];//save current value of alpha

        param[0]=1; //Rmk: alpha replaces param[0] => a bit weird!

        init_mle_vam(false,false);
        model->set_params(param);
        //printf("System %d\n",1);
        model->select_data(0);
        if(model->nb_paramsCov > 0) model->select_current_system(0,true);
        select_leftCensor(0);
        contrast_for_current_system();
        //only if multi-system
        for(int i=1;i<model->nb_system;i++) {
            //printf("System %d\n",i+1);
            model->select_data(i);
            if(model->nb_paramsCov > 0) model->select_current_system(i,true);
            select_leftCensor(i);
            contrast_for_current_system();
        }

        //DEBUG: printf("alpha=%lf,S0=%lf,S1=%lf,S2=%lf,S3=%lf\n",alpha,S0,S1,S2,S3);
        // log-likelihood (with constant +S0*(log(S0)-1))
        if(!alpha_fixed) {
          res[0]=-log(S1) * S0 + S2 +S0*(log(S0)-1)+S3;
        } else {
          res[0]=log(alpha)*S0+S2-alpha*S1+S3;
        }
        if(model->nb_paramsCov>0) res[0] += S4;

        param[0]=alpha;//LD:changed for bayesian
        return res;
        //return res[0]==R_NaN ? R_NegInf : res;
    }

    void gradient_for_current_system() {
        int i,j;
    	init_mle_vam_for_current_system(true,false);
    	int n=(model->time).size() - 1;
    	while(model->k < n) {
    		gradient_update_for_current_system();
    		int type=model->type[model->k + 1 ];
			if(type < 0) type=0;
			//model->indMode = (type < 0 ? 0 : type);
			model->models->at(type)->update(true,false);
    	}
        //model updated for current system: S1,S2,S0,dS1,dS2
        S1 += model->S1;S2 += model->S2; S0 += model->S0; S3 += model->S3;
        if(model->nb_paramsCov>0) S4 += model->S0 * model->sum_cov;
        //precomputation of covariate term to multiply (in fact just exp)
        for(i=0;i<(model->nb_paramsMaintenance);i++) {
            dS1[i] += model->dS1[i] * (model->nb_paramsCov > 0 ? exp(model->sum_cov) : 1.0); dS2[i] += model->dS2[i]; dS3[i] += model->dS3[i];

        }
        for(i=(model->nb_paramsMaintenance);i<(model->nb_paramsMaintenance+model->nb_paramsFamily-1);i++) {
            dS1[i] += model->dS1[i] * (model->nb_paramsCov > 0 ? exp(model->sum_cov) : 1.0); dS2[i] += model->dS2[i];
        }
        for(j=0;j<model->nb_paramsCov;i++,j++) {
            dS1[i] += model->dS1[i] * (model->nb_paramsCov > 0 ? model->get_covariate(j) * exp(model->sum_cov) : 1.0); dS2[i] += model->dS2[i];
            dS4[j] += model->S0 * model->get_covariate(j);
        }
    }

    NumericVector gradient(NumericVector param, bool alpha_fixed=false) {
        NumericVector res(model->nb_paramsMaintenance+model->nb_paramsFamily+model->nb_paramsCov);
        double alpha=param[0];//save current value of alpha

        param[0]=1;

        init_mle_vam(true,false);
        model->set_params(param);
        model->select_data(0);
        if(model->nb_paramsCov > 0) model->select_current_system(0,true);
        select_leftCensor(0);
        gradient_for_current_system();

        //only if multi-system
        for(int i=1;i<model->nb_system;i++) {
            model->select_data(i);
            if(model->nb_paramsCov > 0) model->select_current_system(i,true);
            select_leftCensor(i);
            gradient_for_current_system();
        }

        //compute gradient
        if(!alpha_fixed) {
            res[0] = 0;
            for(int i=0;i<(model->nb_paramsFamily-1);i++) {
                res[i+1] = -dS1[i]/S1 * S0 + dS2[i];
            }
            for(int i=(model->nb_paramsFamily-1);i<(model->nb_paramsMaintenance+model->nb_paramsFamily-1);i++) {
                res[i+1] = -dS1[i]/S1 * S0 + dS2[i]+dS3[i-(model->nb_paramsFamily-1)];
            }
        } else {

            res[0] = S0/alpha-S1;
            for(int i=0;i<(model->nb_paramsFamily-1);i++) {
                res[i+1] = -dS1[i]*alpha + dS2[i];
            }
            for(int i=(model->nb_paramsFamily-1);i<(model->nb_paramsMaintenance+model->nb_paramsFamily-1);i++) {
                res[i+1] = -dS1[i]*alpha + dS2[i]+dS3[i-(model->nb_paramsFamily-1)];
            }
        }
        param[0]=alpha;//LD:changed for bayesian
        return res;
    }

    void hessian_for_current_system() {
        int j;
        init_mle_vam_for_current_system(true,true);
        int n=(model->time).size() - 1;
        while(model->k < n) {
            hessian_update_for_current_system();
            int type=model->type[model->k + 1 ];
            if(type < 0) type=0;
            //model->indMode = (type < 0 ? 0 : type);
            model->models->at(type)->update(true,true);
        }
        //model updated for current system: S1,S2,S0,dS1,dS2,d2S1,d2S2
        S1 += model->S1;S2 += model->S2; S0 += model->S0; S3 += model->S3;
        for(int i=0;i<(model->nb_paramsMaintenance);i++) {
            dS1[i] += model->dS1[i]; dS2[i] += model->dS2[i]; dS3[i] += model->dS3[i];
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                d2S1[i*(i+1)/2+j] += model->d2S1[i*(i+1)/2+j]; d2S2[i*(i+1)/2+j] += model->d2S2[i*(i+1)/2+j]; d2S3[i*(i+1)/2+j] += model->d2S3[i*(i+1)/2+j];
            }
        }
        for(int i=(model->nb_paramsMaintenance);i<(model->nb_paramsMaintenance+model->nb_paramsFamily-1);i++) {
            dS1[i] += model->dS1[i]; dS2[i] += model->dS2[i];
            for(j=0;j<=i;j++) {
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                d2S1[i*(i+1)/2+j] += model->d2S1[i*(i+1)/2+j]; d2S2[i*(i+1)/2+j] += model->d2S2[i*(i+1)/2+j];
            }
        }
    }

    NumericMatrix hessian(NumericVector param, bool alpha_fixed=false) {
        int j;
        NumericMatrix res(model->nb_paramsMaintenance+model->nb_paramsFamily+model->nb_paramsCov,model->nb_paramsMaintenance+model->nb_paramsFamily+model->nb_paramsCov);
        double alpha=param[0];//save current value of alpha

        param[0]=1;

        init_mle_vam(true,true);
        model->set_params(param);
        model->select_data(0);
        if(model->nb_paramsCov > 0) model->select_current_system(0,true);
        select_leftCensor(0);
        hessian_for_current_system();

        //only if multi-system
        for(int i=1;i<model->nb_system;i++) {
            model->select_data(i);
            if(model->nb_paramsCov > 0) model->select_current_system(i,true);
            select_leftCensor(i);
            hessian_for_current_system();
        }

        //compute hessian
        if(!alpha_fixed) {
            res(0,0) = 0;
            for(int i=0;i<(model->nb_paramsFamily-1);i++) {
                res(0,i+1) = 0;
                res(i+1,0) = 0;
                res(i+1,i+1) = pow(dS1[i],2)/pow(S1,2) * S0-d2S1[i*(i+1)/2+i]/S1 * S0 + d2S2[i*(i+1)/2+i];
                for(j=0;j<i;j++) {
                    //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                    res(i+1,j+1) = dS1[i]*dS1[j]/pow(S1,2) * S0-d2S1[i*(i+1)/2+j]/S1 * S0 + d2S2[i*(i+1)/2+j];
                    res(j+1,i+1) = res(i+1,j+1);
                }
            }
            for(int i=(model->nb_paramsFamily-1);i<(model->nb_paramsMaintenance+model->nb_paramsFamily-1);i++) {
                res(0,i+1) = 0;
                res(i+1,0) = 0;
                res(i+1,i+1) = pow(dS1[i],2)/pow(S1,2) * S0-d2S1[i*(i+1)/2+i]/S1 * S0 + d2S2[i*(i+1)/2+i] + d2S3[(i-(model->nb_paramsFamily-1))*(i-(model->nb_paramsFamily-1)+1)/2+i-(model->nb_paramsFamily-1)];
                for(j=0;j<(model->nb_paramsFamily-1);j++) {
                    //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                    res(i+1,j+1) = dS1[i]*dS1[j]/pow(S1,2) * S0-d2S1[i*(i+1)/2+j]/S1 * S0 + d2S2[i*(i+1)/2+j];
                    res(j+1,i+1) = res(i+1,j+1);
                }
                for(j=(model->nb_paramsFamily-1);j<i;j++) {
                    //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                    res(i+1,j+1) = dS1[i]*dS1[j]/pow(S1,2) * S0-d2S1[i*(i+1)/2+j]/S1 * S0 + d2S2[i*(i+1)/2+j]+ d2S3[(i-(model->nb_paramsFamily-1))*(i-(model->nb_paramsFamily-1)+1)/2+j-(model->nb_paramsFamily-1)];
                    res(j+1,i+1) = res(i+1,j+1);
                }
            }
        } else {

            res(0,0) = -S0/pow(alpha,2);
            for(int i=0;i<(model->nb_paramsFamily-1);i++) {
                res(0,i+1) = -dS1[i];
                res(i+1,0) = -dS1[i];
                res(i+1,i+1) = d2S2[i*(i+1)/2+i]-alpha*d2S1[i*(i+1)/2+i];
                for(j=0;j<i;j++) {
                    //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                    res(i+1,j+1) = d2S2[i*(i+1)/2+j]-alpha*d2S1[i*(i+1)/2+j];
                    res(j+1,i+1) = res(i+1,j+1);
                }
            }
            for(int i=(model->nb_paramsFamily-1);i<(model->nb_paramsMaintenance+model->nb_paramsFamily-1);i++) {
                res(0,i+1) = -dS1[i];
                res(i+1,0) = -dS1[i];
                res(i+1,i+1) = d2S2[i*(i+1)/2+i]-alpha*d2S1[i*(i+1)/2+i]+d2S3[(i-(model->nb_paramsFamily-1))*(i- (model->nb_paramsFamily-1)+1)/2+i-(model->nb_paramsFamily-1)];
                for(j=0;j<(model->nb_paramsFamily-1);j++) {
                    //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                    res(i+1,j+1) = d2S2[i*(i+1)/2+j]-alpha*d2S1[i*(i+1)/2+j];
                    res(j+1,i+1) = res(i+1,j+1);
                }
                for(j=(model->nb_paramsFamily-1);j<i;j++) {
                    //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                    res(i+1,j+1) = d2S2[i*(i+1)/2+j]-alpha*d2S1[i*(i+1)/2+j]+d2S3[(i-(model->nb_paramsFamily-1))*(i- (model->nb_paramsFamily-1)+1)/2+j-(model->nb_paramsFamily-1)];
                    res(j+1,i+1) = res(i+1,j+1);
                }
            }
        }
        param[0]=alpha;//LD:changed for bayesian
        return res;
    }

    NumericVector get_alpha_est(NumericVector param) {
        NumericVector res(1);
        contrast(param); //To compute S1 and S0
        res[0]=S0/S1;
        return res;
    }

    VamModel* get_model() {
    	return model;
    }


    //delegate from model cache!
    List get_virtual_age_infos(double by,double from, double to) {
        return model->get_virtual_age_infos(by,from,to);
    }

    DataFrame get_selected_data(int i) {
        return model->get_selected_data(i);
    }

private:

	VamModel* model;

    double S1, S2, S0, S3, S4, *dS1, *dS2, *dS3, *dS4;//Accumulator!
    double *d2S1, *d2S2, *d2S3;//Accumulator!
    int* leftCensors; //CAREFUL: this is a vector of indices!
    int leftCensor; //leftCensor for current system

    void init_mle_vam(bool with_gradient,bool with_hessian) {
        int i;
        int j;
        S1 = 0; S2 = 0; S0 = 0; S3=0; S4=0;
        if(with_hessian) {
            for(i=0;i<(model->nb_paramsMaintenance);i++) {
                dS1[i] = 0; dS2[i] = 0; dS3[i] = 0; dS4[i] = 0;
                for(j=0;j<=i;j++) {
                    //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                    d2S1[i*(i+1)/2+j] = 0; d2S2[i*(i+1)/2+j] = 0; d2S3[i*(i+1)/2+j] = 0;
                }
            }
            for(i=(model->nb_paramsMaintenance);i<(model->nb_paramsMaintenance+model->nb_paramsFamily-1);i++) {
                dS1[i] = 0; dS2[i] = 0;
                for(j=0;j<=i;j++) {
                    //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                    d2S1[i*(i+1)/2+j] = 0; d2S2[i*(i+1)/2+j] = 0;
                }
            }
        }
        else if(with_gradient) {
            for(i=0;i<(model->nb_paramsMaintenance);i++) {
                dS1[i] = 0; dS2[i] = 0; dS3[i] = 0;
            }
            for(i=(model->nb_paramsMaintenance);i<(model->nb_paramsMaintenance+model->nb_paramsFamily-1);i++) {
                dS1[i] = 0; dS2[i] = 0;
            }
        }
    }

    void init_mle_vam_for_current_system(bool with_gradient,bool with_hessian) {
        int i;
        int j;

        for(i=0;i<model->nbPM + 1;i++) model->models->at(i)->init();

    	model->Vright = 0; //100000.;
        model->A=1;
    	model->k=0;
    	model->idMod=0; //id of current model
    	model->S1 = 0;
    	model->S2 = 0;
        model->S3 = 0;
        model->S4 = 0;
    	model->S0 = 0;for(i=0;i<model->type.size();i++) if(model->type[i] < 0) (model->S0) += 1; //TO COMPUTE from model->type
        if(with_hessian) {
            for(i=0;i<(model->nb_paramsMaintenance);i++) {
                model->dS1[i]=0;
                model->dS2[i]=0;
                model->dS3[i]=0; 
                model->dVright[i]=0;
                model->dA[i]=0;
                //for(k=0;k<model->mu;k++) model->dVR_prec[k*model->nb_paramsMaintenance+i]=0;
                for(j=0;j<=i;j++) {
                    //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                    model->d2S1[i*(i+1)/2+j]=0;
                    model->d2S2[i*(i+1)/2+j]=0;
                    model->d2S3[i*(i+1)/2+j]=0;
                    model->d2Vright[i*(i+1)/2+j]=0;
                    model->d2A[i*(i+1)/2+j]=0;
                    //for(k=0;k<model->mu;k++) model->d2VR_prec[k*(model->nb_paramsMaintenance)*(model->nb_paramsMaintenance+1)/2+i*(i+1)/2+j]=0;
                }
            }
            for(i=(model->nb_paramsMaintenance);i<(model->nb_paramsMaintenance+model->nb_paramsFamily-1);i++) {
                model->dS1[i]=0;
                model->dS2[i]=0;
                for(j=0;j<=i;j++) {
                    //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                    model->d2S1[i*(i+1)/2+j]=0;
                    model->d2S2[i*(i+1)/2+j]=0;
                }
            }
        }
        else if(with_gradient) {
            for(i=0;i<(model->nb_paramsMaintenance);i++) {
                model->dS1[i]=0;
                model->dS2[i]=0;
                model->dS3[i]=0;
                model->dVright[i]=0;
                model->dA[i]=0;
            }
            for(i=(model->nb_paramsMaintenance);i<(model->nb_paramsMaintenance+model->nb_paramsFamily-1);i++) {
                model->dS1[i]=0;
                model->dS2[i]=0;
            }

            for(j=0;j<model->nb_paramsCov;i++,j++) {
                model->dS1[i]=0;
                model->dS4[j]=0;
            }
    	}
    }

    void contrast_update_for_current_system(bool with_gradient, bool with_hessian) {
    	model->update_Vleft(with_gradient,with_hessian);
    	model->hVleft=model->family->hazardRate(model->Vleft);
    	model->indType = ((model->type)[(model->k) + 1] < 0 ? 1.0 : 0.0);
    	// printf("HVleft:%d,%lf,%lf\n",model->k,model->Vleft,model->family->cumulative_hazardRate(model->Vleft));
    	// printf("HVright:%lf,%lf\n",model->Vright,model->family->cumulative_hazardRate(model->Vright));
    	// printf("S1:%lf\n",model->S1);
    	// printf("indType,S2,hVleft:%lf,%lf,%lf\n",model->indType,model->S1,model->hVleft);
    	if(model->k >= leftCensor) model->S1 += model->family->cumulative_hazardRate(model->Vleft) - model->family->cumulative_hazardRate(model->Vright);
    	model->S2 += log(model->hVleft)* model->indType;
        model->S3 += log(model->A)* model->indType; 
    	//for(int i=0;i<(model->nbPM)+2;i++) model->dS1[i] += cdVleft[i] - cdVright[i];
    	//model->dS1 += (models->at(0))
    }

    void gradient_update_for_current_system() {
        int i;
    	contrast_update_for_current_system(true,false);

        double *cumhVright_param_derivative=model->family->cumulative_hazardRate_param_derivative(model->Vright,true);
        double *cumhVleft_param_derivative=model->family->cumulative_hazardRate_param_derivative(model->Vleft,false);
        double *hVleft_param_derivative=model->family->hazardRate_param_derivative(model->Vleft,false);
        for(i=0;i<model->nb_paramsFamily-1;i++){
            if(model->k >= leftCensor) model->dS1[i] +=  cumhVleft_param_derivative[i]-cumhVright_param_derivative[i] ;
            model->dS2[i] += hVleft_param_derivative[i]/model->hVleft*model->indType ;
        }
    	double hVright=model->family->hazardRate(model->Vright);
    	double dhVleft=model->family->hazardRate_derivative(model->Vleft);
    	  //printf("k:%d,hVright:%lf,dhVleft:%lf,indType:%lf\n",model->k,hVright,dhVleft,model->indType);
    	for(i=0;i<model->nb_paramsMaintenance;i++) {
    		if(model->k >= leftCensor) model->dS1[i+model->nb_paramsFamily-1] += model->hVleft * model->dVleft[i] - hVright * model->dVright[i];
    		//printf("dS1[%d]=(%lf,%lf,%lf),%lf,",i+1,model->hVleft,model->dVleft[i],model->dVright[i],model->dS1[i+1]);
    		model->dS2[i+model->nb_paramsFamily-1] +=  dhVleft * model->dVleft[i]/model->hVleft * model->indType;
    		//printf("dS2[%d]=%lf,",i+1,model->dS2[i+1]);
            model->dS3[i] +=  model->dA[i]/model->A * model->indType;
    	}
    	//printf("\n");
    }

    void hessian_update_for_current_system() {
        int i;
        int j;
        contrast_update_for_current_system(true,true);

        double *cumhVright_param_derivative=model->family->cumulative_hazardRate_param_derivative(model->Vright,true);
        double *cumhVleft_param_derivative=model->family->cumulative_hazardRate_param_derivative(model->Vleft,false);
        double *hVleft_param_derivative=model->family->hazardRate_param_derivative(model->Vleft,false);
        double *cumhVright_param_2derivative=model->family->cumulative_hazardRate_param_2derivative(model->Vright,true);
        double *cumhVleft_param_2derivative=model->family->cumulative_hazardRate_param_2derivative(model->Vleft,false);
        double *hVleft_param_2derivative=model->family->hazardRate_param_2derivative(model->Vleft);
        for(i=0;i<model->nb_paramsFamily-1;i++){
            if(model->k >= leftCensor) model->dS1[i] +=  cumhVleft_param_derivative[i]-cumhVright_param_derivative[i] ;
            model->dS2[i] += hVleft_param_derivative[i]/model->hVleft*model->indType ;
            for(j=0;j<=i;j++) {
                if(model->k >= leftCensor) model->d2S1[i*(i+1)/2+j] += cumhVleft_param_2derivative[i*(i+1)/2+j]-cumhVright_param_2derivative[i*(i+1)/2+j];
                model->d2S2[i*(i+1)/2+j] += (hVleft_param_2derivative[i*(i+1)/2+j]/model->hVleft -hVleft_param_derivative[i]*hVleft_param_derivative[j]/pow(model->hVleft,2))*model->indType;
                }
            }
        double hVright=model->family->hazardRate(model->Vright);
        double dhVleft=model->family->hazardRate_derivative(model->Vleft);
        double dhVright=model->family->hazardRate_derivative(model->Vright);
        double *hVright_param_derivative=model->family->hazardRate_param_derivative(model->Vright,true);
        double *dhVleft_param_derivative=model->family->hazardRate_derivative_param_derivative(model->Vleft);
        double d2hVleft=model->family->hazardRate_2derivative(model->Vleft);
        //printf("k:%d,hVright:%lf,dhVleft:%lf,indType:%lf\n",model->k,hVright,dhVleft,model->indType);
        for(i=0;i<model->nb_paramsMaintenance;i++) {
            if(model->k >= leftCensor) model->dS1[i+model->nb_paramsFamily-1] += model->hVleft * model->dVleft[i] - hVright * model->dVright[i];
            //printf("dS1[%d]=(%lf,%lf,%lf),%lf,",i+1,model->hVleft,model->dVleft[i],model->dVright[i],model->dS1[i+1]);
            model->dS2[i+model->nb_paramsFamily-1] +=  dhVleft * model->dVleft[i]/model->hVleft * model->indType;
            //printf("dS2[%d]=%lf,",i+1,model->dS2[i+1]);
            //column 0 and i+1 corresponds to the line indice of (inferior diagonal part of) the hessian matrice
            model->dS3[i] +=  model->dA[i]/model->A * model->indType;
            for(j=0;j<model->nb_paramsFamily-1;j++){
                if(model->k >= leftCensor) model->d2S1[(i+model->nb_paramsFamily-1)*(i+model->nb_paramsFamily)/2+j] += hVleft_param_derivative[j] * model->dVleft[i] - hVright_param_derivative[j] * model->dVright[i];
                model->d2S2[(i+model->nb_paramsFamily-1)*(i+model->nb_paramsFamily)/2+j] +=  dhVleft_param_derivative[j] * model->dVleft[i]/model->hVleft * model->indType - hVleft_param_derivative[j]*dhVleft * model->dVleft[i]/pow(model->hVleft,2) * model->indType;
            }
            for(j=0;j<=i;j++){
                //i+1 and j+1(<=i+1) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                if(model->k >= leftCensor) model->d2S1[(i+model->nb_paramsFamily-1)*(i+model->nb_paramsFamily)/2+j+model->nb_paramsFamily-1] += dhVleft*model->dVleft[i]*model->dVleft[j] + model->hVleft * model->d2Vleft[i*(i+1)/2+j] - dhVright*model->dVright[i]*model->dVright[j] - hVright * model->d2Vright[i*(i+1)/2+j];
                model->d2S2[(i+model->nb_paramsFamily-1)*(i+model->nb_paramsFamily)/2+j+model->nb_paramsFamily-1] += ( model->dVleft[i]*model->dVleft[j]*(d2hVleft/model->hVleft-pow(dhVleft/model->hVleft,2)) + dhVleft * model->d2Vleft[i*(i+1)/2+j]/model->hVleft )* model->indType;
                model->d2S3[i*(i+1)/2+j] += (model->d2A[i*(i+1)/2+j]/model->A -model->dA[i]*model->dA[j]/pow(model->A,2))* model->indType;
            }
        }
        //printf("\n");
    }

};

#endif //RCPP_MLE_VAM_H

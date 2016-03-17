#ifndef RCPP_MLE_VAM_H
#define RCPP_MLE_VAM_H
#include <Rcpp.h>
#include "rcpp_vam_model.h"

using namespace Rcpp ;

class MLEVam { 

public:

    MLEVam(List model_,List data_) {
        model=new VamModel(model_,data_);
        dS1=new double[model->nb_paramsMaintenance+model->nb_paramsFamily-1];//LD3
        dS2=new double[model->nb_paramsMaintenance+model->nb_paramsFamily-1];//LD3
        d2S1=new double[(model->nb_paramsMaintenance+model->nb_paramsFamily-1)*(model->nb_paramsMaintenance+model->nb_paramsFamily)/2];//inferior diagonal part of the hessian matrice by lines//LD//LD3
        d2S2=new double[(model->nb_paramsMaintenance+model->nb_paramsFamily-1)*(model->nb_paramsMaintenance+model->nb_paramsFamily)/2];//inferior diagonal part of the hessian matrice by lines//LD//LD3
    }

    ~MLEVam() {
        //DEBUG: printf("MLEVAM: %p, %p, %p\n",model,dS1,dS2);
        delete model;
        delete[] dS1;
        delete[] dS2;
        delete[] d2S1;//LD
        delete[] d2S2;//LD
    };

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
    	init_mle_vam_for_current_system(false,false);//LD:second false
		int n=(model->time).size() - 1;
		while(model->k < n) {
            //printf("  Time=%f, Type=%d\n",model->time[model->k+1],model->type[model->k+1]);
			contrast_update_for_current_system(false,false);//LD:second false
			// previous model for the next step
			int type=model->type[model->k + 1 ];
			if(type < 0) type=0;
			//model->indMode = (type < 0 ? 0 : type);
			model->models->at(type)->update(false,false);//LD:second false
		}
        //model updated for current system: S1,S2,S3
        S1 += model->S1;S2 += model->S2; S3 += model->S3;
        //printf("Conclusion : S1=%f, S2=%f, S3=%f\n",model->S1,model->S2,model->S3);

    }

    NumericVector contrast(NumericVector param, bool alpha_fixed=false) {//LD2
        NumericVector res(1);
        double alpha=param[0];//LD2:save current value of alpha

        param[0]=1;//LD2

        init_mle_vam(false,false);//LD:second false
        model->set_params(param);
        //printf("System %d\n",1);
        model->select_data(0);
        contrast_for_current_system();
        //only if multi-system
        for(int i=1;i<model->nb_system;i++) {
            //printf("System %d\n",i+1);
            model->select_data(i);
            contrast_for_current_system();
        }

        // log-likelihood (with constant +S3*(log(S3)-1))
        if(!alpha_fixed)//LD2
            {res[0]=-log(S1) * S3 + S2 +S3*(log(S3)-1);}
        else//LD2
            {res[0]=log(alpha)*S3+S2-alpha*S1;}//LD2

        return res;
        //return res[0]==R_NaN ? R_NegInf : res;
    }

    void gradient_for_current_system() {
    	init_mle_vam_for_current_system(true,false);//LD:false
    	int n=(model->time).size() - 1;
    	while(model->k < n) {
    		gradient_update_for_current_system();
    		int type=model->type[model->k + 1 ];
			if(type < 0) type=0;
			//model->indMode = (type < 0 ? 0 : type);
			model->models->at(type)->update(true,false);//LD:false
    	}
        //model updated for current system: S1,S2,S3,dS1,dS2
        S1 += model->S1;S2 += model->S2; S3 += model->S3;
        for(int i=0;i<(model->nb_paramsMaintenance+model->nb_paramsFamily-1);i++) {//LD3
            dS1[i] += model->dS1[i]; dS2[i] += model->dS2[i];
        }
    }

    NumericVector gradient(NumericVector param, bool alpha_fixed=false) {//LD2
        NumericVector res(model->nb_paramsMaintenance+model->nb_paramsFamily);//LD2//LD3
        double alpha=param[0];//LD2:save current value of alpha

        param[0]=1;//LD2

        init_mle_vam(true,false);//LD:false
        model->set_params(param);
        model->select_data(0);
        gradient_for_current_system();
        
        //only if multi-system
        for(int i=1;i<model->nb_system;i++) {
            model->select_data(i);
            gradient_for_current_system();
        }

        //compute gradient
        if(!alpha_fixed) {//LD2
            res[0] = 0;//LD2
            for(int i=0;i<(model->nb_paramsMaintenance+model->nb_paramsFamily-1);i++) {//LD3
                res[i+1] = -dS1[i]/S1 * S3 + dS2[i];//LD2
            }
        } else {//LD2

            res[0] = S3/alpha-S1;//LD2
            for(int i=0;i<(model->nb_paramsMaintenance+model->nb_paramsFamily-1);i++) {//LD2//LD3
                res[i+1] = dS2[i]-alpha*dS1[i];//LD2
            }//LD2
        }//LD2
        
        return res;
    }

    void hessian_for_current_system() {//LD
        int j;//LD
        init_mle_vam_for_current_system(true,true);//LD
        int n=(model->time).size() - 1;//LD
        while(model->k < n) {//LD
            hessian_update_for_current_system();//LD
            int type=model->type[model->k + 1 ];//LD
            if(type < 0) type=0;//LD
            //model->indMode = (type < 0 ? 0 : type);//LD
            model->models->at(type)->update(true,true);//LD
        }
        //model updated for current system: S1,S2,S3,dS1,dS2,d2S1,d2S2
        S1 += model->S1;S2 += model->S2; S3 += model->S3;//LD
        for(int i=0;i<(model->nb_paramsMaintenance+model->nb_paramsFamily-1);i++) {//LD//LD3
            dS1[i] += model->dS1[i]; dS2[i] += model->dS2[i];//LD
            for(j=0;j<=i;j++) {//LD
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                d2S1[i*(i+1)/2+j] += model->d2S1[i*(i+1)/2+j]; d2S2[i*(i+1)/2+j] += model->d2S2[i*(i+1)/2+j];//LD
            }//LD
        }//LD
    }//LD

    NumericMatrix hessian(NumericVector param, bool alpha_fixed=false) {//LD//LD2
        int j;//LD
        NumericMatrix res(model->nb_paramsMaintenance+model->nb_paramsFamily,model->nb_paramsMaintenance+model->nb_paramsFamily);//LD3//LD2
        double alpha=param[0];//LD2:save current value of alpha

        param[0]=1;//LD2

        init_mle_vam(true,true);//LD
        model->set_params(param);//LD
        model->select_data(0);
        hessian_for_current_system();//LD

        //only if multi-system
        for(int i=1;i<model->nb_system;i++) {//LD
            model->select_data(i);//LD
            hessian_for_current_system();//LD
        }//LD

        //compute hessian
        if(!alpha_fixed) {//LD2
            res(0,0) = 0;//LD2
            for(int i=0;i<(model->nb_paramsMaintenance+model->nb_paramsFamily-1);i++) {//LD//LD3
                res(0,i+1) = 0;//LD2
                res(i+1,0) = 0;//LD2
                res(i+1,i+1) = pow(dS1[i],2)/pow(S1,2) * S3-d2S1[i*(i+1)/2+i]/S1 * S3 + d2S2[i*(i+1)/2+i];//LD//LD2
                for(j=0;j<i;j++) {//LD
                    //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                    res(i+1,j+1) = dS1[i]*dS1[j]/pow(S1,2) * S3-d2S1[i*(i+1)/2+j]/S1 * S3 + d2S2[i*(i+1)/2+j];//LD//LD2
                    res(j+1,i+1) = res(i+1,j+1);//LD//LD2
                }//LD
            }//LD
        } else {//LD2
            
            res(0,0) = -S3/pow(alpha,2);//LD2
            for(int i=0;i<(model->nb_paramsMaintenance+model->nb_paramsFamily-1);i++) {//LD2//LD3
                res(0,i+1) = -dS1[i];//LD2
                res(i+1,0) = -dS1[i];//LD2
                res(i+1,i+1) = d2S2[i*(i+1)/2+i]-alpha*d2S1[i*(i+1)/2+i];//LD2
                for(j=0;j<i;j++) {//LD2
                    //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                    res(i+1,j+1) = d2S2[i*(i+1)/2+j]-alpha*d2S1[i*(i+1)/2+j];//LD2
                    res(j+1,i+1) = res(i+1,j+1);//LD2
                }//LD2
            }//LD2
        }//LD2
        
        return res;//LD
    }//LD

    NumericVector get_alpha_est(NumericVector param) {
        NumericVector res(1);
        contrast(param); //To compute S1 and S3
        res[0]=S3/S1;
        return res;
    }

    VamModel* get_model() {
    	return model;
    }


    //delegate from model cache!
    List get_virtual_age_infos(double by) {
        return model->get_virtual_age_infos(by);
    }

    DataFrame get_selected_data(int i) {
        return model->get_selected_data(i);
    }

private:

	VamModel* model;

    double S1, S2, S3, *dS1, *dS2;//Accumulator!
    double *d2S1, *d2S2;//LD:Accumulator!

    void init_mle_vam(bool with_gradient,bool with_hessian) {//LD
        int i;//LD
        int j;//LD
        S1 = 0; S2 = 0; S3 = 0;
        if(with_hessian) {//LD
            for(i=0;i<(model->nb_paramsMaintenance+model->nb_paramsFamily-1);i++) {//LD
                dS1[i] = 0; dS2[i] = 0;//LD
                for(j=0;j<=i;j++) {//LD
                    //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                    d2S1[i*(i+1)/2+j] = 0; d2S2[i*(i+1)/2+j] = 0;//LD
                }//LD
            }//LD
        }//LD
        else if(with_gradient) {//LD
            for(i=0;i<(model->nb_paramsMaintenance+model->nb_paramsFamily-1);i++) {//LD:enlever déclaration de int i
                dS1[i] = 0; dS2[i] = 0;
            }
        }
    }

    void init_mle_vam_for_current_system(bool with_gradient,bool with_hessian) {//LD)
        int i;//LD
        int j;//LD
    	model->Vright = 0; //100000.;
    	model->k=0;
    	model->idMod=0; //id of current model
    	model->S1 = 0;
    	model->S2 = 0;
    	model->S3 = 0;for(i=0;i<model->type.size();i++) if(model->type[i] < 0) (model->S3) += 1; //TO COMPUTE from model->type
        if(with_hessian) {//LD
            for(i=0;i<(model->nb_paramsMaintenance+model->nb_paramsFamily-1);i++) {//LD
                model->dS1[i]=0;//LD
                model->dS2[i]=0;//LD
                for(j=0;j<=i;j++) {//LD
                    //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                    model->d2S1[i*(i+1)/2+j]=0;//LD
                    model->d2S2[i*(i+1)/2+j]=0;//LD
                }
            }
            for(i=0;i<(model->nb_paramsMaintenance);i++) {//LD
                model->dVright[i]=0;//LD
                for(j=0;j<=i;j++) {//LD
                    //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                    model->d2Vright[i*(i+1)/2+j]=0;//LD
                }
            }
        }//LD
        else if(with_gradient) {//LD
            for(i=0;i<(model->nb_paramsMaintenance+model->nb_paramsFamily-1);i++) {
                model->dS1[i]=0;
                model->dS2[i]=0;
            }
            for(i=0;i<(model->nb_paramsMaintenance);i++) {//LD
                model->dVright[i]=0;//LD
            }
    	}
    }

    void contrast_update_for_current_system(bool with_gradient, bool with_hessian) {//LD
    	model->update_Vleft(with_gradient,with_hessian);//LD
    	model->hVleft=model->family->hazardRate(model->Vleft);
    	model->indType = ((model->type)[(model->k) + 1] < 0 ? 1.0 : 0.0);
    	// printf("HVleft:%d,%lf,%lf\n",model->k,model->Vleft,model->family->cumulative_hazardRate(model->Vleft));
    	// printf("HVright:%lf,%lf\n",model->Vright,model->family->cumulative_hazardRate(model->Vright));
    	// printf("S1:%lf\n",model->S1);
    	// printf("indType,S2,hVleft:%lf,%lf,%lf\n",model->indType,model->S1,model->hVleft);
    	model->S1 += model->family->cumulative_hazardRate(model->Vleft) - model->family->cumulative_hazardRate(model->Vright);
    	model->S2 += log(model->hVleft)* model->indType;
    	//for(int i=0;i<(model->nbPM)+2;i++) model->dS1[i] += cdVleft[i] - cdVright[i];
    	//model->dS1 += (models->at(0))
    }

    void gradient_update_for_current_system() {
        int i;
    	contrast_update_for_current_system(true,false);//LD
        
        NumericVector cumhVright_param_derivative=model->family->cumulative_hazardRate_param_derivative(model->Vright);//LD3
        NumericVector cumhVleft_param_derivative=model->family->cumulative_hazardRate_param_derivative(model->Vleft);//LD3
        NumericVector hVleft_param_derivative=model->family->hazardRate_param_derivative(model->Vleft);//LD3
        for(i=0;i<model->nb_paramsFamily-1;i++){//LD3
            model->dS1[i] +=  cumhVleft_param_derivative[i]-cumhVright_param_derivative[i] ;//LD3
            model->dS2[i] += hVleft_param_derivative[i]/model->hVleft*model->indType ;//LD3
        }//LD3
    	double hVright=model->family->hazardRate(model->Vright);
    	double dhVleft=model->family->hazardRate_derivative(model->Vleft);
    	  //printf("k:%d,hVright:%lf,dhVleft:%lf,indType:%lf\n",model->k,hVright,dhVleft,model->indType);
    	for(i=0;i<model->nb_paramsMaintenance;i++) {//LD3
    		model->dS1[i+model->nb_paramsFamily-1] += model->hVleft * model->dVleft[i] - hVright * model->dVright[i];//LD3
    		//printf("dS1[%d]=(%lf,%lf,%lf),%lf,",i+1,model->hVleft,model->dVleft[i],model->dVright[i],model->dS1[i+1]);
    		model->dS2[i+model->nb_paramsFamily-1] +=  dhVleft * model->dVleft[i]/model->hVleft * model->indType;//LD3
    		//printf("dS2[%d]=%lf,",i+1,model->dS2[i+1]);
    	}
    	//printf("\n");
    }

    void hessian_update_for_current_system() {//LD
        int i;
        int j;//LD
        contrast_update_for_current_system(true,true);//LD

        NumericVector cumhVright_param_derivative=model->family->cumulative_hazardRate_param_derivative(model->Vright);//LD3
        NumericVector cumhVleft_param_derivative=model->family->cumulative_hazardRate_param_derivative(model->Vleft);//LD3
        NumericVector hVleft_param_derivative=model->family->hazardRate_param_derivative(model->Vleft);//LD3
        NumericVector cumhVright_param_2derivative=model->family->cumulative_hazardRate_param_2derivative(model->Vright);//LD3
        NumericVector cumhVleft_param_2derivative=model->family->cumulative_hazardRate_param_2derivative(model->Vleft);//LD3
        NumericVector hVleft_param_2derivative=model->family->hazardRate_param_2derivative(model->Vleft);//LD3
        for(i=0;i<model->nb_paramsFamily-1;i++){//LD3
            model->dS1[i] +=  cumhVleft_param_derivative[i]-cumhVright_param_derivative[i] ;//LD3
            model->dS2[i] += hVleft_param_derivative[i]/model->hVleft*model->indType ;//LD3
            for(j=0;j<=i;j++) {//LD
                model->d2S1[i*(i+1)/2+j] += cumhVleft_param_2derivative[i*(i+1)/2+j]-cumhVright_param_2derivative[i*(i+1)/2+j];//LD
                model->d2S2[i*(i+1)/2+j] += (hVleft_param_2derivative[i*(i+1)/2+j]/model->hVleft -hVleft_param_derivative[i]*hVleft_param_derivative[j]/pow(model->hVleft,2))*model->indType;//LD
                }//LD3
            }//LD3
        double hVright=model->family->hazardRate(model->Vright);//LD
        double dhVleft=model->family->hazardRate_derivative(model->Vleft);//LD
        double dhVright=model->family->hazardRate_derivative(model->Vright);//LD
        NumericVector hVright_param_derivative=model->family->hazardRate_param_derivative(model->Vright);//LD
        NumericVector dhVleft_param_derivative=model->family->hazardRate_derivative_param_derivative(model->Vleft);//LD
        double d2hVleft=model->family->hazardRate_2derivative(model->Vleft);//LD
        //printf("k:%d,hVright:%lf,dhVleft:%lf,indType:%lf\n",model->k,hVright,dhVleft,model->indType);
        for(i=0;i<model->nb_paramsMaintenance;i++) {//LD
            model->dS1[i+model->nb_paramsFamily-1] += model->hVleft * model->dVleft[i] - hVright * model->dVright[i];//LD
            //printf("dS1[%d]=(%lf,%lf,%lf),%lf,",i+1,model->hVleft,model->dVleft[i],model->dVright[i],model->dS1[i+1]);
            model->dS2[i+model->nb_paramsFamily-1] +=  dhVleft * model->dVleft[i]/model->hVleft * model->indType;//LD
            //printf("dS2[%d]=%lf,",i+1,model->dS2[i+1]);
            //column 0 and i+1 corresponds to the line indice of (inferior diagonal part of) the hessian matrice
            for(j=0;j<model->nb_paramsFamily-1;j++){//LD3
                model->d2S1[(i+model->nb_paramsFamily-1)*(i+model->nb_paramsFamily)/2+j] += hVleft_param_derivative[j] * model->dVleft[i] - hVright_param_derivative[j] * model->dVright[i];//LD
                model->d2S2[(i+model->nb_paramsFamily-1)*(i+model->nb_paramsFamily)/2+j] +=  dhVleft_param_derivative[j] * model->dVleft[i]/model->hVleft * model->indType - hVleft_param_derivative[j]*dhVleft * model->dVleft[i]/pow(model->hVleft,2) * model->indType;//LD
            }//LD3
            for(j=0;j<=i;j++){//LD
                //i+1 and j+1(<=i+1) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2S1[(i+model->nb_paramsFamily-1)*(i+model->nb_paramsFamily)/2+j+model->nb_paramsFamily-1] += dhVleft*model->dVleft[i]*model->dVleft[j] + model->hVleft * model->d2Vleft[i*(i+1)/2+j] - dhVright*model->dVright[i]*model->dVright[j] - hVright * model->d2Vright[i*(i+1)/2+j];//LD
                model->d2S2[(i+model->nb_paramsFamily-1)*(i+model->nb_paramsFamily)/2+j+model->nb_paramsFamily-1] += ( model->dVleft[i]*model->dVleft[j]*(d2hVleft/model->hVleft-pow(dhVleft/model->hVleft,2)) + dhVleft * model->d2Vleft[i*(i+1)/2+j]/model->hVleft )* model->indType;//LD
            }//LD
        }//LD
        //printf("\n");
    }//LD

};

#endif //RCPP_SIM_VAM_H
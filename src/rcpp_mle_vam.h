#ifndef RCPP_MLE_VAM_H
#define RCPP_MLE_VAM_H
#include <Rcpp.h>
#include "rcpp_vam_model.h"

using namespace Rcpp ;

class MLEVam { 

public:

    MLEVam(List model_,List data_) {
        model=new VamModel(model_,data_);
        dS1=new double[model->nbPM+2];
        dS2=new double[model->nbPM+2];
        d2S1=new double[(model->nbPM+2)*(model->nbPM+3)/2];//inferior diagonal part of the hessian matrice by lines//LD
        d2S2=new double[(model->nbPM+2)*(model->nbPM+3)/2];//inferior diagonal part of the hessian matrice by lines//LD
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
			contrast_update_for_current_system(false,false);//LD:second false
			// previous model for the next step
			int type=model->type[model->k + 1 ];
			if(type < 0) type=0;
			//model->indMode = (type < 0 ? 0 : type);
			model->models->at(type)->update(false,false);//LD:second false
		}
        //model updated for current system: S1,S2,S3
        S1 += model->S1;S2 += model->S2; S3 += model->S3;

    }

    NumericVector contrast(NumericVector param) {
        NumericVector res(1);

        init_mle_vam(false,false);//LD:second false
        model->set_params(param);
        contrast_for_current_system();
        //only if multi-system
        for(int i=1;i<model->nb_system;i++) {
            //printf("system %d\n",i);
            model->select_data(i);
            contrast_for_current_system();
        }

        // log-likelihood (with constant +S3*(log(S3)-1))
        res[0]=-log(S1) * S3 + S2 +S3*(log(S3)-1);

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
        for(int i=0;i<model->nbPM + 2;i++) {
            dS1[i] += model->dS1[i]; dS2[i] += model->dS2[i];
        }
    }

    NumericVector gradient(NumericVector param) {
        NumericVector res(model->nbPM + 2);

        init_mle_vam(true,false);//LD:false
        model->set_params(param);
        gradient_for_current_system();
        
        //only if multi-system
        for(int i=1;i<model->nb_system;i++) {
            model->select_data(i);
            gradient_for_current_system();
        }

        //compute gradient
        for(int i=0;i<model->nbPM + 2;i++) {
            res[i] = -dS1[i]/S1 * S3 + dS2[i];
        }
        
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
        for(int i=0;i<model->nbPM + 2;i++) {//LD
            dS1[i] += model->dS1[i]; dS2[i] += model->dS2[i];//LD
            for(j=0;j<=i;j++) {//LD
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                d2S1[i*(i+1)/2+j] += model->d2S1[i*(i+1)/2+j]; d2S2[i*(i+1)/2+j] += model->d2S2[i*(i+1)/2+j];//LD
            }//LD
        }//LD
    }//LD

    NumericMatrix hessian(NumericVector param) {//LD
        int j;//LD
        NumericMatrix res(model->nbPM + 2,model->nbPM + 2);//LD

        init_mle_vam(true,true);//LD
        model->set_params(param);//LD
        hessian_for_current_system();//LD
        
        //only if multi-system
        for(int i=1;i<model->nb_system;i++) {//LD
            model->select_data(i);//LD
            hessian_for_current_system();//LD
        }//LD

        //compute hessian
        for(int i=0;i<model->nbPM + 2;i++) {//LD
            res(i,i) = pow(dS1[i],2)/pow(S1,2) * S3-d2S1[i*(i+1)/2+i]/S1 * S3 + d2S2[i*(i+1)/2+i];//LD
            for(j=0;j<i;j++) {//LD
                //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                res(i,j) = dS1[i]*dS1[j]/pow(S1,2) * S3-d2S1[i*(i+1)/2+j]/S1 * S3 + d2S2[i*(i+1)/2+j];//LD
                res(j,i) = res(i,j);//LD
            }//LD
        }//LD
        
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
            for(i=0;i<model->nbPM + 2;i++) {//LD
                dS1[i] = 0; dS2[i] = 0;//LD
                for(j=0;j<=i;j++) {//LD
                    //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                    d2S1[i*(i+1)/2+j] = 0; d2S2[i*(i+1)/2+j] = 0;//LD
                }//LD
            }//LD
        }//LD
        else if(with_gradient) {//LD
            for(i=0;i<model->nbPM + 2;i++) {//LD:enlever déclaration de int i
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
            model->dS1[0]=0;//LD
            model->dS2[0]=0;//LD
            model->d2S1[0]=0;//LD
            model->d2S2[0]=0;//LD
            for(i=0;i<(model->nbPM)+1;i++) {//LD
                model->dVright[i]=0;//LD
                model->dS1[i+1]=0;//LD
                model->dS2[i+2]=0;//LD
                model->d2S1[(i+1)*(i+2)/2]=0;//LD
                model->d2S2[(i+1)*(i+2)/2]=0;//LD
                for(j=0;j<=i;j++) {//LD
                    //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                    model->d2Vright[i*(i+1)/2+j]=0;//LD
                    model->d2S1[(i+1)*(i+2)/2+j+1]=0;//LD
                    model->d2S2[(i+1)*(i+2)/2+j+1]=0;//LD
                }//LD
            }//LD
        }//LD
        else if(with_gradient) {//LD
    		model->dVright[0]=0;
    		model->dS1[0]=0;model->dS1[1]=0;
    		model->dS2[0]=0;model->dS2[1]=0;
    		for(i=0;i<(model->nbPM);i++) {
    			model->dVright[i+1]=0;
    			model->dS1[i+2]=0;
    			model->dS2[i+2]=0;
    		}
    	}
        
    }

    void contrast_update_for_current_system(bool with_gradient, bool with_hessian) {//LD
    	model->update_Vleft(with_gradient,with_hessian);//LD
    	model->hVleft=model->family->density(model->Vleft);
    	model->indType = ((model->type)[(model->k) + 1] < 0 ? 1.0 : 0.0);
    	// printf("HVleft:%d,%lf,%lf\n",model->k,model->Vleft,model->family->cumulative_density(model->Vleft));
    	// printf("HVright:%lf,%lf\n",model->Vright,model->family->cumulative_density(model->Vright));
    	// printf("S1:%lf\n",model->S1);
    	// printf("indType,S2,hVleft:%lf,%lf,%lf\n",model->indType,model->S1,model->hVleft);
    	model->S1 += model->family->cumulative_density(model->Vleft) - model->family->cumulative_density(model->Vright);
    	model->S2 += log(model->hVleft)* model->indType;
    	//for(int i=0;i<(model->nbPM)+2;i++) model->dS1[i] += cdVleft[i] - cdVright[i];
    	//model->dS1 += (models->at(0))
    }

    void gradient_update_for_current_system() {
    	contrast_update_for_current_system(true,false);//LD
    	model->dS1[0] += model->family->cumulative_density_param_derivative(model->Vleft) - model->family->cumulative_density_param_derivative(model->Vright);
    	model->dS2[0] += model->family->density_param_derivative(model->Vleft)/model->hVleft*model->indType ;
    	double hVright=model->family->density(model->Vright);
    	double dhVleft=model->family->density_derivative(model->Vleft);
    	//printf("k:%d,hVright:%lf,dhVleft:%lf,indType:%lf\n",model->k,hVright,dhVleft,model->indType);
    	for(int i=0;i<(model->nbPM)+1;i++) {
    		model->dS1[i+1] += model->hVleft * model->dVleft[i] - hVright * model->dVright[i];
    		//printf("dS1[%d]=(%lf,%lf,%lf),%lf,",i+1,model->hVleft,model->dVleft[i],model->dVright[i],model->dS1[i+1]);
    		model->dS2[i+1] +=  dhVleft * model->dVleft[i]/model->hVleft * model->indType;
    		//printf("dS2[%d]=%lf,",i+1,model->dS2[i+1]);
    	}
    	//printf("\n");
    }

    void hessian_update_for_current_system() {//LD
        int j;//LD
        contrast_update_for_current_system(true,true);//LD
        model->dS1[0] += model->family->cumulative_density_param_derivative(model->Vleft) - model->family->cumulative_density_param_derivative(model->Vright);//LD
        model->dS2[0] += model->family->density_param_derivative(model->Vleft)/model->hVleft*model->indType ;//LD
        model->d2S1[0] += model->family->cumulative_density_param_2derivative(model->Vleft) - model->family->cumulative_density_param_2derivative(model->Vright);//LD
        model->d2S2[0] += model->family->density_param_2derivative(model->Vleft)/model->hVleft*model->indType - pow(model->family->density_param_derivative(model->Vleft)/model->hVleft,2)*model->indType;//LD
        double hVright=model->family->density(model->Vright);//LD
        double dhVleft=model->family->density_derivative(model->Vleft);//LD
        double dhVright=model->family->density_derivative(model->Vright);//LD
        double hVright_param_derivative=model->family->density_param_derivative(model->Vright);//LD
        double hVleft_param_derivative=model->family->density_param_derivative(model->Vleft);//LD
        double dhVleft_param_derivative=model->family->density_derivative_param_derivative(model->Vleft);//LD
        double d2hVleft=model->family->density_2derivative(model->Vleft);//LD
        //printf("k:%d,hVright:%lf,dhVleft:%lf,indType:%lf\n",model->k,hVright,dhVleft,model->indType);
        for(int i=0;i<(model->nbPM)+1;i++) {//LD
            model->dS1[i+1] += model->hVleft * model->dVleft[i] - hVright * model->dVright[i];//LD
            //printf("dS1[%d]=(%lf,%lf,%lf),%lf,",i+1,model->hVleft,model->dVleft[i],model->dVright[i],model->dS1[i+1]);
            model->dS2[i+1] +=  dhVleft * model->dVleft[i]/model->hVleft * model->indType;//LD
            //printf("dS2[%d]=%lf,",i+1,model->dS2[i+1]);
            //column 0 and i+1 corresponds to the line indice of (inferior diagonal part of) the hessian matrice
            model->d2S1[(i+1)*(i+2)/2] += hVleft_param_derivative * model->dVleft[i] - hVright_param_derivative * model->dVright[i];//LD
            model->d2S2[(i+1)*(i+2)/2] +=  dhVleft_param_derivative * model->dVleft[i]/model->hVleft * model->indType - hVleft_param_derivative*dhVleft * model->dVleft[i]/pow(model->hVleft,2) * model->indType;//LD
            for(j=0;j<=i;j++){//LD
                //i+1 and j+1(<=i+1) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model->d2S1[(i+1)*(i+2)/2+j+1] += dhVleft*model->dVleft[i]*model->dVleft[j] + model->hVleft * model->d2Vleft[i*(i+1)/2+j] - dhVright*model->dVright[i]*model->dVright[j] - hVright * model->d2Vright[i*(i+1)/2+j];//LD
                model->d2S2[(i+1)*(i+2)/2+j+1] += ( model->dVleft[i]*model->dVleft[j]*(d2hVleft/model->hVleft-pow(dhVleft/model->hVleft,2)) + dhVleft * model->d2Vleft[i*(i+1)/2+j]/model->hVleft )* model->indType;//LD
            }//LD
        }//LD
        //printf("\n");
    }//LD

};

#endif //RCPP_SIM_VAM_H
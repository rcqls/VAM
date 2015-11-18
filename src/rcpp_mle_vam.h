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
    }

    ~MLEVam() {
        //DEBUG: printf("MLEVAM: %p, %p, %p\n",model,dS1,dS2);
        delete model;
        delete[] dS1;
        delete[] dS2;
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
    	init_mle_vam_for_current_system(false);
		int n=(model->time).size() - 1;
		while(model->k < n) {
			contrast_update_for_current_system(false);
			// previous model for the next step
			int type=model->type[model->k + 1 ];
			if(type < 0) type=0;
			//model->indMode = (type < 0 ? 0 : type);
			model->models->at(type)->update(false);
		}
        //model updated for current system: S1,S2,S3
        S1 += model->S1;S2 += model->S2; S3 += model->S3;

    }

    NumericVector contrast(NumericVector param) {
        NumericVector res(1);

        init_mle_vam(false);
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
    }

    void gradient_for_current_system() {
    	init_mle_vam_for_current_system(true);
    	int n=(model->time).size() - 1;
    	while(model->k < n) {
    		gradient_update_for_current_system();
    		int type=model->type[model->k + 1 ];
			if(type < 0) type=0;
			//model->indMode = (type < 0 ? 0 : type);
			model->models->at(type)->update(true);
    	}
        //model updated for current system: S1,S2,S3,dS1,dS2
        S1 += model->S1;S2 += model->S2; S3 += model->S3;
        for(int i=0;i<model->nbPM + 2;i++) {
            dS1[i] += model->dS1[i]; dS2[i] += model->dS2[i];
        }
    }

    NumericVector gradient(NumericVector param) {
        NumericVector res(model->nbPM + 2);

        init_mle_vam(true);
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

    void init_mle_vam(bool with_gradient) {
        S1 = 0; S2 = 0; S3 = 0;
        if(with_gradient) {
            for(int i=0;i<model->nbPM + 2;i++) {
                dS1[i] = 0; dS2[i] = 0;
            }
        }
    }

    void init_mle_vam_for_current_system(bool with_gradient) {
    	model->Vright = 0; //100000.;
    	model->k=0;
    	model->idMod=0; //id of current model
    	model->S1 = 0;
    	model->S2 = 0;
    	model->S3 = 0;for(int i=0;i<model->type.size();i++) if(model->type[i] < 0) (model->S3) += 1; //TO COMPUTE from model->type
    	if(with_gradient) {
    		model->dVright[0]=0;
    		model->dS1[0]=0;model->dS1[1]=0;
    		model->dS2[0]=0;model->dS2[1]=0;
    		for(int i=0;i<(model->nbPM)+1;i++) {
    			model->dVright[i+1]=0;
    			model->dS1[i+2]=0;
    			model->dS2[i+2]=0;
    		}
    	}
    }

    void contrast_update_for_current_system(bool with_gradient) {
    	model->update_Vleft(with_gradient);
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
    	contrast_update_for_current_system(true);
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

};

#endif //RCPP_SIM_VAM_H
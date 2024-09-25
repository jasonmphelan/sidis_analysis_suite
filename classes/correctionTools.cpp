#include "correctionTools.h"
#include "TFile.h"
#include "TF1.h"
#define CORR_PATH _DATA


correctionTools::correctionTools(int _mode){
    SetCorrectionMode (_mode);
}

correctionTools::~correctionTools(){}
	
void loadParameters(){
	loadWeightParameters();
	loadkTopiParameters();
	loadpiTokParameters();
}	

void loadWeightParameters(){
	TFile weightFile((TString) _Data + "/correctionFiles/"+ weight_name);
	for( int i = 0; i < 4; i++ ){
		for( int j = 0; j < 4; j++ ){
			TF1 * fp = (TF1*)weightFile.Get();
			weight_Parameters[0][i][j][0] = fp->GetParameter(0);
			weight_Parameters[0][i][j][1] = fp->GetParameter(1);
			weight_Parameters[0][i][j][2] = fp->GetParameter(2);
			weight_Parameters[0][i][j][3] = fp->GetParameter(3);
			
			TF1 * fm = (TF1*)weightFile.Get();
			weight_Parameters[1][i][j][0] = fm->GetParameter(0);
			weight_Parameters[1][i][j][1] = fm->GetParameter(1);
			weight_Parameters[1][i][j][2] = fm->GetParameter(2);
			weight_Parameters[1][i][j][3] = fm->GetParameter(3);
		}
	}
}

void loadkTopiParameters(){
	TFile weightFile((TString) _Data + "/correctionFiles/"+ kaon_To_pi_name);
	for( int i = 0; i < 4; i++ ){
		for( int j = 0; j < 4; j++ ){
			for( int k = 0; k < 4; k++ ){
				TF1 * fp = (TF1*)weightFile.Get();
				k_To_pi_Parameters[0][k][i][j][0] = fp->GetParameter(0);
				k_To_pi_Parameters[0][k][i][j][1] = fp->GetParameter(1);
				k_To_pi_Parameters[0][k][i][j][2] = fp->GetParameter(2);
				k_To_pi_Parameters[0][k][i][j][3] = fp->GetParameter(3);
			
				TF1 * fm = (TF1*)weightFile.Get();
				k_To_pi_Parameters[1][k][i][j][0] = fm->GetParameter(0);
				k_To_pi_Parameters[1][k][i][j][1] = fm->GetParameter(1);
				k_To_pi_Parameters[1][k][i][j][2] = fm->GetParameter(2);
				k_To_pi_Parameters[1][k][i][j][3] = fm->GetParameter(3);
			
			}
		}
	}
}

void loadpiTokParameters(){
	TFile weightFile((TString) _Data + "/correctionFiles/"+ pi_To_kaon_name);
	for( int i = 0; i < 4; i++ ){
		for( int j = 0; j < 4; j++ ){
			for( int k = 0; k < 4; k++ ){
				TF1 * fp = (TF1*)weightFile.Get();
				pi_To_k_Parameters[0][k][i][j][0] = fp->GetParameter(0);
				pi_To_k_Parameters[0][k][i][j][1] = fp->GetParameter(1);
				pi_To_k_Parameters[0][k][i][j][2] = fp->GetParameter(2);
				pi_To_k_Parameters[0][k][i][j][3] = fp->GetParameter(3);
			
				TF1 * fm = (TF1*)weightFile.Get();
				pi_To_k_Parameters[1][k][i][j][0] = fm->GetParameter(0);
				pi_To_k_Parameters[1][k][i][j][1] = fm->GetParameter(1);
				pi_To_k_Parameters[1][k][i][j][2] = fm->GetParameter(2);
				pi_To_k_Parameters[1][k][i][j][3] = fm->GetParameter(3);
			
			}
		}
	}
}

void setKinematics( double x, double q, double z, double p ){ 
	kin[0] = x;
	kin[1] = q;
	kin[2] = z;
	kin[3] = p;
}

double getCorrectionFactor( int type, int charge ){
	if( charge != 0 || charge != 1 ){
		std::cout<<"Invalid charge\n";
		return 0;
	}
	double weight = 0;
	double pi_To_k = 0;	
	double k_To_pi = 0;	
	
	for( int i = 0; i < 4; i++ ){
		double b = 0;
		for( int j = 0; j < 4; j++ ){
			double c = 0;
		
			for( int k = 0; k < 4; k++ ){
				c += weight_Parameters[type][i][j][k]*pow( x, k );
			}
			
			b+= c*pow(q, j);
		}
		weight+= b*pow(z, i);
	}
			
	if( type == 1 ) return weight;
	
	for( int i = 0; i < 4; i++ ){
		double b = 0;
		for( int j = 0; j < 4; j++ ){
			double c = 0;
		
			for( int k = 0; k < 4; k++ ){
				c += pi_To_k_Parameters[type][i][j][k]*pow( x, k );
			}
			
			b+= c*pow(q, j);
		}
		pi_To_k+= b*pow(z, i);
	}

	if( type == 2 ) return pi_To_k;
	
	for( int i = 0; i < 4; i++ ){
		double b = 0;
		for( int j = 0; j < 4; j++ ){
			double c = 0;
		
			for( int k = 0; k < 4; k++ ){
				c += k_To_pi_Parameters[type][i][j][k]*pow( x, k );
			}
			
			b+= c*pow(q, j);
		}
		k_To_pi+= b*pow(z, i);
	}

	if( type == 3 ) return k_To_pi;

	if( type > 3 ){
		std::cout<<"Invalid argument... returning 0\n";
		return 0;
	}
	
	return weight*pi_To_k;
}

			

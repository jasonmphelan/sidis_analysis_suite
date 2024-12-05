#include "correctionTools.h"
#include "TFile.h"
#include "TF1.h"
#include "TH3F.h"
#include "constants.h"
#include "cut_values.h"
#define CORR_PATH _DATA

using namespace constants;
using namespace cutVals;

correctionTools::correctionTools(int _mode){
    setCorrectionLevel (_mode);
}

correctionTools::~correctionTools(){}
	
void correctionTools::printFilePaths(){
	std::cout<< (TString) _DATA + "/correctionFiles/"+ weight_name<<std::endl;
	std::cout<< (TString) _DATA + "/correctionFiles/"+ kaon_To_pi_name<<std::endl;
	std::cout<< (TString) _DATA + "/correctionFiles/"+ pi_To_kaon_name<<std::endl;
}

void correctionTools::loadHistograms(){
//	TFile weightFile((TString) _DATA + "/correctionFiles/"+ weight_name);
//	TFile pi2kFile((TString) _DATA + "/correctionFiles/"+ kaon_To_pi_name);
//	TFile k2piFile((TString) _DATA + "/correctionFiles/"+ pi_To_kaon_name);
 	weightFile = new TFile((TString) _DATA + "/correctionFiles/"+ weight_name);
	pi2kFile = new TFile((TString) _DATA + "/correctionFiles/"+ kaon_To_pi_name);
	k2piFile = new TFile((TString) _DATA + "/correctionFiles/"+ pi_To_kaon_name);


	accCorrection[0] = (TH3F*)weightFile->Get( "hAccCorrectionP" );
	accCorrection[1] = (TH3F*)weightFile->Get( "hAccCorrectionM" );
	binMigration[0] = (TH3F*)weightFile->Get( "hBinMigrationP" );
	binMigration[1] = (TH3F*)weightFile->Get( "hBinMigrationM" );

	for( int i = 0; i < 4; i++ ){
		k_to_pi_Correction[0][i] = (TH3F*)pi2kFile->Get( Form( "hKaonCorrP_%i", i) );
		k_to_pi_Correction[1][i] = (TH3F*)pi2kFile->Get( Form( "hKaonCorrM_%i", i) );
		pi_to_k_Correction[0][i] = (TH3F*)k2piFile->Get( Form( "hKaonCorrP_%i", i) );
		pi_to_k_Correction[1][i] = (TH3F*)k2piFile->Get( Form( "hKaonCorrM_%i", i) );
	}

}

void correctionTools::loadFits(){
 	weightFile = new TFile((TString) _DATA + "/correctionFiles/"+ weight_name);
	pi2kFile = new TFile((TString) _DATA + "/correctionFiles/"+ kaon_To_pi_name);
	k2piFile = new TFile((TString) _DATA + "/correctionFiles/"+ pi_To_kaon_name);
	
	for( int x = 0; x < bins_xB; x++ ){
		for( int q = 0; q < bins_Q2; q++ ){
			weightFit[0][x][q] = (TF1*)weightFile->Get( Form("fitPip_%i_%i", q, x) );
			weightFit[1][x][q] = (TF1*)weightFile->Get( Form("fitPim_%i_%i", q, x) );

			for( int p = 0; p < 4; p++ ){
				k2piFit[0][x][q][p] = (TF1*)k2piFile->Get( Form("fitPip_%i_%i_%i", p, q, x) );
				k2piFit[1][x][q][p] = (TF1*)k2piFile->Get( Form("fitPip_%i_%i_%i", p, q, x) );
				
				pi2kFit[0][x][q][p] = (TF1*)pi2kFile->Get( Form("fitPip_%i_%i_%i", p, q, x) );
				pi2kFit[1][x][q][p] = (TF1*)pi2kFile->Get( Form("fitPip_%i_%i_%i", p, q, x) );
			}
		}
	}
}

void correctionTools::loadParameters(){
	printFilePaths();
	loadWeightParameters();
	std::cout<<"Loaded Weight Parameters\n";
	loadkTopiParameters();
	std::cout<<"Loaded k-to-pi Parameters\n";
	loadpiTokParameters();
	//std::cout<<"Loaded pi-to-k Parameters\n";
}	

void correctionTools::loadWeightParameters(){
	TFile weightFile((TString) _DATA + "/correctionFiles/"+ weight_name);
	std::cout<<"Opened weight file\n";
	for( int i = 0; i < 4; i++ ){
		for( int j = 0; j < 3; j++ ){
			
			TF1 * fp = (TF1*)weightFile.Get( Form("fitPipQ_f_%i_%i", i, j));
			weight_Parameters[0][i][j][0] = fp->GetParameter(0);
			weight_Parameters[0][i][j][1] = fp->GetParameter(1);
			weight_Parameters[0][i][j][2] = fp->GetParameter(2);
			weight_Parameters[0][i][j][3] = fp->GetParameter(3);
			
			TF1 * fm = (TF1*)weightFile.Get(Form("fitPimQ_f_%i_%i", i, j));
			weight_Parameters[1][i][j][0] = fm->GetParameter(0);
			weight_Parameters[1][i][j][1] = fm->GetParameter(1);
			weight_Parameters[1][i][j][2] = fm->GetParameter(2);
			weight_Parameters[1][i][j][3] = fm->GetParameter(3);
		}
	}
}

void correctionTools::loadkTopiParameters(){
	TFile weightFile((TString) _DATA + "/correctionFiles/"+ kaon_To_pi_name);
	std::cout<<"Opened p2k file\n";
	for( int i = 0; i < 4; i++ ){
		for( int j = 0; j < 4; j++ ){
			for( int k = 2; k < 4; k++ ){
				TF1 * fp = (TF1*)weightFile.Get(Form("fitPipX_%i_%i_%i", k, i, j));
				k_To_pi_Parameters[0][k][i][j][0] = fp->GetParameter(0);
				k_To_pi_Parameters[0][k][i][j][1] = fp->GetParameter(1);
				k_To_pi_Parameters[0][k][i][j][2] = fp->GetParameter(2);
				k_To_pi_Parameters[0][k][i][j][3] = fp->GetParameter(3);
				//k_To_pi_Parameters[0][k][i][j][4] = fp->GetParameter(4);
				
				TF1 * fm = (TF1*)weightFile.Get(Form("fitPimX_%i_%i_%i", k, i, j));
				k_To_pi_Parameters[1][k][i][j][0] = fm->GetParameter(0);
				k_To_pi_Parameters[1][k][i][j][1] = fm->GetParameter(1);
				k_To_pi_Parameters[1][k][i][j][2] = fm->GetParameter(2);
				k_To_pi_Parameters[1][k][i][j][3] = fm->GetParameter(3);
				//k_To_pi_Parameters[1][k][i][j][4] = fm->GetParameter(4);
				
				//std::cout<<"Filled Parameters : i = "<<i<<" , j = "<<j<<std::endl;
			}
		}
	}
}

void correctionTools::loadpiTokParameters(){
	TFile weightFile((TString) _DATA + "/correctionFiles/"+ pi_To_kaon_name);
	for( int i = 0; i < 4; i++ ){//momentum bins
		for( int j = 0; j < 4; j++ ){//loop through fits for each a(x, q)
			for( int k = 2; k < 4; k++ ){//loop through fits for each b(x)
				TF1 * fp = (TF1*)weightFile.Get(Form("fitPipX_%i_%i_%i", k, i, j));
				pi_To_k_Parameters[0][k][i][j][0] = fp->GetParameter(0);
				pi_To_k_Parameters[0][k][i][j][1] = fp->GetParameter(1);
				pi_To_k_Parameters[0][k][i][j][2] = fp->GetParameter(2);
				pi_To_k_Parameters[0][k][i][j][3] = fp->GetParameter(3);
				//pi_To_k_Parameters[0][k][i][j][4] = fp->GetParameter(4);
			
				TF1 * fm = (TF1*)weightFile.Get(Form("fitPimX_%i_%i_%i", k, i, j));
				pi_To_k_Parameters[1][k][i][j][0] = fm->GetParameter(0);
				pi_To_k_Parameters[1][k][i][j][1] = fm->GetParameter(1);
				pi_To_k_Parameters[1][k][i][j][2] = fm->GetParameter(2);
				pi_To_k_Parameters[1][k][i][j][3] = fm->GetParameter(3);
				//pi_To_k_Parameters[1][k][i][j][4] = fm->GetParameter(4);
			
			}
		}
	}
}

void correctionTools::setKinematics( double x, double q, double z, double p ){ 
	kin[0] = x;
	kin[1] = q;
	kin[2] = z;
	kin[3] = p;
}


double correctionTools::getCorrectionFactor( int type, int charge ){
	if( charge != 0 && charge != 1 ){
		std::cout<<"Invalid charge\n";
		return 0;
	}

	double x = kin[0];
	double q = kin[1];
	double z = kin[2];
	double p = kin[3];

	int this_bin_p = -1;
			
	for( int j= 0; j < bins_p; j++ ){
		if( p > p_bin_edges[j] && p <= p_bin_edges[j+1] ){
			this_bin_p = j;
		}
	}

	
	if(mode == 0){
		double weight = 0;
		double pi_To_k = 0;	
		double k_To_pi = 0;	
		
		for( int i = 0; i < 4; i++ ){
			double b = 0;
			for( int j = 0; j < 3; j++ ){
				double c = 0;
			
				for( int k = 0; k < 4; k++ ){

			
					c += weight_Parameters[charge][i][j][k]*pow( x, k );
				}
				
				b+= c*pow(q, j);
			}
			weight+= b*pow(z, i);
		}
				
		if( type == 1 ) return weight;

		if( type != 1 && this_bin_p < 2 ){
			std::cout<<"Entered momentum p = "<<p<<std::endl;
			std::cout<<"Returned bin : " <<this_bin_p<<std::endl;
			std::cout<<"Invalid momentum bin\n";
			return 1; 
		}

		for( int i = 0; i < 4; i++ ){
			double b = 0;
			for( int j = 0; j < 4; j++ ){
				double c = 0;
			
				for( int k = 0; k < 4; k++ ){
					c += pi_To_k_Parameters[charge][this_bin_p][i][j][k]*pow( x, k );
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
					c += k_To_pi_Parameters[charge][this_bin_p][i][j][k]*pow( x, k );
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

	else if(mode == 1){
		int this_bin_Q2 = (int)( ( (q - Q2_min)/(Q2_max-Q2_min) )*bins_Q2);
                int this_bin_xB = (int)( ( (x - xB_min)/(xB_max-xB_min) )*bins_xB);
		
		double mcWeight = weightFit[charge][this_bin_xB][this_bin_Q2]->Eval(z);
		double pi2kWeight = pi2kFit[charge][this_bin_xB][this_bin_Q2][this_bin_p]->Eval(z); 
		double k2piWeight = k2piFit[charge][this_bin_xB][this_bin_Q2][this_bin_p]->Eval(z); 
		double kWeight = 0;

		if( type == 0 ) return 1.;
		if( type == 1 ) return mcWeight;
		if( type == 2 ) kWeight = pi2kWeight;
		if( type == 3 ) kWeight = k2piWeight;
		if( type > 3 ){ return 0; }	
		
		if( type == 2  && this_bin_p < 2 ) return 1;
		if( type == 3  && this_bin_p < 2 ) return 0;
		//Check if in correction phase space       
		if( kWeight <= 0 ){return 0;}
		if( kWeight >= 1 ){return 1;}
		/*
		else{
			double kWeightTemp;
			int tempCharge;
			if( charge == 0 ){ tempCharge = 1; }
			if( charge == 1 ){ tempCharge = 0; }

			if(type == 2){
				kWeightTemp = pi2kWeight[tempCharge][this_bin_xB][this_bin_Q2][this_bin_p]->Eval(z); 
			}
			if( type == 3 ){
				kWeightTemp = k_to_pi_Correction[tempCharge][this_bin_p]->GetBinContent( x_bin_k, y_bin_k, z_bin_k );
			}

			if( kWeightTemp <= 0 || kWeightTemp > 1 ){
			       //if( this_bin_p >= 2 ){ return 1; }
			       //else{ return 0; }
				return 1;
			}
			return kWeight;
		}
		*/
		return 1;

	}
	else{
		int x_bin_w = binMigration[1]->GetXaxis()->FindBin(x);
		int y_bin_w = binMigration[1]->GetYaxis()->FindBin(q);
		int z_bin_w = binMigration[1]->GetZaxis()->FindBin(z);
		int x_bin_k = pi_to_k_Correction[charge][0]->GetXaxis()->FindBin(x);
		int y_bin_k = pi_to_k_Correction[charge][0]->GetYaxis()->FindBin(q);
		int z_bin_k = pi_to_k_Correction[charge][0]->GetZaxis()->FindBin(z);
		double binWeight = binMigration[charge]->GetBinContent( x_bin_w, y_bin_w, z_bin_w);
		double accWeight = accCorrection[charge]->GetBinContent( x_bin_w, y_bin_w, z_bin_w);
		double pi2kWeight = pi_to_k_Correction[charge][this_bin_p]->GetBinContent( x_bin_k, y_bin_k, z_bin_k );
		double k2piWeight = k_to_pi_Correction[charge][this_bin_p]->GetBinContent( x_bin_k, y_bin_k, z_bin_k );
		double kWeight = 0;

		if( type == 0 ) return binWeight;
		if( type == 1 ) return accWeight;
		if( type == 2 ) kWeight = pi2kWeight;
		if( type == 3 ) kWeight = k2piWeight;
		if( type > 3 ){ return 0; }	
		
		//Check if in correction phase space       
		if( kWeight <= 0 || kWeight > 1 ){return 1;}
		else{
			double kWeightTemp;
			int tempCharge;
			if( charge == 0 ){ tempCharge = 1; }
			if( charge == 1 ){ tempCharge = 0; }

			if(type == 2){
				kWeightTemp = pi_to_k_Correction[tempCharge][this_bin_p]->GetBinContent( x_bin_k, y_bin_k, z_bin_k );
			}
			if( type == 3 ){
				kWeightTemp = k_to_pi_Correction[tempCharge][this_bin_p]->GetBinContent( x_bin_k, y_bin_k, z_bin_k );
			}

			if( kWeightTemp <= 0 || kWeightTemp > 1 ){
			       //if( this_bin_p >= 2 ){ return 1; }
			       //else{ return 0; }
				return 1;
			}
			return kWeight;
		}
		
		return 1;

	}

}		

double correctionTools::getCorrectionError( int type, int charge ){
	if( charge != 0 && charge != 1 ){
		std::cout<<"Invalid charge\n";
		return 0;
	}

	double x = kin[0];
	double q = kin[1];
	double z = kin[2];
	double p = kin[3];

	int this_bin_p = -1;
			
	for( int j= 0; j < bins_p; j++ ){
		if( p > p_bin_edges[j] && p <= p_bin_edges[j+1] ){
			this_bin_p = j;
		}
	}
	int this_bin_Q2 = (int)( ( (q - Q2_min)/(Q2_max-Q2_min) )*bins_Q2);
        int this_bin_xB = (int)( ( (x - xB_min)/(xB_max-xB_min) )*bins_xB);

	if( mode == 0 ){ return 0; }
	if( mode == 1 ){ 

		double mcErr = 0;
		double pi2kErr = 0;
		double k2piErr = 0;

		if( type == 0 ) return 1.;
	
		for( int i = 0; i < 4; i++ ){
			mcErr += pow(pow(x, i)* weightFit[charge][this_bin_xB][this_bin_Q2]->GetParError(i), 2);
			pi2kErr += pow(pow(x, i)* pi2kFit[charge][this_bin_xB][this_bin_Q2][this_bin_p]->GetParError(i), 2);
			k2piErr += pow(pow(x, i)* k2piFit[charge][this_bin_xB][this_bin_Q2][this_bin_p]->GetParError(i), 2);
		}
		
		if( type == 1 ) return mcErr;
		if( type == 2  && this_bin_p < 2 ) return 1;
		if( type == 2 ) return pi2kErr;
		if( type == 3  && this_bin_p < 2 ) return 0;
		if( type == 3 ) return  k2piErr;
		if( type > 3 ){ return 0; }	
		
		//Check if in correction phase space       

	}
	else{
		int x_bin_w = binMigration[1]->GetXaxis()->FindBin(x);
		int y_bin_w = binMigration[1]->GetYaxis()->FindBin(q);
		int z_bin_w = binMigration[1]->GetZaxis()->FindBin(z);
		int x_bin_k = pi_to_k_Correction[charge][0]->GetXaxis()->FindBin(x);
		int y_bin_k = pi_to_k_Correction[charge][0]->GetYaxis()->FindBin(q);
		int z_bin_k = pi_to_k_Correction[charge][0]->GetZaxis()->FindBin(z);
		
		double binErr = binMigration[charge]->GetBinContent( x_bin_w, y_bin_w, z_bin_w);
		double accErr = accCorrection[charge]->GetBinContent( x_bin_w, y_bin_w, z_bin_w);
		double pi2kErr = pi_to_k_Correction[charge][this_bin_p]->GetBinContent( x_bin_k, y_bin_k, z_bin_k );
		double k2piErr = k_to_pi_Correction[charge][this_bin_p]->GetBinContent( x_bin_k, y_bin_k, z_bin_k );

		if( type == 0 ) return binErr;
		if( type == 1 ) return accErr;
		if( type == 2 ) return pi2kErr;
		if( type == 3 ) return k2piErr;
		if( type > 3 ){ return 0; }	
		
	}

	return 0;
}		

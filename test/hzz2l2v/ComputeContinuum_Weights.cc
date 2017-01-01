#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdio.h>

#include <TH1D.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TFile.h>
#include <TPaveStats.h>
#include <TPad.h>

using namespace std;

void ComputeContinuum_Weights(){


	string File="plotter_2016_12_30_ComputeWeights.root";
	string Dir="/storage_mnt/storage/user/amagitte/HZZ_2L2Nu_Analysis_2016/HiggsCoupling_2016/GreenLight_9Novembre2016/CMSSW_8_0_14/src/UserCode/llvv_fwk/test/hzz2l2v/";

	string Mass[13]={"200","300","400","500","600","700","800","900","1000","1500","2000","2500","3000"};
	string Cp[4]={"1.00","0.60","0.30","0.10"};
	string Channel[2]={"ggH","qqH"};

	FILE *outfile = fopen("Continuum_sF.txt","w");

	//Loop on the channel
	for(unsigned int ch=0; ch<2; ch++){
		//Loop on the Cprime
		for(unsigned int cp=0; cp<4; cp++){
			//Loop on the Mass
			for(unsigned int m=0; m<13; m++){

				double sF=0.0;				

				string NameHistoSig="all_higgsMass_shape_cp"+Cp[cp]+"_brn0.00";
				string NameHistoBckg="all_higgsMass_raw";	
			
				string DirSig=Channel[ch]+"("+Mass[m]+")_ContinuumOnly";
				string DirBckg="ZZ_filt1113";

				std::cout << DirSig << ", " << DirBckg << std::endl;
				std::cout << NameHistoSig << ", " << NameHistoBckg << std::endl;
				TFile *RFile=new TFile((Dir+File).c_str());
				if( RFile->GetDirectory(DirBckg.c_str())!=0 && RFile->GetDirectory(DirSig.c_str())!=0 ){ 
					RFile->cd(DirBckg.c_str());
					TH1F *HBckg=(TH1F*)gDirectory->Get(NameHistoBckg.c_str());
					RFile->cd(DirSig.c_str());
					TH1F *HSig=(TH1F*)gDirectory->Get(NameHistoSig.c_str());
	
					sF=HBckg->Integral(99,310)/HSig->Integral(99,310);
				} else { sF=1.0; }

				fprintf(outfile, "Channel: %s Cprime: %s Mass: %s sF: %20.19f \n", Channel[ch].c_str(), Cp[cp].c_str(), Mass[m].c_str(), sF);
			
			}
		}
	}

}

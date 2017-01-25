#include "UserCode/llvv_fwk/interface/HxswgUtils.h"
#include "TGraphErrors.h"

namespace Hxswg{

  namespace utils{

     TGraph* makeGraphFromColXandY(std::string dataFile, int colX, int colY){
        FILE* pFile = fopen(dataFile.c_str(), "r");
        if(!pFile){  printf("Couldn't open file %s to read Higgs width values\n", dataFile.c_str()); return NULL; }
        TGraph* toReturn = new TGraph(9999);  int N=0;
        char line [4096];
        while(fgets(line, 4096, pFile)){
           if(std::string(line).find("//")==0)continue; //skip line starting by //
           char* pch=strtok(line,"\t"); int Arg=0; double x; double y;
           while (pch!=NULL){ 
              if(Arg==colX){        sscanf(pch, "%lf", &x); 
              }else if(Arg==colY){  sscanf(pch, "%lf", &y);
              }              
              pch=strtok(NULL,"\t");Arg++; 
           }
           toReturn->SetPoint(N, x, y); N++; 
        }fclose(pFile);
        toReturn->Set(N);
        return toReturn;
   }

    TGraph* getXSec(std::string Name){
       if(Name.find("SM")!=std::string::npos){
          if(Name.find("VBF")!=std::string::npos){
	     //getHWidthExtended gives the BR values for all the mass points above 1000 GeV
             if(Name.find("13TeV")!=std::string::npos){return multiplyGraph( getVBFXSec13TeV(), getBRHtoZZ());}
             if(Name.find("8TeV" )!=std::string::npos){return multiplyGraph(  getVBFXSec8TeV(), getBRHtoZZ());}
             if(Name.find("7TeV" )!=std::string::npos){return multiplyGraph(  getVBFXSec7TeV(), getBRHtoZZ());}
             return NULL;
          }else{ //GGF
             if(Name.find("13TeV")!=std::string::npos){return multiplyGraph( getGGFXSec13TeV(), getBRHtoZZ());}
             if(Name.find("8TeV" )!=std::string::npos){return multiplyGraph(  getGGFXSec8TeV(), getBRHtoZZ());}
             if(Name.find("7TeV" )!=std::string::npos){return multiplyGraph(  getGGFXSec7TeV(), getBRHtoZZ());}
             return NULL;
          }
       }else if(Name.find("RsGrav")!=std::string::npos){
          if(Name.find("13TeV" )!=std::string::npos){return multiplyGraph(   getRsGravXSec13TeV(),  getBRRsGravtoZZ());}
          return NULL;
       }else if(Name.find("BulkGrav")!=std::string::npos){
          if(Name.find("13TeV" )!=std::string::npos){return multiplyGraph( getBulkGravXSec13TeV(), getBRBulkGravtoZZ());}
          return NULL;
       }else if(Name.find("Rad")!=std::string::npos){
          if(Name.find("13TeV" )!=std::string::npos){return multiplyGraph( getRadXSec13TeV(), getBRRadtoZZ());}
          return NULL;
       }
     
       return NULL;
    }

    TGraph* getVBFoverGGF(std::string Name){
       if(Name.find("13TeV")!=std::string::npos){return divideGraph(getVBFXSec13TeV(), getGGFXSec13TeV());}
       if(Name.find("8TeV" )!=std::string::npos){return divideGraph(getVBFXSec8TeV(),  getGGFXSec8TeV() );}
       if(Name.find("7TeV" )!=std::string::npos){return divideGraph(getVBFXSec7TeV(),  getGGFXSec7TeV() );}
       return NULL;
    }

    TGraph *getXSecMELA( float cprime){

	double    mass[13] = {200,300,400,500,600,700,800,900,1000,1500,2000,2500,3000};
	double XSec100[13] = {19.167002,77.414873,56.635416,32.492678,12.572093,5.528836,2.265198,1.107643,0.539356,0.087179,0.027123,0.009779,0.007671};
	double XSec060[13] = {52.497654,205.635846,165.082019,98.609912,38.749277,17.467784,7.370453,3.786566,1.891071,0.338315,0.095826,0.030301,0.007671};
	double XSec030[13] = {210.000136,807.658034,682.012702,409.669810,160.025744,73.459488,31.378913,16.290753,8.081482,1.475269,0.367997,0.096478,0.007671};
	double XSec010[13] = {1928.504989,7136.708697,6186.640072,3880.054859,1429.841422,650.550864,285.01729,148.040109,73.428032165,13.407377,3.14365,0.77877,0.007670};	
	double XSec[13];

	std::cout << "Cprime: " << cprime << std::endl;

	if( Hxswg::utils::Equal(cprime, 1.0, 0.000001) ){
		std::cout << "Inside C'=1.0" << std::endl;
		for(unsigned int k=0;k<13;k++){ XSec[k]=XSec100[k]; }
	}else if( Hxswg::utils::Equal(cprime, 0.6, 0.000001) ){ 
		std::cout << "Inside C'=0.6" << std::endl;
		for(unsigned int k=0;k<13;k++){ XSec[k]=XSec060[k]; }
	}else if( Hxswg::utils::Equal(cprime, 0.3, 0.000001) ){
		std::cout << "Inside C'=0.3" << std::endl;
		for(unsigned int k=0;k<13;k++){ XSec[k]=XSec030[k]; }
	}else if( Hxswg::utils::Equal(cprime, 0.1, 0.000001) ){
		std::cout << "Inside C'=0.1" << std::endl;
		for(unsigned int k=0;k<13;k++){ XSec[k]=XSec010[k]; }
	}

	TGraph *XSecMELA= new TGraph( 13, mass, XSec);

	return XSecMELA;
    }


  }
}

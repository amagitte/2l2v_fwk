#include "UserCode/llvv_fwk/interface/HiggsUtils.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "UserCode/llvv_fwk/interface/th1fmorph.h"

namespace higgs{

  namespace utils{

    //
    EventCategory::EventCategory(int mode):mode_(mode) { }

    EventCategory::~EventCategory() { }
    
    TString EventCategory::GetCategory(pat::JetCollection &jets, LorentzVector &boson)
    {
      //jet multiplicity
      int NJets(0);
      for(size_t ijet=0; ijet<jets.size(); ijet++){
	if(jets[ijet].pt()<=30)continue;
	NJets++;
      }
      
      //VBF tag
      bool isVBF(false);
      if(NJets>=2){
	LorentzVector VBFSyst = jets[0].p4() + jets[1].p4();
	double j1eta=jets[0].eta() ;
	double j2eta=jets[1].eta();
	double dEta = fabs(j1eta-j2eta);
	
	int NCentralJet(0), NCentralBoson(0);
	double MaxEta, MinEta;
	if(j1eta<j2eta) { MinEta=j1eta; MaxEta=j2eta;}
	else            { MinEta=j2eta; MaxEta=j1eta;}
	for(size_t ijet=2; ijet<jets.size(); ijet++){
	  float jpt=jets[ijet].pt();
	  float jeta=jets[ijet].eta();
	  if(jpt<30)continue; 
	  if(jeta>MinEta && jeta<MaxEta) NCentralJet++;  
	}
	
	if(boson.eta()>MinEta && boson.eta()<MaxEta) NCentralBoson=1;
	isVBF=( (dEta>4.0) && (VBFSyst.M()>500) && (NCentralJet==0) && (NCentralBoson==1) );
      }
      
      //build classification
      TString cat("");
      switch(mode_)
	{
	case EXCLUSIVEVBF:
	  {
	    cat= isVBF ? "vbf":"novbf";
	    break;
	  }
	case EXCLUSIVE3JETS:
	  {
	    if(NJets==0)      cat="eq0jets";
	    else if(NJets==1) cat="eq1jets";
	    else              cat="geq2jets";
	    break;
	  }
	case EXCLUSIVE3JETSVBF:
	  {
	    if(isVBF)         cat="vbf";
	    else if(NJets==0) cat="eq0jets";
	    else if(NJets==1) cat="eq1jets";
	    else              cat="geq2jets";
	    break;
	  }
	case EXCLUSIVE2JETS:
	  {
	    if(NJets==0)      cat="eq0jets";
	    else              cat="geq1jets";
	    break;
	  }
	case  EXCLUSIVE2JETSVBF:
	  {
	    if(isVBF)         cat="vbf";
	    else if(NJets==0) cat="eq0jets";
	    else              cat="geq1jets";
	    break;
	  }
	default:
	  break;
	}
	
      return cat;
    }


    
    //
    double transverseMass(const LorentzVector &visible, const LorentzVector &invisible, bool assumeSameMass){
      if(assumeSameMass){
	LorentzVector sum=visible+invisible;
	double tMass = pow(sqrt(pow(visible.pt(),2)+pow(visible.mass(),2))+sqrt(pow(invisible.pt(),2)+pow(91.188,2)),2);
	tMass-=pow(sum.pt(),2);
	return sqrt(tMass);
      }else{
	double dphi=fabs(deltaPhi(invisible.phi(),visible.phi()));
	return sqrt(2*invisible.pt()*visible.pt()*(1-cos(dphi)));
      }
      return -1;
    }
    
    //    
    double weightToH125Interference(double mass,double width,TFile *intFile, TString var)
    {
      if(width==0 || intFile==0) return 1.0;
      TString name("weights_ceq"); name+=Int_t(width); 
      if(var!="") name += "_"+var;
      TGraphErrors *gr=(TGraphErrors *)intFile->Get(name);
      if(gr==0) return 1.0;
      return gr->Eval(mass);
    }


    TGraph* getWeightGraphFromShapes(TH1D* newLineShape, TH1D* originalLineShape, double mH){
      double RMS = originalLineShape->GetRMS();
      int RebinFactor = 1; while((originalLineShape->GetBinWidth(1)*RebinFactor) < RMS/ 8.0){RebinFactor*=2;}

      TH1D* weightsHDenominator = (TH1D*)originalLineShape->Clone("weightsHDenominator");    weightsHDenominator->Rebin(RebinFactor);
//    TH1D* weightsH   = mH>=400?(TH1D*)hSI_nnlo  ->Clone("weightsH")  : (TH1D*)hS_nnlo->Clone("weightsH");   weightsH  ->Rebin(RebinFactor);  weightsH->  Divide(weightsHDenominator);
      TH1D* weightsH   = (TH1D*)newLineShape->Clone("weightsH");   weightsH  ->Rebin(RebinFactor);  weightsH->  Divide(weightsHDenominator);
      for(int x=1;x<=weightsH->GetNbinsX();x++){
         weightsH  ->SetBinContent(x, std::max(0.0,weightsH  ->GetBinContent(x))  );   weightsH  ->SetBinError(x, 0.0  ); //truncate to positive values
         if(fabs(weightsH->GetXaxis()->GetBinCenter(x) - mH)>3*RMS){ weightsH  ->SetBinContent(x, 0.0  );   weightsH  ->SetBinError(x, 0.0  );} //truncate to 3sigma
      }
      TGraph* toReturn=new TGraph(weightsH);
      delete weightsH; 
      delete weightsHDenominator;
      return toReturn;     
   }

    //    
    TH1D* getHistoFromNRfile(std::string histoName, double mass, double Cprime, double BRnew, TFile *nrLineShapesFile){
      //We only have lineshape for Cprime=X and BRnew=0, so we need to find the lineshape equivalent to (Cprime, BRnew) to (Csecond, 0)
      //that is easy because the signal width for (Cprime, Brnew) is SM width * cprimeÂ / (1-BRnew), so the equivalent with can be taken from the pair (Cprime/sqrt(1-BRnew), 0) in order to get the same width
      //some care is needed for the cross-section because it does not scale the same way, but this does not matter since we have the xsection normalized to SM afterward
      //do not allow Csecond to be larger than 1 though

      char nrShapeBuf[100];


      double Csecond = Cprime;
      double BRnew2 = BRnew;
      if(BRnew!=0){
         Csecond = BRnew<=0?Cprime:std::min(1.0, Cprime / sqrt(1-BRnew));
         BRnew2  = 0;        
         printf("BRnew is different than 0, so apply the change of variable (cprime=%f, Brnew=%f) --> (csecond=%f, 0)\n", Cprime, BRnew, Csecond);
      }

      //1st check if is in the file
      if(!nrLineShapesFile){
         printf("LineShapeFile not ok for narrow resonnance\n");  fflush(stdout);
         return NULL;
      }

      //very likely Csecond is not a multiple of 0.1 anymore,
      //we need to morph between to 0.1 multiple
      //first check if Csecond is multiple of 0.1  (in steps of 0.01)
      if( (int(Csecond*10000)/10)%100==0 ){  
  	     sprintf(nrShapeBuf,"%d_%3.1f_%3.1f/%s",int(mass),Csecond,BRnew2, histoName.c_str());
             return (TH1D*)nrLineShapesFile->Get(nrShapeBuf);

     //Csecond is NOT a multiple of 0.1 --> need to morph
      }else{ 
         printf("Csecond %f is NOT a multiple of 0.1 --> %i --> %i\n", Csecond, int(Csecond*1000), int(Csecond*1000)%100);

         //identify neighboring values that are multiple of 0.1
         double CsecondL = (int(Csecond*1000)/100)/10.0;
         double CsecondR = CsecondL+0.1;
         printf("morph the lineshape between %f and %f\n", CsecondL, CsecondR);
         TH1D* HL=getHistoFromNRfile(histoName, mass, CsecondL, BRnew2, nrLineShapesFile);
         TH1D* HR=getHistoFromNRfile(histoName, mass, CsecondR, BRnew2, nrLineShapesFile);
         if(!HL || !HR){printf("Left and Right NR shapes can not be found in file\n"); return NULL;}          
         return th1fmorph("hC","hC", HL, HR, CsecondL, CsecondR, Csecond, 1.0, 0);
      }

      return NULL;
  }


    //    
    TGraph* weightNarrowResonnance(bool isVBF, double mass, double Cprime, double BRnew, TFile *nrLineShapesFile, double& Norm, TString pf){
      if((Cprime<0 || BRnew<0) || (Cprime==0 && BRnew==0)){
         TGraph* g = new TGraph(2);
         g->SetPoint(0,    0, 1.0);
         g->SetPoint(1, 9999, 1.0);
         return g;
      }


      //We only have lineshape for Cprime=X and BRnew=0, so we need to find the lineshape equivalent to (Cprime, BRnew) to (Csecond, 0)
      //that is easy because the signal width for (Cprime, Brnew) is SM width * cprimeÂ / (1-BRnew), so the equivalent with can be taken from the pair (Cprime/sqrt(1-BRnew), 0) in order to get the same width
      //some care is needed for the cross-section because it does not scale the same way, but this does not matter since we have the xsection normalized to SM afterward
      //do not allow Csecond to be larger than 1 though
      double Csecond = BRnew<=0?Cprime:std::min(1.0, Cprime / sqrt(1-BRnew));
      double BRnew2  = 0;

      if(BRnew!=0)printf("BRnew is different than 0, so apply the change of variable (cprime=%f, Brnew=%f) --> (csecond=%f, 0)\n", Cprime, BRnew, Csecond);

      //1st check if is in the file
      if(!nrLineShapesFile){
         printf("LineShapeFile not ok for narrow resonnance\n");  fflush(stdout);
         return NULL;
      }

      TH1D* original  = getHistoFromNRfile("mH_S_NR"                  , mass, 1.0, 0.0, nrLineShapesFile);  // SM signal only
//      TH1D* lineshape = getHistoFromNRfile(TString("mH_SI_NR_nnlo")+pf, mass, Csecond, BRnew2, nrLineShapesFile);  // (signal + interference) * NNLOKFactors   
      TH1D* lineshape = NULL;
      if(isVBF){
         //Signal+Interference[h2-h1, h2-Bckg, h1-Bckg] (LO)
         lineshape = getHistoFromNRfile((TString("mH_SI_NR")+pf).Data(), mass, Csecond, BRnew2, nrLineShapesFile);              
      } else {
	 //Signal+Interference[h2-h1, h2-Bckg, h1-Bckg] (LO*KFactor NNLO)
         lineshape = getHistoFromNRfile((TString("mH_SI_NR_nnlo")+pf).Data(), mass, Csecond, BRnew2, nrLineShapesFile);     
      }
      if(!original || !lineshape)return NULL;

      Norm = lineshape->Integral()/original->Integral();

      TGraph* nrGr = getWeightGraphFromShapes(lineshape, original, mass);
      for(int i=0;i<nrGr->GetN();i++){double x, y; nrGr->GetPoint(i, x, y);  if(y<0 || y>1000){y=0;} nrGr->SetPoint(i, x, y);} //make sure weights are not crazy

      //add 20% uncertainty on VBF lineshape
      if(isVBF){
         if(pf.Contains("up"  )){        for(int i=0;i<nrGr->GetN();i++){double x, y; nrGr->GetPoint(i, x, y);  nrGr->SetPoint(i, x, y*1.2);}      }
         if(pf.Contains("down")){        for(int i=0;i<nrGr->GetN();i++){double x, y; nrGr->GetPoint(i, x, y);  nrGr->SetPoint(i, x, y*0.8);}      }
      }


      return nrGr;
   }


    TGraph* weightGGZZContinuum(TFile *nrLineShapesFile, double& Norm, TString pf){
      //1st check if is in the file
      if(!nrLineShapesFile){
         printf("LineShapeFile not ok for ggZZ Continuum\n");  fflush(stdout);
         return NULL;
      }

      char nrShapeBuf[100];
      sprintf(nrShapeBuf,"%d_%3.1f_%3.1f/%s%s",int(1000) ,1.0, 0.0, "kFactors", pf.Data());
      TGraph* nrGr = (TGraph*)nrLineShapesFile->Get(nrShapeBuf);
      for(int i=0;i<nrGr->GetN();i++){double x, y; nrGr->GetPoint(i, x, y);  if(y<0 || y>1000){y=0;} nrGr->SetPoint(i, x, y);} //make sure weights are not crazy

      Norm = 1.0;

      return nrGr;
   }

   double Get_NNLO_kFactors( double mass){

     double sF=1;
     TString nnlosf_FileUrl(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
     gSystem->ExpandPathName(nnlosf_FileUrl);
     TFile *nnlosf_File = TFile::Open(nnlosf_FileUrl);
     if(!nnlosf_File){
         printf("nnlo k-Factors for Signal is are not found\n");  fflush(stdout);
     }
     TGraph *sFGr = (TGraph*) nnlosf_File->Get("kfactor_Nominal");
     sF = sFGr->Eval(mass);
     return sF;

   }

   float ComputeInterfWeight( Mela& mela, bool isVBF, TString MelaMode, double heavyMass, double heavyWidth, SimpleParticleCollection_t& daughters, SimpleParticleCollection_t& associated, SimpleParticleCollection_t& mothers){

	float Interf_weight=0;

	if(isVBF){
		mela.setInputEvent(&daughters, &associated, &mothers, true);
                if( MelaMode.Contains("Interf") && MelaMode.Contains("h2") && MelaMode.Contains("Continuum") ){

			float Bckg_wg=0; float Sigh2_wg=0; float All_wg=0;

			//Bckg Only
                	mela.setProcess( TVar::bkgZZ, TVar::MCFM, TVar::JJVBF_S);
               		mela.computeProdDecP( Bckg_wg, false);

			//Heavy Boson Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);	
			mela.computeProdDecP( Sigh2_wg, false);	


			//Bckg+Heavy Boson
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);	
                        mela.computeProdDecP( All_wg, false); 

			Interf_weight = All_wg - Sigh2_wg - Bckg_wg;

		} else if( MelaMode.Contains("Interf") && MelaMode.Contains("h1") && MelaMode.Contains("Continuum") ){

                        float Bckg_wg=0; float Sigh1_wg=0; float All_wg=0;

			//Bckg Only
                        mela.setProcess( TVar::bkgZZ, TVar::MCFM, TVar::JJVBF_S);
                        mela.computeProdDecP( Bckg_wg, false); 

                        //Light Higgs Only 
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.computeProdDecP( Sigh1_wg, false); 

                        //Bckg+Light Higgs 
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.computeProdDecP( All_wg, false); 
                        
                        Interf_weight = All_wg - Sigh1_wg - Bckg_wg;
                
		} else if( MelaMode.Contains("Interf") && MelaMode.Contains("Full") ){

                        float Bckg_wg=0; float Sigh1_wg=0; float Sigh2_wg=0; float All_wg=0;

			//Bckg Only
                        mela.setProcess( TVar::bkgZZ, TVar::MCFM, TVar::JJVBF_S);
                        mela.computeProdDecP( Bckg_wg, false); 

			//Light Higgs Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.computeProdDecP( Sigh1_wg, false); 

			//Heavy Higgs Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
                        mela.computeProdDecP( Sigh2_wg, false); 

			//Bckg+Heavy Higgs+Light Higgs
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 1);
			mela.computeProdDecP( All_wg, false);

			Interf_weight = All_wg - Sigh1_wg - Sigh2_wg - Bckg_wg;

                } else if( MelaMode.Contains("Interf") && MelaMode.Contains("h2") && MelaMode.Contains("h1") ){

                        float Bckg_wg=0; float Sigh1_wg=0; float Sigh2_wg=0; float All_wg=0;
                        float All_h2Bckg=0; float All_h1Bckg=0; float h2Bckg=0; float h1Bckg=0;

			//Bckg Only
                        mela.setProcess( TVar::bkgZZ, TVar::MCFM, TVar::JJVBF_S);
                        mela.computeProdDecP( Bckg_wg, false);

			//Light Higgs Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.computeProdDecP( Sigh1_wg, false);

			//Heavy Higgs Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
                        mela.computeProdDecP( Sigh2_wg, false);

			//Bckg+Light Higgs+Heavy Higgs
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);	
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 1); 
                        mela.computeProdDecP( All_wg, false);

			//Bckg+Light Higgs
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.computeProdDecP( All_h1Bckg, false);

			//Bckg+Heavy Higgs
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
                        mela.computeProdDecP( All_h2Bckg, false);

                        h2Bckg = All_h2Bckg - Sigh2_wg - Bckg_wg;
                        h1Bckg = All_h1Bckg - Sigh1_wg - Bckg_wg;
                        Interf_weight = All_wg - Sigh1_wg - Sigh2_wg - Bckg_wg - h2Bckg - h1Bckg;

                } else if(MelaMode.Contains("Interf") && MelaMode.Contains("h2") && MelaMode.Contains("h1") && MelaMode.Contains("Continuum")){

			float All_wg=0; float All_h1Bckg=0; float Sigh2_wg=0;	
			//Bckg+Light Higgs
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.computeProdDecP( All_h1Bckg, false);

			//Heavy Higgs Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
                        mela.computeProdDecP( Sigh2_wg, false);

			//Light Higgs+Heavy Higgs+Continuum
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 1);
                        mela.computeProdDecP( All_wg, false);

			Interf_weight = All_wg - All_h1Bckg - Sigh2_wg;
		}

	} else {
                mela.setInputEvent(&daughters, 0, &mothers, true);
                if( MelaMode.Contains("Interf") && MelaMode.Contains("h2") && MelaMode.Contains("Continuum") ){

                        float Bckg_wg=0; float Sigh2_wg=0; float All_wg=0;
                       
			//Bckg Only 
                        mela.setProcess( TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
                        mela.computeP( Bckg_wg, false);
                       
			//Heavy Higgs Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
                        //mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
			mela.computeP( Sigh2_wg, false);
			  
			//Bckg+Heavy Higgs
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
                        //mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        mela.computeP( All_wg, false); 

                        Interf_weight = All_wg - Sigh2_wg - Bckg_wg;

                } else if( MelaMode.Contains("Interf") && MelaMode.Contains("h1") && MelaMode.Contains("Continuum") ){

                        float Bckg_wg=0; float Sigh1_wg=0; float All_wg=0;
                       
			//Bckg Only 
                        mela.setProcess( TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
                        mela.computeP( Bckg_wg, false);
                       
			//Light Higgs Only 
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        //mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        mela.computeP( Sigh1_wg, false);
                       
			//Bckg+Light Higgs Only 
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        //mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        mela.computeP( All_wg, false); 
                        
                        Interf_weight = All_wg - Sigh1_wg - Bckg_wg;

                } else if( MelaMode.Contains("Interf") && MelaMode.Contains("Full") ){

                        float Bckg_wg=0; float Sigh1_wg=0; float Sigh2_wg=0; float All_wg=0;

			//Bckg Only
                        mela.setProcess( TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
                        mela.computeP( Bckg_wg, false);

			//Light Higgs Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        //mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        mela.computeP( Sigh1_wg, false);

			//Heavy Higgs Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
                        //mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        mela.computeP( Sigh2_wg, false);

			//Bckg+Light Higgs+Heavy Higgs
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 1);
                        //mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHbbcoupl[1][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHttcoupl[1][gHIGGS_KAPPA][0]=1.0;
                        mela.computeP( All_wg, false);

                        Interf_weight = All_wg - Sigh1_wg - Sigh2_wg - Bckg_wg;

                } else if( MelaMode.Contains("Interf") && MelaMode.Contains("h2") && MelaMode.Contains("h1") ){

                        float Bckg_wg=0; float Sigh1_wg=0; float Sigh2_wg=0; float All_wg=0; 
			float All_h2Bckg=0; float All_h1Bckg=0; float h2Bckg=0; float h1Bckg=0;

			//Bckg Only
                        mela.setProcess( TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
                        mela.computeP( Bckg_wg, false);

			//Light Higgs Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        //mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        mela.computeP( Sigh1_wg, false); 

			//Heavy Higgs Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
                        //mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        mela.computeP( Sigh2_wg, false); 

			//Bckg+Light Higgs+Heavy Higgs
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 1);
                        //mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHbbcoupl[1][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHttcoupl[1][gHIGGS_KAPPA][0]=1.0;
                        mela.computeP( All_wg, false);

			//Bckg+Light Higgs
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        //mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        mela.computeP( All_h1Bckg, false);

			//Bckg+Heavy Higgs
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
                        //mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        mela.computeP( All_h2Bckg, false);

			h2Bckg = All_h2Bckg - Sigh2_wg - Bckg_wg;
			h1Bckg = All_h1Bckg - Sigh1_wg - Bckg_wg;
                        Interf_weight = All_wg - Sigh1_wg - Sigh2_wg - Bckg_wg - h2Bckg - h1Bckg;

		} else if(MelaMode.Contains("Interf") && MelaMode.Contains("h2") && MelaMode.Contains("h1") && MelaMode.Contains("Continuum")){

                        float All_wg=0; float All_h1Bckg=0; float Sigh2_wg=0;

			//Light higgs+Continuum
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG); 
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        //mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        mela.computeP( All_h1Bckg, false);

			//Heavy Higgs
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
                        //mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        mela.computeP( Sigh2_wg, false);

			//Light Higgs+Heavy Higgs+Continuum
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
                        //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 1);
                        //mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHbbcoupl[1][gHIGGS_KAPPA][0]=1.0;
                        //mela.selfDHttcoupl[1][gHIGGS_KAPPA][0]=1.0;
                        mela.computeP( All_wg, false);

                        Interf_weight = All_wg - All_h1Bckg - Sigh2_wg;

		} 

	}

	return Interf_weight;

   }

   float weightNarrowResonnance_MELA( Mela& mela, bool isVBF, TString MelaMode, double CP, double heavyMass, fwlite::Event& eV){

        //Mela mela( 13, heavyMass, TVar::DEBUG);

        fwlite::Handle< LHEEventProduct > lheEv;
        lheEv.getByLabel(eV, "externalLHEProducer");

        //Weight to reweight the MELA shape to the real cross-section
        double continuum_weight=1;
        continuum_weight = weightContinuum_MELA(isVBF,CP,heavyMass);
	
	//Fill a Map with Mass and Width SM Like
	double heavyWidth=0; float weightSM=1; float weightMELA=1; float finalweight=1; //float kF; 

	std::map< double, double>  SM_Info;
	SM_Info[200]=1.43;  SM_Info[300]=8.43;  SM_Info[400]=29.3; 
	SM_Info[500]=68;    SM_Info[600]=123;   SM_Info[700]=199;
	SM_Info[800]=304;   SM_Info[900]=499;   SM_Info[1000]=647;
	SM_Info[1500]=1500; SM_Info[2000]=2000; SM_Info[2500]=2500;

	heavyWidth=SM_Info[heavyMass];

        SimpleParticleCollection_t daughters, mothers, associated; // associated;
        TLorentzVector Higgs;

	//Loop on particles and fill SimpleParticleCollection_t 
        for(int k=0; k<lheEv->hepeup().NUP; k++){

	    double PdgId=0.; double Status=0.;
	    PdgId=lheEv->hepeup().IDUP.at(k); 
	    Status=lheEv->hepeup().ISTUP.at(k);
	    double Px=lheEv->hepeup().PUP.at(k)[0]; double Py=lheEv->hepeup().PUP.at(k)[1]; 
	    double Pz=lheEv->hepeup().PUP.at(k)[2]; double  E=lheEv->hepeup().PUP.at(k)[3];
	    double Pt=std::sqrt(std::pow(Px,2)+std::pow(Py,2));
	    printf("Particle: %4.1f Status: %4.1f Pt: %10.5f \n", PdgId, Status, Pt);
            if( (abs(PdgId)<7.0 || PdgId==21.0) && Status<0.0 ){
                TLorentzVector partons( Px, Py, Pz, E);
		if (abs(PdgId)<7.0 && isVBF) mothers.push_back( SimpleParticle_t( PdgId, partons)); //Filling Infos
                else mothers.push_back(SimpleParticle_t(0, partons)); //Else fill gluons as 0 (unknown parton) in case the initial state is qg in ggF, or qg or gg in VBF
	    } else if ( (abs(PdgId)<7.0 || PdgId==21.0) && Status>0.0){
                TLorentzVector extra_partons( Px, Py, Pz, E);
                if (abs(PdgId)<7.0 && isVBF) associated.push_back( SimpleParticle_t( PdgId, extra_partons));
                else if(abs(PdgId)==21.0 && isVBF) associated.push_back(SimpleParticle_t(0, extra_partons)); 
	    } else if ( abs(PdgId)==11.0 || abs(PdgId)==12.0 || abs(PdgId)==13.0 || abs(PdgId)==14.0 || abs(PdgId)==15.0 || abs(PdgId)==16.0 ){
		TLorentzVector lepP( Px, Py, Pz, E);
		daughters.push_back( SimpleParticle_t( PdgId, lepP)); //Filling Infos
	    } else if (  abs(PdgId)==25.0 ){
                Higgs.SetPxPyPzE( Px, Py, Pz, E);
	    }

	}

	for(unsigned int m=0; m<mothers.size(); m++){
		printf("Mother collection=> Particle: %10i Pt: %10.5f \n", mothers.at(m).first, mothers.at(m).second.Pt());
	}
	printf(" \n");
        for(unsigned int m=0; m<associated.size(); m++){
                printf("Associated collection=> Particle: %10i Pt: %10.5f \n", associated.at(m).first, associated.at(m).second.Pt());
        }
        printf(" \n");
	std::sort( associated.begin(), associated.end(), utils::sort_CandidatesByPt_V2);
        for(unsigned int m=0; m<associated.size(); m++){
                printf("Associated collection post Soterd function=> Particle: %10i Pt: %10.5f \n", associated.at(m).first, associated.at(m).second.Pt());
        }
        printf(" \n");	

        mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ); //Mela Candidate mode initialized
	if(isVBF){ 
		mela.setInputEvent(&daughters, &associated, &mothers, true);
		mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S);
        	mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
		mela.computeProdDecP( weightSM, false);
	}else{ 
		mela.setInputEvent(&daughters, 0, &mothers, true);
		mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
        	mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
                mela.computeP( weightSM, false);
	} 
	
        //BSM reweighiting
        mela.resetInputEvent();

	heavyWidth=heavyWidth*CP*CP;

	if( !MelaMode.Contains("Interf") ){
        	if(isVBF){
                	mela.setInputEvent(&daughters, &associated, &mothers, true);
                	if(MelaMode.Contains("Continuum")){	
				mela.setProcess( TVar::bkgZZ, TVar::MCFM, TVar::JJVBF_S);
				mela.computeProdDecP( weightSM, false);
			} else if(MelaMode.Contains("Bckg")){
                                mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF_S);
                                //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                                mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
				mela.computeProdDecP( weightSM, false);
			} else if(MelaMode.Contains("Sigh2")){
                        	mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S);
        			//mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
        			mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
			        mela.computeProdDecP( weightSM, false);
			} else if(MelaMode.Contains("Sigh1")){
                        	mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        	//mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        	mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);	
	                	mela.computeProdDecP( weightSM, false);
			} else if(MelaMode.Contains("All")){
                        	mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        	//mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        	mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        	mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 1);
	                	mela.computeProdDecP( weightSM, false);
			}
        	}else{
                	mela.setInputEvent(&daughters, 0, &mothers, true);	
			if(MelaMode.Contains("Continuum")){
                		mela.setProcess( TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
			} else if(MelaMode.Contains("Bckg")){
                                mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
                                //mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                                mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
				//mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
				//mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0; 
			} else if(MelaMode.Contains("Sigh2")){ 
                        	mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
                        	//mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses	
                        	mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
                                //mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
                                //mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
			} else if(MelaMode.Contains("Sigh1")){ 
                        	mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
                        	//mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        	mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                                //mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
                                //mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
			} else if(MelaMode.Contains("All")){	
                        	mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
                        	//mela.setMelaHiggsMassWidth( -1, 0, 0); //Clean all the masses
                        	mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        	mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 1);
                                //mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
                                //mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
                                //mela.selfDHbbcoupl[1][gHIGGS_KAPPA][0]=1.0;
                                //mela.selfDHttcoupl[1][gHIGGS_KAPPA][0]=1.0;
                	}
        	}
		
		mela.computeP( weightMELA, false);

	} else if( MelaMode.Contains("Interf") ){
		weightMELA = ComputeInterfWeight( mela, isVBF, MelaMode, heavyMass, heavyWidth, daughters, associated, mothers);
	}

        //if(!isVBF) kF=Get_NNLO_kFactors(Higgs.M());

        //printf("Mela Weight BSM: %20.18f Mela Weigth SM: %20.18f Continuum Weigth: %20.18f \n", weightMELA, weightSM, continuum_weight);
	mela.resetInputEvent();
        finalweight=(weightMELA/weightSM)*continuum_weight;

        return finalweight;

    }

 float weightContinuum_MELA( bool isVBF, double CP, double heavyMass){

	double continuumWeight=1;
	std::map< double, std::map<double,double> > cpWeight_MapVBF;
        std::map< double, std::map<double,double> > cpWeight_MapggH;

	//Filling ContinuumWeights GGH
	cpWeight_MapggH[1.0][200] = 1.0; //0.0174916426606125; 
	cpWeight_MapggH[1.0][300] = 1.0; //0.151788185833618; 
	cpWeight_MapggH[1.0][400] = 1.0; //0.133535059751491; 
	cpWeight_MapggH[1.0][500] = 1.0; //0.0635195783662186;	
	cpWeight_MapggH[1.0][600] = 1.0; //0.0294003411175982;
	cpWeight_MapggH[1.0][700] = 1.0; //0.0134870105407229;
	cpWeight_MapggH[1.0][800] = 1.0; //0.00714088114955624;
	cpWeight_MapggH[1.0][900] = 1.0; //0.00310337236376096;	
	cpWeight_MapggH[1.0][1000] = 1.0; //0.00149478972041437;
	cpWeight_MapggH[1.0][1500] = 1.0; //0.00020851824925088;
	cpWeight_MapggH[1.0][2000] = 1.0;
	cpWeight_MapggH[1.0][2500] = 1.0;
	cpWeight_MapggH[1.0][3000] = 1.0;

        cpWeight_MapggH[0.6][200] = 1.0;
        cpWeight_MapggH[0.6][300] = 1.0;
        cpWeight_MapggH[0.6][400] = 1.0;
        cpWeight_MapggH[0.6][500] = 1.0;
        cpWeight_MapggH[0.6][600] = 1.0;
        cpWeight_MapggH[0.6][700] = 1.0;
        cpWeight_MapggH[0.6][800] = 1.0;
        cpWeight_MapggH[0.6][900] = 1.0;
        cpWeight_MapggH[0.6][1000] = 1.0;
        cpWeight_MapggH[0.6][1500] = 1.0;
        cpWeight_MapggH[0.6][2000] = 1.0;
        cpWeight_MapggH[0.6][2500] = 1.0;
        cpWeight_MapggH[0.6][3000] = 1.0;

        cpWeight_MapggH[0.3][200] = 1.0;
        cpWeight_MapggH[0.3][300] = 1.0;
        cpWeight_MapggH[0.3][400] = 1.0;
        cpWeight_MapggH[0.3][500] = 1.0;
        cpWeight_MapggH[0.3][600] = 1.0;
        cpWeight_MapggH[0.3][700] = 1.0;
        cpWeight_MapggH[0.3][800] = 1.0;
        cpWeight_MapggH[0.3][900] = 1.0;
        cpWeight_MapggH[0.3][1000] = 1.0;
        cpWeight_MapggH[0.3][1500] = 1.0;
        cpWeight_MapggH[0.3][2000] = 1.0;
        cpWeight_MapggH[0.3][2500] = 1.0;
        cpWeight_MapggH[0.3][3000] = 1.0;

        cpWeight_MapggH[0.1][200] = 1.0;
        cpWeight_MapggH[0.1][300] = 1.0;
        cpWeight_MapggH[0.1][400] = 1.0;
        cpWeight_MapggH[0.1][500] = 1.0;
        cpWeight_MapggH[0.1][600] = 1.0;
        cpWeight_MapggH[0.1][700] = 1.0;
        cpWeight_MapggH[0.1][800] = 1.0;
        cpWeight_MapggH[0.1][900] = 1.0;
        cpWeight_MapggH[0.1][1000] = 1.0;
        cpWeight_MapggH[0.1][1500] = 1.0;
        cpWeight_MapggH[0.1][2000] = 1.0;
        cpWeight_MapggH[0.1][2500] = 1.0;
        cpWeight_MapggH[0.1][3000] = 1.0;

	//Filling ContinuumWeights VBF
        cpWeight_MapVBF[1.0][200] = 1; 
        cpWeight_MapVBF[1.0][300] = 1; 
        cpWeight_MapVBF[1.0][400] = 1;
        cpWeight_MapVBF[1.0][500] = 1;
        cpWeight_MapVBF[1.0][600] = 1;
        cpWeight_MapVBF[1.0][700] = 1;
        cpWeight_MapVBF[1.0][800] = 1;
        cpWeight_MapVBF[1.0][900] = 1;
        cpWeight_MapVBF[1.0][1000] = 1;
        cpWeight_MapVBF[1.0][1500] = 1;
        cpWeight_MapVBF[1.0][2000] = 1;
        cpWeight_MapVBF[1.0][2500] = 1;
        cpWeight_MapVBF[1.0][3000] = 1;

        cpWeight_MapVBF[0.6][200] = 1;
        cpWeight_MapVBF[0.6][300] = 1;
        cpWeight_MapVBF[0.6][400] = 1;
        cpWeight_MapVBF[0.6][500] = 1;
        cpWeight_MapVBF[0.6][600] = 1;
        cpWeight_MapVBF[0.6][700] = 1;
        cpWeight_MapVBF[0.6][800] = 1;
        cpWeight_MapVBF[0.6][900] = 1;
        cpWeight_MapVBF[0.6][1000] = 1;
        cpWeight_MapVBF[0.6][1500] = 1;
        cpWeight_MapVBF[0.6][2000] = 1;
        cpWeight_MapVBF[0.6][2500] = 1;
        cpWeight_MapVBF[0.6][3000] = 1;

        cpWeight_MapVBF[0.3][200] = 1;
        cpWeight_MapVBF[0.3][300] = 1;
        cpWeight_MapVBF[0.3][400] = 1;
        cpWeight_MapVBF[0.3][500] = 1;
        cpWeight_MapVBF[0.3][600] = 1;
        cpWeight_MapVBF[0.3][700] = 1;
        cpWeight_MapVBF[0.3][800] = 1;
        cpWeight_MapVBF[0.3][900] = 1;
        cpWeight_MapVBF[0.3][1000] = 1;
        cpWeight_MapVBF[0.3][1500] = 1;
        cpWeight_MapVBF[0.3][2000] = 1;
        cpWeight_MapVBF[0.3][2500] = 1;
        cpWeight_MapVBF[0.3][3000] = 1;

        cpWeight_MapVBF[0.1][200] = 1;
        cpWeight_MapVBF[0.1][300] = 1;
        cpWeight_MapVBF[0.1][400] = 1;
        cpWeight_MapVBF[0.1][500] = 1;
        cpWeight_MapVBF[0.1][600] = 1;
        cpWeight_MapVBF[0.1][700] = 1;
        cpWeight_MapVBF[0.1][800] = 1;
        cpWeight_MapVBF[0.1][900] = 1;
        cpWeight_MapVBF[0.1][1000] = 1;
        cpWeight_MapVBF[0.1][1500] = 1;
        cpWeight_MapVBF[0.1][2000] = 1;
        cpWeight_MapVBF[0.1][2500] = 1;
        cpWeight_MapVBF[0.1][3000] = 1;

	if(isVBF){ continuumWeight = cpWeight_MapVBF[CP][heavyMass]; } else { continuumWeight = cpWeight_MapggH[CP][heavyMass]; }
	return continuumWeight;

 }

  }
}

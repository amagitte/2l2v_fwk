#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "TH1F.h"
#include "TSystem.h"

namespace utils
{
  namespace cmssw
  {
    //
    FactorizedJetCorrector *getJetCorrector(TString baseDir, bool isMC)
    {
      gSystem->ExpandPathName(baseDir);
      TString pf(isMC ? "MC" : "Data");
      
      //order matters: L1 -> L2 -> L3 (-> Residuals)
      std::vector<std::string> jetCorFiles;
      std::cout << baseDir+"/"+pf+"_L1FastJet_AK5PFchs.txt" << std::endl;
      jetCorFiles.push_back((baseDir+"/"+pf+"_L1FastJet_AK5PFchs.txt").Data());
      jetCorFiles.push_back((baseDir+"/"+pf+"_L2Relative_AK5PFchs.txt").Data());
      jetCorFiles.push_back((baseDir+"/"+pf+"_L3Absolute_AK5PFchs.txt").Data());
      if(!isMC) jetCorFiles.push_back((baseDir+"/"+pf+"_L2L3Residual_AK5PFchs.txt").Data());
      
      //init the parameters for correction
      std::vector<JetCorrectorParameters> corSteps;
      for(size_t i=0; i<jetCorFiles.size(); i++) corSteps.push_back(JetCorrectorParameters(jetCorFiles[i]));
      
      //return the corrector
      return new FactorizedJetCorrector(corSteps);
    }
 
    //
    const reco::Candidate *getGeneratorFinalStateFor(const reco::Candidate *p, bool isSherpa)
    {
      if(p==0) return 0;
      
      const reco::Candidate *prevState=p;
      do{	
	const reco::Candidate *nextState=0;
	int nDaughters = prevState->numberOfDaughters();
	for(int iDaughter=0; iDaughter<nDaughters; iDaughter++)
	  {
	    const reco::Candidate *dau = prevState->daughter(iDaughter);
	    if(dau==0) continue;
	    if(dau->pdgId()!= p->pdgId()) continue;
	    nextState=dau;	   
	    break;
	  }
	if(nextState==0) break;
	if(nextState==prevState) break;
	prevState=nextState;
      }while(1);
      return prevState;
    }

    //
    bool isBhadron(int pdgId)
    {
      int absid=abs(pdgId);
      return ( (absid>= 5122 && absid<= 5554) ||    //baryons
	       (absid>=20513 && absid<=20543) ||    //mesons
	       (absid>=10511 && absid<=10543) || 
	       (absid>=  511 && absid<=  545) ||
	       (absid>=  551 && absid<=  557) ||    //bbar mesons
	       (absid>=10551 && absid<=10557) ||
	       (absid>=100551 && absid<=100557) ||
	       (absid>=110551 && absid<=110557) ||
	       (absid>=200551 && absid<=200557) ||
	       (absid>=210551 && absid<=210557) );
    }

    //cf. https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    std::vector<float> smearJER(float pt, float eta, float genPt)
    {
      std::vector<float> toReturn(3,pt);
      if(genPt<=0) return toReturn;
      
      //
      eta=fabs(eta);
      double ptSF(1.0), ptSF_err(0.06);
      if(eta<0.5)                  { ptSF=1.052; ptSF_err=sqrt(pow(0.012,2)+pow(0.5*(0.062+0.061),2)); }
      else if(eta>=0.5 && eta<1.1) { ptSF=1.057; ptSF_err=sqrt(pow(0.012,2)+pow(0.5*(0.056+0.055),2)); }
      else if(eta>=1.1 && eta<1.7) { ptSF=1.096; ptSF_err=sqrt(pow(0.017,2)+pow(0.5*(0.063+0.062),2)); }
      else if(eta>=1.7 && eta<2.3) { ptSF=1.134; ptSF_err=sqrt(pow(0.035,2)+pow(0.5*(0.087+0.085),2)); }
      else if(eta>=2.3 && eta<5.0) { ptSF=1.288; ptSF_err=sqrt(pow(0.127,2)+pow(0.5*(0.155+0.153),2)); }
      
      toReturn[0]=TMath::Max(0.,(genPt+ptSF*(pt-genPt)));
      toReturn[1]=TMath::Max(0.,(genPt+(ptSF+ptSF_err)*(pt-genPt)));
      toReturn[2]=TMath::Max(0.,(genPt+(ptSF-ptSF_err)*(pt-genPt)));
      return toReturn;
    }

    //
    std::vector<float> smearJES(float pt, float eta, JetCorrectionUncertainty *jecUnc)
    {
      jecUnc->setJetEta(eta);
      jecUnc->setJetPt(pt);
      float relShift=fabs(jecUnc->getUncertainty(true));
      std::vector<float> toRet;
      toRet.push_back((1.0+relShift)*pt);
      toRet.push_back((1.0-relShift)*pt);
      return toRet;
    }


    //
    PuShifter_t getPUshifters(std::vector< float > &Lumi_distr, float puUnc)
    {
      Int_t NBins = Lumi_distr.size();
      TH1F *pu=new TH1F("putmp","",NBins,-0.5,float(NBins)-0.5);
      TH1F *puup=(TH1F *)pu->Clone("puuptmp");
      TH1F *pudown=(TH1F *)pu->Clone("pudowntmp");
      for(size_t i=0; i<Lumi_distr.size(); i++)  pu->SetBinContent(i+1,Lumi_distr[i]);
      
      for(int ibin=1; ibin<=pu->GetXaxis()->GetNbins(); ibin++)
	{
	  Double_t xval=pu->GetBinCenter(ibin);
	  TGraph *gr = new TGraph;
	  for(int ishift=-3; ishift<3; ishift++)
	    {
	      if(ibin+ishift<0) continue;
	      if(ibin+ishift>pu->GetXaxis()->GetNbins()) continue;
	      
	      gr->SetPoint(gr->GetN(),xval+ishift,pu->GetBinContent(ibin+ishift));
	    }
	  if(gr->GetN()>1)
	    {
	      Double_t newval(gr->Eval(xval*(1+puUnc)));
	      pudown->SetBinContent(ibin,newval>0?newval:0.0);
	      newval=gr->Eval(xval*(1-puUnc));
	      puup->SetBinContent(ibin,newval>0?newval:0.0);
	    }
	  delete gr;
	}
      puup->Scale(pu->Integral()/puup->Integral());
      pudown->Scale(pu->Integral()/pudown->Integral());
      std::cout << "getPUshifts will shift average PU by " << puup->GetMean()-pu->GetMean() << " / " << pudown->GetMean()-pu->GetMean() << std::endl; 
      
      puup->Divide(pu);    TGraph *puupWgt = new TGraph(puup);
      pudown->Divide(pu);  TGraph *pudownWgt = new TGraph(pudown);
      delete puup;
      delete pudown;  
      delete pu;
      
      PuShifter_t res(2);
      res[PUDOWN] = pudownWgt;
      res[PUUP]   = puupWgt;
      return res;
    }
    

    //
    Float_t getEffectiveArea(int id,float eta,int cone,TString isoSum)
    {
      Float_t Aeff(1.0);
      if(abs(id)==11)
	{
	  if(fabs(eta)<1.0)                         Aeff=(cone==3? 0.100 : 0.180);
	  else if(fabs(eta)>1.0 && fabs(eta)<1.479) Aeff=(cone==3? 0.120 : 0.200);
	  else if(fabs(eta)>1.479 && fabs(eta)<2.0) Aeff=(cone==3? 0.085 : 0.150);
	  else if(fabs(eta)>2.0 && fabs(eta)<2.2)   Aeff=(cone==3? 0.110 : 0.190);
	  else if(fabs(eta)>2.2 && fabs(eta)<2.3)   Aeff=(cone==3? 0.120 : 0.210);
	  else if(fabs(eta)>2.3 && fabs(eta)<2.4)   Aeff=(cone==3? 0.120 : 0.220);
	  else Aeff=0.14;
	}
      else if(abs(id)==22){
	if(isoSum=="chIso"){
	  if(fabs(eta)<1.0)                         Aeff=0.012;
          else if(fabs(eta)>1.0 && fabs(eta)<1.479) Aeff=0.010;
          else if(fabs(eta)>1.479 && fabs(eta)<2.0) Aeff=0.014;
          else if(fabs(eta)>2.0 && fabs(eta)<2.2)   Aeff=0.012;
          else if(fabs(eta)>2.2 && fabs(eta)<2.3)   Aeff=0.016;
          else if(fabs(eta)>2.3 && fabs(eta)<2.4)   Aeff=0.020;
          else                                      Aeff=0.012;
	}
	if(isoSum=="nhIso"){
	  if(fabs(eta)<1.0)                         Aeff=0.030;
          else if(fabs(eta)>1.0 && fabs(eta)<1.479) Aeff=0.057;
          else if(fabs(eta)>1.479 && fabs(eta)<2.0) Aeff=0.039;
          else if(fabs(eta)>2.0 && fabs(eta)<2.2)   Aeff=0.015;
          else if(fabs(eta)>2.2 && fabs(eta)<2.3)   Aeff=0.024;
          else if(fabs(eta)>2.3 && fabs(eta)<2.4)   Aeff=0.039;
          else                                      Aeff=0.072;
	}
	if(isoSum=="gIso"){
	  if(fabs(eta)<1.0)                         Aeff=0.148;
          else if(fabs(eta)>1.0 && fabs(eta)<1.479) Aeff=0.130;
          else if(fabs(eta)>1.479 && fabs(eta)<2.0) Aeff=0.112;
          else if(fabs(eta)>2.0 && fabs(eta)<2.2)   Aeff=0.216;
          else if(fabs(eta)>2.2 && fabs(eta)<2.3)   Aeff=0.262;
          else if(fabs(eta)>2.3 && fabs(eta)<2.4)   Aeff=0.260;
          else                                      Aeff=0.266;
	}
      }
      return Aeff;
    }


   double relIso(llvvLepton lep, double rho){
      if(abs(lep.id)==11){
          return (TMath::Max(lep.nhIso03+lep.gIso03-rho*utils::cmssw::getEffectiveArea(11,lep.electronInfoRef->sceta),double(0.))+lep.chIso03)/lep.pt();
      }else if(abs(lep.id)==13){
          return (TMath::Max(lep.nhIso04+lep.gIso04-0.5*lep.puchIso04,double(0.))+lep.chIso04)/lep.pt();
      }else{
          return -1;
      }
   }
    
    
    void getSingleMuTrigEff(const double& pt, const double& abseta, double& muontriggerefficiency){
      // Muon trigger/ID/Iso scale factors for efficiency are taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs                                                                                                                                           
      if(abseta>=0. && abseta <0.9){ // ABCD
	if(pt>=140. /*&& pt<500.*/) muontriggerefficiency=0.98041749810533507;
	if(pt>=25.  && pt<30.)      muontriggerefficiency=0.98372524384334614;
	if(pt>=30.  && pt<35.)      muontriggerefficiency=0.98406344315477012;
	if(pt>=35.  && pt<40.)      muontriggerefficiency=0.98391658181685537;
	if(pt>=40.  && pt<50.)      muontriggerefficiency=0.98345252700570363;
	if(pt>=50.  && pt<60.)      muontriggerefficiency=0.98429177039157478;
	if(pt>=60.  && pt<90.)      muontriggerefficiency=0.98467201842489449;
	if(pt>=90.  && pt<140.)     muontriggerefficiency=0.98091711658069591;
      }
      if(abseta>=0.9 && abseta <1.2){ // ABCD
	if(pt>=140. /*&& pt<500.*/) muontriggerefficiency=0.97127896196175556;
	if(pt>=25.  && pt<30.)      muontriggerefficiency=0.96838127559931908;
	if(pt>=30.  && pt<35.)      muontriggerefficiency=0.96538054889610103;
	if(pt>=35.  && pt<40.)      muontriggerefficiency=0.96696514151670487;
	if(pt>=40.  && pt<50.)      muontriggerefficiency=0.96667958160832501;
	if(pt>=50.  && pt<60.)      muontriggerefficiency=0.96273957552501865;
	if(pt>=60.  && pt<90.)      muontriggerefficiency=0.95952416834753307;
	if(pt>=90.  && pt<140.)     muontriggerefficiency=0.96444182461126438;
      }
      if(abseta>=1.2 && abseta <2.1){ // ABCD
	if(pt>=140. /*&& pt<500.*/) muontriggerefficiency=0.99416866829048334;
	if(pt>=25.  && pt<30.)      muontriggerefficiency=1.0051991254438037;
	if(pt>=30.  && pt<35.)      muontriggerefficiency=1.0013781590159485;
	if(pt>=35.  && pt<40.)      muontriggerefficiency=0.99616640424792002;
	if(pt>=40.  && pt<50.)      muontriggerefficiency=0.99425410141043047;
	if(pt>=50.  && pt<60.)      muontriggerefficiency=0.99054467301217797;
	if(pt>=60.  && pt<90.)      muontriggerefficiency=0.98829374192885855;
	if(pt>=90.  && pt<140.)     muontriggerefficiency=0.98187598993908232;
      }
    }
  
  }
  //
  std::string toLatexRounded(double value, double error, double systError,bool doPowers)
  {
    using namespace std;

    if(value==0.0 && error==0.0)return string("");
    
    if(!doPowers){
      char tmpchar[255];
      if(systError<0)
	sprintf(tmpchar,"$%.0f\\pm%.0f$",value,error);
      else
	sprintf(tmpchar,"$%.0f\\pm%.0f\\pm%.0f$",value,error,systError);
      return string(tmpchar);
    }
    
    double power = floor(log10(value));
    if(power<=-3)     {power=power+3;}
    else if(power>=2) {power=power-2;}
    else              {power=0;}
    
    value = value / pow(10,power);
    error = error / pow(10,power);
    if(systError>=0)systError = systError / pow(10,power);
    int ValueFloating;
    if(systError<0){
      ValueFloating = 1 + std::max(-1*log10(error),0.0);
    }else{
      ValueFloating = 1 + std::max(-1*log10(systError), std::max(-1*log10(error),0.0));
    }
    int ErrorFloating = ValueFloating;
    
    char tmpchar[255];
    if(power!=0){
      if(systError<0){
        sprintf(tmpchar,"$(%.*f\\pm%.*f)\\times 10^{%g}$",ValueFloating,value,ErrorFloating,error,power);
      }else{
        sprintf(tmpchar,"$(%.*f\\pm%.*f\\pm%.*f)\\times 10^{%g}$",ValueFloating,value,ErrorFloating,error,ErrorFloating,systError,power);
      }
      
    }else{
      if(systError<0){
        sprintf(tmpchar,"$%.*f\\pm%.*f$",ValueFloating,value,ErrorFloating,error);
      }else{
        sprintf(tmpchar,"$%.*f\\pm%.*f\\pm%.*f$",ValueFloating,value,ErrorFloating,error,ErrorFloating,systError);
      }
    }
    return string(tmpchar);
  }

  //
  void TLatexToTex(TString &expr)
  {
    expr = "$"+expr;
    expr += "$";
    expr.ReplaceAll("mu","\\mu"); 
    expr.ReplaceAll("_"," "); 
    expr.ReplaceAll("#","\\");
  }







	// loop on all the lumi blocks for an EDM file in order to count the number of events that are in a sample
	// this is useful to determine how to normalize the events (compute weight)
	unsigned long getMergeableCounterValue(const std::vector<std::string>& urls, std::string counter)
	{
	   unsigned long Total = 0;
	   for(unsigned int f=0;f<urls.size();f++){
	      TFile *file = TFile::Open(urls[f].c_str());      
	      fwlite::LuminosityBlock ls( file );
	      for(ls.toBegin(); !ls.atEnd(); ++ls){
		 fwlite::Handle<edm::MergeableCounter> nEventsTotalCounter;
		 nEventsTotalCounter.getByLabel(ls,counter.c_str());
		 if(!nEventsTotalCounter.isValid()){printf("Invalid nEventsTotalCounterH\n");continue;}
		 Total+= nEventsTotalCounter->value;
	      }
	   }
	   return Total;
	}



  void getMCPileupDistributionFromMiniAOD(fwlite::ChainEvent& ev, unsigned int Npu, std::vector<float>& mcpileup)
  {
    mcpileup.clear();
    mcpileup.resize(Npu);
    for(Long64_t ientry=0;ientry<ev.size();ientry++){
      ev.to(ientry);
      
      fwlite::Handle< std::vector<PileupSummaryInfo> > puInfoH;
      puInfoH.getByLabel(ev, "addPileupInfo");
      if(!puInfoH.isValid()){printf("collection PileupSummaryInfos with name addPileupInfo does not exist\n"); exit(0);}
      unsigned int ngenITpu = 0;
      for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++){
         if(it->getBunchCrossing()==0)      { ngenITpu += it->getPU_NumInteractions(); }
      }
      if(ngenITpu>=Npu){printf("ngenITpu is larger than vector size... vector is being resized, but you should check that all is ok!"); mcpileup.resize(ngenITpu+1);}
      mcpileup[ngenITpu]++;
    }
  }


  void getMCPileupDistribution(fwlite::ChainEvent& ev, unsigned int Npu, std::vector<float>& mcpileup)
  {
    mcpileup.clear();
    mcpileup.resize(Npu);
    for(Long64_t ientry=0;ientry<ev.size();ientry++){
      ev.to(ientry);
      
      fwlite::Handle< llvvGenEvent > genEventHandle;
      genEventHandle.getByLabel(ev, "llvvObjectProducersUsed");
//      if(!genEventHandle.isValid()){printf("llvvGenEvent Object NotFound\n");continue;}
      if(!genEventHandle.isValid()){continue;} //remove the warning... CAREFULL!
      unsigned int ngenITpu = (int)genEventHandle->ngenITpu;
      if(ngenITpu>=Npu){printf("ngenITpu is larger than vector size... vector is being resized, but you should check that all is ok!"); mcpileup.resize(ngenITpu+1);}
      mcpileup[ngenITpu]++;
    }
  }
  
  void getPileupNormalization(std::vector<float>& mcpileup, double* PUNorm, edm::LumiReWeighting* LumiWeights, utils::cmssw::PuShifter_t PuShifters){
    PUNorm[0]=0; PUNorm[1]=0; PUNorm[2]=0;
    double NEvents=0;
    for(unsigned int i=0;i<mcpileup.size();i++){
      NEvents+=mcpileup[i];
      double puWeight = LumiWeights->weight((int)i);
      PUNorm[0]+=mcpileup[i]*puWeight;
      PUNorm[1]+=mcpileup[i]*puWeight*PuShifters[utils::cmssw::PUDOWN]->Eval(i);
      PUNorm[2]+=mcpileup[i]*puWeight*PuShifters[utils::cmssw::PUUP  ]->Eval(i);
    }
    PUNorm[0]/=NEvents;
    PUNorm[1]/=NEvents;
    PUNorm[2]/=NEvents;
  }



  bool passTriggerPatternsAndGetName(edm::TriggerResultsByName& tr, std::string& pathName, std::string pattern){
     if(edm::is_glob(pattern)){
        std::vector< std::vector<std::string>::const_iterator > matches = edm::regexMatch(tr.triggerNames(), pattern);
        for(size_t t=0;t<matches.size();t++){
           if(tr.accept( matches[t]->c_str() ) ){pathName = *matches[t]; return true;}
        }
     }else{
        if(tr.accept( pattern.c_str() ) ) { pathName = pattern; return true;}
     }
     return false;
  }


  bool passTriggerPatterns(edm::TriggerResultsByName& tr, std::string pattern){
     if(edm::is_glob(pattern)){
        std::vector< std::vector<std::string>::const_iterator > matches = edm::regexMatch(tr.triggerNames(), pattern);
        for(size_t t=0;t<matches.size();t++){
           if(tr.accept( matches[t]->c_str() ) )return true;
        }
     }else{
        if(tr.accept( pattern.c_str() ) ) return true;
     }
     return false;
  }

  bool passTriggerPatterns(edm::TriggerResultsByName& tr, std::string pattern1, std::string pattern2, std::string pattern3, std::string pattern4){
     if(pattern1!="" && passTriggerPatterns(tr, pattern1))return true;
     if(pattern2!="" && passTriggerPatterns(tr, pattern2))return true;
     if(pattern3!="" && passTriggerPatterns(tr, pattern3))return true;
     if(pattern4!="" && passTriggerPatterns(tr, pattern4))return true;
     return false;
  }

  bool passTriggerPatterns(edm::TriggerResultsByName& tr, std::vector<std::string>& patterns){
     for(size_t p=0;p<patterns.size();p++){
        if(passTriggerPatterns(tr, patterns[p]))return true;
     }
     return false;
  }



}

#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "ZJetSkimTree.h"
#include <string>
#include <vector>

const Bool_t isDebug = false;

const Float_t lepPtCut = 20;

const Float_t muEtaCut = 2.4;
const Float_t muPtCut = 20;

const Float_t muChi2NDFCut = 10;
const Float_t muInnerD0Cut = 0.2;
const Float_t muInnerDzCut = 0.5;
const Int_t muMuonHitsCut = 0;
const Int_t muStationsCut = 1;
const Int_t muTrkLayersCut = 5;
const Int_t muPixelHitsCut = 0;

//see https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Spring15_selection_25ns
const Float_t eleEtaCut = 2.4;
const Float_t elePtCut = 20;

const Float_t barrelEndcapEta = 1.479;
const Int_t nBarrelEndcap = 2; // barrel = 0, endcap = 1
const Float_t eleSigmaIEtaIEta_VetoCut[nBarrelEndcap] = {0.0114, 0.0352};
const Float_t eleDEtaIn_VetoCut[nBarrelEndcap] = {0.0152, 0.0113};
const Float_t eleDPhiIn_VetoCut[nBarrelEndcap] = {0.216, 0.237};
const Float_t eleHOverE_VetoCut[nBarrelEndcap] = {0.181, 0.116};
const Float_t eleRelIsoWithEA_VetoCut[nBarrelEndcap] = {0.126, 0.144};
const Float_t eleOOEmooP_VetoCut[nBarrelEndcap] = {0.207, 0.174};
const Float_t eleD0_VetoCut[nBarrelEndcap] = {0.0564, 0.222};
const Float_t eleDZ_VetoCut[nBarrelEndcap] = {0.472, 0.921};
const Float_t eleMissingInnerHits_VetoCut[nBarrelEndcap] = {2, 3};

void makeZJetSkimTree(const std::string inFileName1 = "", const std::string inFileName2 = "")
{
  if(!strcmp(inFileName1.c_str(), "")){
    std::cout << "No inputs specified. return" << std::endl;
    return;
  }

  std::string outFileName = inFileName1;
  const std::string inString = ".root";
  TDatime* date = new TDatime();
  const std::string outString = Form("_%d_ZJetSkim.root", date->GetDate());
  std::size_t strIndex = 0;

  std::string tempOutName = outFileName;
  Int_t totOutNum = 0;

  if(isDebug) std::cout << __LINE__ << std::endl;

  while(tempOutName.find("/") != std::string::npos){
    tempOutName.replace(0, tempOutName.find("/")+1, "");
    totOutNum++;
  }

  Int_t outNum = 0;

  if(isDebug) std::cout << __LINE__ << std::endl;

  while(outNum < totOutNum-1){
    outFileName.replace(0, outFileName.find("/")+1, "");
    outNum++;
  }
  
  if(isDebug) std::cout << __LINE__ << std::endl;

  outFileName.replace(outFileName.find("/"), 1, "_");

  strIndex = outFileName.find(inString);
  if(!(strIndex == std::string::npos)){
    outFileName.replace(strIndex, inString.length(), outString);
  }

  if(isDebug) std::cout << __LINE__ << std::endl;

  TFile* outFile_p = new TFile(outFileName.c_str(), "UPDATE");
  skimTree_p = new TTree("skimTree", "skimTree");
  BookTree();

  std::vector<std::string>* inFileNames_p = new std::vector<std::string>;
  inFileNames_p->push_back(inFileName1);
  if(strcmp(inFileName2.c_str(), "") != 0) inFileNames_p->push_back(inFileName2);

  const Int_t nFiles = (Int_t)inFileNames_p->size();


  for(Int_t fileIter = 0; fileIter < nFiles; fileIter++){
    std::cout << "On file: " << fileIter << "/" << nFiles << std::endl;

    TFile* inFile_p = new TFile(inFileNames_p->at(fileIter).c_str(), "READ");
    TTree* lepTree_p = (TTree*)inFile_p->Get("ggHiNtuplizer/EventTree");
    TTree* jetTree_p = (TTree*)inFile_p->Get("akPu3PFJetAnalyzer/t");
    TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");

    std::vector<float>* muPt_p = 0;
    std::vector<float>* muPhi_p = 0;
    std::vector<float>* muEta_p = 0;
    std::vector<int>* muChg_p = 0;
    std::vector<float>* muChi2NDF_p = 0;
    std::vector<float>* muInnerD0_p = 0;
    std::vector<float>* muInnerDz_p = 0;
    std::vector<int>* muMuonHits_p = 0;
    std::vector<int>* muStations_p = 0;
    std::vector<int>* muTrkLayers_p = 0;
    std::vector<int>* muPixelHits_p = 0;

    std::vector<float>* elePt_p = 0;
    std::vector<float>* elePhi_p = 0;
    std::vector<float>* eleEta_p = 0;
    std::vector<int>* eleChg_p = 0;
    std::vector<float>* eleSigmaIEtaIEta_p = 0;
    std::vector<float>* eleDEtaAtVtx_p = 0;
    std::vector<float>* eleDPhiAtVtx_p = 0;
    std::vector<float>* eleHOverE_p = 0;
    std::vector<float>* eleD0_p = 0;
    std::vector<float>* eleDz_p = 0;
    

    lepTree_p->SetBranchStatus("*", 0);
    lepTree_p->SetBranchStatus("muPt", 1);
    lepTree_p->SetBranchStatus("muPhi", 1);
    lepTree_p->SetBranchStatus("muEta", 1);
    lepTree_p->SetBranchStatus("muCharge", 1);
    lepTree_p->SetBranchStatus("muChi2NDF", 1);
    lepTree_p->SetBranchStatus("muInnerD0", 1);
    lepTree_p->SetBranchStatus("muInnerDz", 1);
    lepTree_p->SetBranchStatus("muMuonHits", 1);
    lepTree_p->SetBranchStatus("muStations", 1);
    lepTree_p->SetBranchStatus("muTrkLayers", 1);
    lepTree_p->SetBranchStatus("muPixelHits", 1);
    
    lepTree_p->SetBranchAddress("muPt", &muPt_p);
    lepTree_p->SetBranchAddress("muPhi", &muPhi_p);
    lepTree_p->SetBranchAddress("muEta", &muEta_p);
    lepTree_p->SetBranchAddress("muCharge", &muChg_p);
    lepTree_p->SetBranchAddress("muChi2NDF", &muChi2NDF_p);
    lepTree_p->SetBranchAddress("muInnerD0", &muInnerD0_p);
    lepTree_p->SetBranchAddress("muInnerDz", &muInnerDz_p);
    lepTree_p->SetBranchAddress("muMuonHits", &muMuonHits_p);
    lepTree_p->SetBranchAddress("muStations", &muStations_p);
    lepTree_p->SetBranchAddress("muTrkLayers", &muTrkLayers_p);
    lepTree_p->SetBranchAddress("muPixelHits", &muPixelHits_p);    

    lepTree_p->SetBranchStatus("elePt", 1);
    lepTree_p->SetBranchStatus("elePhi", 1);
    lepTree_p->SetBranchStatus("eleEta", 1);
    lepTree_p->SetBranchStatus("eleCharge", 1);
    lepTree_p->SetBranchStatus("eleSigmaIEtaIEta", 1);
    lepTree_p->SetBranchStatus("eledEtaAtVtx", 1);
    lepTree_p->SetBranchStatus("eledPhiAtVtx", 1);
    lepTree_p->SetBranchStatus("eleHoverE", 1);
    lepTree_p->SetBranchStatus("eleD0", 1);
    lepTree_p->SetBranchStatus("eleDz", 1);
    
    lepTree_p->SetBranchAddress("elePt", &elePt_p);
    lepTree_p->SetBranchAddress("elePhi", &elePhi_p);
    lepTree_p->SetBranchAddress("eleEta", &eleEta_p);
    lepTree_p->SetBranchAddress("eleCharge", &eleChg_p);
    lepTree_p->SetBranchAddress("eleSigmaIEtaIEta", &eleSigmaIEtaIEta_p);
    lepTree_p->SetBranchAddress("eledEtaAtVtx", &eleDEtaAtVtx_p);
    lepTree_p->SetBranchAddress("eledPhiAtVtx", &eleDPhiAtVtx_p);
    lepTree_p->SetBranchAddress("eleHoverE", &eleHOverE_p);
    lepTree_p->SetBranchAddress("eleD0", &eleD0_p);
    lepTree_p->SetBranchAddress("eleDz", &eleDz_p);
    
    jetTree_p->SetBranchStatus("*", 0);
    jetTree_p->SetBranchStatus("nref", 1);
    jetTree_p->SetBranchStatus("jtpt", 1);
    jetTree_p->SetBranchStatus("jtphi", 1);
    jetTree_p->SetBranchStatus("jteta", 1);
    
    jetTree_p->SetBranchAddress("nref", &nJt_);
    jetTree_p->SetBranchAddress("jtpt", jtPt_);
    jetTree_p->SetBranchAddress("jtphi", jtPhi_);
    jetTree_p->SetBranchAddress("jteta", jtEta_);
    
    hiTree_p->SetBranchStatus("*", 0);
    hiTree_p->SetBranchStatus("run", 1);
    hiTree_p->SetBranchStatus("evt", 1);
    hiTree_p->SetBranchStatus("lumi", 1);
    hiTree_p->SetBranchStatus("hiBin", 1);
    hiTree_p->SetBranchStatus("vz", 1);
    
    hiTree_p->SetBranchAddress("run", &run_);
    hiTree_p->SetBranchAddress("evt", &evt_);
    hiTree_p->SetBranchAddress("lumi", &lumi_);
    hiTree_p->SetBranchAddress("hiBin", &hiBin_);
    hiTree_p->SetBranchAddress("vz", &vz_);
    
    if(isDebug) std::cout << __LINE__ << std::endl;
    
    if(isDebug) std::cout << __LINE__ << std::endl;
    
    const Int_t nEntries = (Int_t)lepTree_p->GetEntries();
    Int_t entryDiv = ((Int_t)(nEntries/20));
    
    if(isDebug) std::cout << __LINE__ << std::endl;
    
    for(Int_t entry = 0; entry < nEntries; entry++){
      if(isDebug) std::cout << __LINE__ << std::endl;
      
      if(entry%entryDiv == 0 && nEntries >= 10000) std::cout << "Entry # " << entry << "/" << nEntries << std::endl;
      if(isDebug) std::cout << __LINE__ << std::endl;
      
      lepTree_p->GetEntry(entry);
      jetTree_p->GetEntry(entry);
      hiTree_p->GetEntry(entry);
      
      if(TMath::Abs(vz_) > 15) continue;
      
      if(isDebug) std::cout << __LINE__ << std::endl;
      
      Float_t tempMuPt_[nLep];
      Float_t tempMuPhi_[nLep];
      Float_t tempMuEta_[nLep];
      Int_t tempMuChg_[nLep];
      
      Float_t tempElePt_[nLep];
      Float_t tempElePhi_[nLep];
      Float_t tempEleEta_[nLep];
      Int_t tempEleChg_[nLep];
      
      lepID_ = -999;
      lepM_ = -999;
      
      zPt_ = -999;
      zPhi_ = -999;
      zEta_ = -999;
      zM_ = -999;

      for(Int_t lepIter = 0; lepIter < nLep; lepIter++){
	lepPt_[lepIter] = -999;
	lepPhi_[lepIter] = -999;
	lepEta_[lepIter] = -999;
	lepChg_[lepIter] = -999;
	
	tempMuPt_[lepIter] = -999;
	tempMuPhi_[lepIter] = -999;
	tempMuEta_[lepIter] = -999;
	tempMuChg_[lepIter] = -999;
	
	tempElePt_[lepIter] = -999;
	tempElePhi_[lepIter] = -999;
	tempEleEta_[lepIter] = -999;
	tempEleChg_[lepIter] = -999;
      }
      
      if(isDebug) std::cout << __LINE__ << std::endl;
      
      const Int_t nMu = (Int_t)muPt_p->size();
      
      const Int_t muTrkLayersCut = 5;
      const Int_t muPixelHitsCut = 0;


      for(Int_t muIter = 0; muIter < nMu; muIter++){
	if(TMath::Abs(muEta_p->at(muIter)) > muEtaCut) continue;
	if(muPt_p->at(muIter) < muPtCut) continue;
	
	if(muChi2NDF_p->at(muIter) >= muChi2NDFCut) continue;
	if(TMath::Abs(muInnerD0_p->at(muIter)) > muInnerD0Cut) continue;
	if(TMath::Abs(muInnerDz_p->at(muIter)) > muInnerDzCut) continue;
	if(muMuonHits_p->at(muIter) <= muMuonHitsCut) continue;
	if(muStations_p->at(muIter) <= muStationsCut) continue;
	if(muTrkLayers_p->at(muIter) <= muTrkLayersCut) continue;
	if(muPixelHits_p->at(muIter) <= muPixelHitsCut) continue;

	if(muPt_p->at(muIter) > lepPt_[0]){
	  tempMuPt_[1] = tempMuPt_[0];
	  tempMuPhi_[1] = tempMuPhi_[0];
	  tempMuEta_[1] = tempMuEta_[0];
	  tempMuChg_[1] = tempMuChg_[0];
	  
	  tempMuPt_[0] = muPt_p->at(muIter);
	  tempMuPhi_[0] = muPhi_p->at(muIter);
	  tempMuEta_[0] = muEta_p->at(muIter);
	  tempMuChg_[0] = muChg_p->at(muIter);
	}
	else if(muPt_p->at(muIter) > tempMuPt_[1]){
	  tempMuPt_[1] = muPt_p->at(muIter);
	  tempMuPhi_[1] = muPhi_p->at(muIter);
	  tempMuEta_[1] = muEta_p->at(muIter);
	  tempMuChg_[1] = muChg_p->at(muIter);
	}
      }
      
      const Int_t nEle = (Int_t)elePt_p->size();

      for(Int_t eleIter = 0; eleIter < nEle; eleIter++){
	if(TMath::Abs(eleEta_p->at(eleIter)) > eleEtaCut) continue;
	if(elePt_p->at(eleIter) < elePtCut) continue;	

	Int_t eleEtaCutPos = 0;
	if(TMath::Abs(eleEta_p->at(eleIter)) > barrelEndcapEta) eleEtaCutPos = 1;
	
	if(eleSigmaIEtaIEta_p->at(eleIter) > eleSigmaIEtaIEta_VetoCut[eleEtaCutPos]) continue;
	if(TMath::Abs(eleDEtaAtVtx_p->at(eleIter)) > eleDEtaIn_VetoCut[eleEtaCutPos]) continue;
	if(TMath::Abs(eleDPhiAtVtx_p->at(eleIter)) > eleDPhiIn_VetoCut[eleEtaCutPos]) continue;
	if(eleHOverE_p->at(eleIter) > eleHOverE_VetoCut[eleEtaCutPos]) continue;
	if(TMath::Abs(eleD0_p->at(eleIter)) > eleD0_VetoCut[eleEtaCutPos]) continue;
	if(TMath::Abs(eleDz_p->at(eleIter)) > eleDZ_VetoCut[eleEtaCutPos]) continue;

	if(elePt_p->at(eleIter) > tempElePt_[0]){
	  tempElePt_[1] = tempElePt_[0];
	  tempElePhi_[1] = tempElePhi_[0];
	  tempEleEta_[1] = tempEleEta_[0];
	  tempEleChg_[1] = tempEleChg_[0];
	  
	  tempElePt_[0] = elePt_p->at(eleIter);
	  tempElePhi_[0] = elePhi_p->at(eleIter);
	  tempEleEta_[0] = eleEta_p->at(eleIter);
	  tempEleChg_[0] = eleChg_p->at(eleIter);
	}
	else if(elePt_p->at(eleIter) > tempElePt_[1]){
	  tempElePt_[1] = elePt_p->at(eleIter);
	  tempElePhi_[1] = elePhi_p->at(eleIter);
	  tempEleEta_[1] = eleEta_p->at(eleIter);
	  tempEleChg_[1] = eleChg_p->at(eleIter);
	}
      }
      
      if(tempMuPt_[0] > muPtCut && tempMuPt_[1] > muPtCut && tempElePt_[0] > elePtCut && tempElePt_[1] > elePtCut){
	if(tempMuPt_[0] > tempElePt_[0] && tempMuPt_[1] > tempElePt_[1]){
	  for(Int_t lepIter = 0; lepIter < nLep; lepIter++){
	    lepPt_[lepIter] = tempMuPt_[lepIter];
	    lepPhi_[lepIter] = tempMuPhi_[lepIter];
	    lepEta_[lepIter] = tempMuEta_[lepIter];
	    lepChg_[lepIter] = tempMuChg_[lepIter];
	}
	  lepID_ = muID;
	  lepM_ = muM;
	}
	else if(tempElePt_[0] > tempMuPt_[0] && tempElePt_[1] > tempMuPt_[1]){
	  for(Int_t lepIter = 0; lepIter < nLep; lepIter++){
	    lepPt_[lepIter] = tempElePt_[lepIter];
	    lepPhi_[lepIter] = tempElePhi_[lepIter];
	    lepEta_[lepIter] = tempEleEta_[lepIter];
	    lepChg_[lepIter] = tempEleChg_[lepIter];
	  }
	  lepID_ = eleID;
	  lepM_ = eleM;
	}
	else std::cout << "We got a tie, how do you want to handle?" << std::endl;
      }
      else if(tempMuPt_[0] > muPtCut && tempMuPt_[1] > muPtCut){
	for(Int_t lepIter = 0; lepIter < nLep; lepIter++){
	  lepPt_[lepIter] = tempMuPt_[lepIter];
	  lepPhi_[lepIter] = tempMuPhi_[lepIter];
	  lepEta_[lepIter] = tempMuEta_[lepIter];
	  lepChg_[lepIter] = tempMuChg_[lepIter];
	}
	lepID_ = muID;
	lepM_ = muM;
      }
      else if(tempElePt_[0] > elePtCut && tempElePt_[1] > elePtCut){
	for(Int_t lepIter = 0; lepIter < nLep; lepIter++){
	  lepPt_[lepIter] = tempElePt_[lepIter];
	  lepPhi_[lepIter] = tempElePhi_[lepIter];
	  lepEta_[lepIter] = tempEleEta_[lepIter];
	  lepChg_[lepIter] = tempEleChg_[lepIter];
	}
	lepID_ = eleID;
	lepM_ = eleM;
      }
      
      if(lepPt_[0] > lepPtCut && lepPt_[1] > lepPtCut){
	TLorentzVector lep1, lep2;
	lep1.SetPtEtaPhiM(lepPt_[0], lepEta_[0], lepPhi_[0], lepM_);
	lep2.SetPtEtaPhiM(lepPt_[1], lepEta_[1], lepPhi_[1], lepM_);

	TLorentzVector z = lep1+lep2;
	zPt_ = z.Pt();
	zPhi_ = z.Phi();
	zEta_ = z.Eta();
	zM_ = z.M();

	skimTree_p->Fill();
      }
    }

    inFile_p->Close();
    delete inFile_p;    
  }

  outFile_p->cd();
  TNamed pathStr1("pathStr1", inFileName1.c_str());
  pathStr1.Write("", TObject::kOverwrite);
  TNamed pathStr2("pathStr2", inFileName2.c_str());
  pathStr2.Write("", TObject::kOverwrite);
  skimTree_p->Write("", TObject::kOverwrite);

  delete date;
  
  outFile_p->Close();
  delete outFile_p;

  return;
}

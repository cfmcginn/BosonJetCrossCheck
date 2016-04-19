#ifndef ZJetSkimTree_h
#define ZJetSkimTree_h

#include "TTree.h"

#include <iostream>

TTree* skimTree_p = 0;

const Int_t nLep = 2;
const Int_t eleID = 11;
const Int_t muID = 13;
const Float_t eleM = .000510998910;
const Float_t muM = .1056583715;

UInt_t run_, lumi_;
ULong64_t evt_;
Int_t hiBin_;
Float_t vz_;

Int_t lepID_;
Float_t lepPt_[nLep];
Float_t lepPhi_[nLep];
Float_t lepEta_[nLep];
Int_t lepChg_[nLep];
Float_t lepM_;

Float_t zPt_;
Float_t zPhi_;
Float_t zEta_;
Float_t zM_;

const Int_t nMaxJets = 500;
Int_t nJt_;
Float_t jtPt_[nMaxJets];
Float_t jtPhi_[nMaxJets];
Float_t jtEta_[nMaxJets];

void BookTree()
{
  if(skimTree_p == NULL){
    std::cout << "BOOKTREE error; skimTree_p is NULL. return" << std::endl;
    return;
  }

  skimTree_p->Branch("run", &run_, "run/i");
  skimTree_p->Branch("evt", &evt_, "evt/l");
  skimTree_p->Branch("lumi", &lumi_, "lumi/i");
  skimTree_p->Branch("hiBin", &hiBin_, "hiBin/I");
  skimTree_p->Branch("vz", &vz_, "vz/F");

  skimTree_p->Branch("lepID", &lepID_, "lepID/I");
  skimTree_p->Branch("lepPt", lepPt_, Form("lepPt[%d]/F", nLep));
  skimTree_p->Branch("lepPhi", lepPhi_, Form("lepPhi[%d]/F", nLep));
  skimTree_p->Branch("lepEta", lepEta_, Form("lepEta[%d]/F", nLep));
  skimTree_p->Branch("lepChg", lepChg_, Form("lepChg[%d]/I", nLep));
  skimTree_p->Branch("lepM", &lepM_, "lepM/F");

  skimTree_p->Branch("zPt", &zPt_, "zPt/F");
  skimTree_p->Branch("zPhi", &zPhi_, "zPhi/F");
  skimTree_p->Branch("zEta", &zEta_, "zEta/F");
  skimTree_p->Branch("zM", &zM_, "zM/F");

  skimTree_p->Branch("nJt", &nJt_, "nJt/I");
  skimTree_p->Branch("jtPt", jtPt_, "jtPt[nJt]/F");
  skimTree_p->Branch("jtPhi", jtPhi_, "jtPhi[nJt]/F");
  skimTree_p->Branch("jtEta", jtEta_, "jtEta[nJt]/F");

  return;
}


void ReadTree()
{
  if(skimTree_p == NULL){
    std::cout << "READTREE error; skimTree_p is NULL. return" << std::endl;
    return;
  }


  skimTree_p->SetBranchStatus("*", 0);
  skimTree_p->SetBranchStatus("run", 1);
  skimTree_p->SetBranchStatus("evt", 1);
  skimTree_p->SetBranchStatus("lumi", 1);
  skimTree_p->SetBranchStatus("hiBin", 1);
  skimTree_p->SetBranchStatus("vz", 1);
  skimTree_p->SetBranchStatus("lepID", 1);
  skimTree_p->SetBranchStatus("lepPt", 1);
  skimTree_p->SetBranchStatus("lepPhi", 1);
  skimTree_p->SetBranchStatus("lepEta", 1);
  skimTree_p->SetBranchStatus("lepChg", 1);
  skimTree_p->SetBranchStatus("lepM", 1);
  skimTree_p->SetBranchStatus("zPt", 1);
  skimTree_p->SetBranchStatus("zPhi", 1);
  skimTree_p->SetBranchStatus("zEta", 1);
  skimTree_p->SetBranchStatus("zM", 1);
  skimTree_p->SetBranchStatus("nJt", 1);
  skimTree_p->SetBranchStatus("jtPt", 1);
  skimTree_p->SetBranchStatus("jtPhi", 1);
  skimTree_p->SetBranchStatus("jtEta", 1);

  skimTree_p->SetBranchAddress("run", &run_);
  skimTree_p->SetBranchAddress("evt", &evt_);
  skimTree_p->SetBranchAddress("lumi", &lumi_);
  skimTree_p->SetBranchAddress("hiBin", &hiBin_);
  skimTree_p->SetBranchAddress("vz", &vz_);
  skimTree_p->SetBranchAddress("lepID", &lepID_);
  skimTree_p->SetBranchAddress("lepPt", lepPt_);
  skimTree_p->SetBranchAddress("lepPhi", lepPhi_);
  skimTree_p->SetBranchAddress("lepEta", lepEta_);
  skimTree_p->SetBranchAddress("lepChg", lepChg_);
  skimTree_p->SetBranchAddress("lepM", &lepM_);
  skimTree_p->SetBranchAddress("zPt", &zPt_);
  skimTree_p->SetBranchAddress("zPhi", &zPhi_);
  skimTree_p->SetBranchAddress("zEta", &zEta_);
  skimTree_p->SetBranchAddress("zM", &zM_);
  skimTree_p->SetBranchAddress("nJt", &nJt_);
  skimTree_p->SetBranchAddress("jtPt", jtPt_);
  skimTree_p->SetBranchAddress("jtPhi", jtPhi_);
  skimTree_p->SetBranchAddress("jtEta", jtEta_);

  return;
}

#endif

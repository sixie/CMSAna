//================================================================================================
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2F.h>                   // 2D histograms
#include <TGraphAsymmErrors.h>     
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings



// define structures to read in ntuple
#include "CMSAna/DataTree/interface/DataTreeDefs.hh"
#include "CMSAna/DataTree/interface/TEventInfo.hh"
#include "CMSAna/DataTree/interface/TElectron.hh"
#include "CMSAna/DataTree/interface/TPhoton.hh"
#include "CMSAna/DataTree/interface/TMuon.hh"
#include "CMSAna/DataTree/interface/TJet.hh"
#include "CMSAna/DataTree/interface/TPFCandidate.hh"
#include "CMSAna/DataTree/interface/TGenParticle.hh"
#include "CMSAna/DataTree/interface/TGenJet.hh"

#include "CMSAna/ObjectStudies/interface/ElectronTree.h"
#include "CMSAna/Utils/LeptonTools.hh"



#endif

//*************************************************************************************************
//Convert int to string
//*************************************************************************************************
string IntToString(int i) {
  char temp[100];
  sprintf(temp, "%d", i);
  string str = temp;
  return str;
}

//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, const char* tname)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = TFile::Open(infname,"read");
  assert(inf);

  TTree* t = (TTree*)inf->Get(tname);
  
  if (!t) {
    TDirectory *dir = (TDirectory*)inf->FindObjectAny("HwwNtuplerMod");
    if (!dir) {
      cout << "Cannot get Directory HwwNtuplerMod from file " << infname << endl;
      assert(dir);
    }
    t = (TTree*)dir->Get(tname);
  }

  if (!t) {
    cout << "Cannot get Tree with name " << tname << " from file " << infname << endl;
  }
  assert(t);


  if (verbose) {
    cout << "---\tRecovered tree " << t->GetName()
	 << " with "<< t->GetEntries() << " entries" << endl;
  }
  
  return t;
}


//=== FUNCTION DECLARATIONS ======================================================================================



// print event dump
void MakeNtuple(const string inputFilename,  const string outputFilename);



//=== MAIN MACRO =================================================================================================
void MakeHZZFakeEleNtuple(Int_t Sample = 0) {

  //default
  if (Sample == 100) {
    MakeNtuple("/data1/sixie/ntuples/BACON/HZZ4L/V1/BACONNtuple_100.root","ElectronNtuple.HZZ4L.100.root");
  }
} 


void MakeHZZFakeEleNtuple(const string inputFilename, const string outputFilename) {
  MakeNtuple(inputFilename,outputFilename);
}



void MakeNtuple(const string inputFilename, const string outputFilename)
{  
  gBenchmark->Start("WWTemplate");

  //*****************************************************************************************
  //Define Data Era
  //*****************************************************************************************
  UInt_t DataEra = kDataEra_NONE;
  DataEra = kDataEra_2012_MC;

  cout << "using DataEra = " << DataEra << endl;


  //*****************************************************************************************
  //Setup Output Tree
  //*****************************************************************************************
  TFile *outputFile = new TFile(outputFilename.c_str(), "RECREATE");
  cmsana::ElectronTree *eleTree = new cmsana::ElectronTree;
  eleTree->CreateTree(cmsana::ElectronTree::kEleTreeLight);
  eleTree->tree_->SetAutoFlush(0);

  UInt_t NElectronsFilled = 0;
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  cmsana::TEventInfo *info    = new cmsana::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("cmsana::TElectron");
  TClonesArray *muonArr = new TClonesArray("cmsana::TMuon");
  TClonesArray *genjetArr = new TClonesArray("cmsana::TGenJet");
  TClonesArray *genparticleArr = new TClonesArray("cmsana::TGenParticle");

  Int_t NEvents = 0;

  //********************************************************
  // Get Tree
  //********************************************************
  eventTree = getTreeFromFile(inputFilename.c_str(),"Events"); 
  TBranch *infoBr;
  TBranch *electronBr;
  TBranch *muonBr;
  TBranch *genjetBr;
  TBranch *genparticleBr;


  //*****************************************************************************************
  //Loop over Data Tree
  //*****************************************************************************************
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("Muon", &muonArr); muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("GenJet", &genjetArr); genjetBr = eventTree->GetBranch("GenJet");
  eventTree->SetBranchAddress("GenParticle", &genparticleArr); genparticleBr = eventTree->GetBranch("GenParticle");

  cout << "InputFile " << inputFilename << " --- Total Events : " << eventTree->GetEntries() << endl;
  for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
    infoBr->GetEntry(ientry);
		
    if (ientry % 10000 == 0) cout << "Event " << ientry << endl;

    Double_t eventweight = info->eventweight;

    //********************************************************
    // Load the branches
    //********************************************************
    electronArr->Clear(); 
    muonArr->Clear();
    genparticleArr->Clear(); 
    genjetArr->Clear();    
    electronBr->GetEntry(ientry);
    muonBr->GetEntry(ientry);
    genjetBr->GetEntry(ientry);
    genparticleBr->GetEntry(ientry);

    //********************************************************
    // Pileup Energy Density
    //********************************************************
    Double_t rho = info->RhoKt6PFJets;

    //********************************************************
    // Only use events where z->ee or z->mumu was reconstructed
    //********************************************************
    int NElectrons  = 0;
    int NMuons  = 0;
    for(Int_t i=0; i<electronArr->GetEntries(); i++) {
      const cmsana::TElectron *ele = (cmsana::TElectron*)((*electronArr)[i]);

      bool match = false;
      for(Int_t j=0; j<genparticleArr->GetEntries(); j++) {
        const cmsana::TGenParticle *gen = (cmsana::TGenParticle*)((*genparticleArr)[j]);
        if (abs(gen->pdgid) == 11 && gen->status == 1 && (abs(gen->motherPdgID) == 23 || gen->motherPdgID == gen->pdgid)) {
          if (cmsana::deltaR(ele->eta,ele->phi,gen->eta,gen->phi) < 0.1) match = true;
        }
      }
      if (match) NElectrons++;
    }
    for(Int_t i=0; i<muonArr->GetEntries(); i++) {
      const cmsana::TMuon *mu = (cmsana::TMuon*)((*muonArr)[i]);

      bool match = false;
      for(Int_t j=0; j<genparticleArr->GetEntries(); j++) {
        const cmsana::TGenParticle *gen = (cmsana::TGenParticle*)((*genparticleArr)[j]);
        if (abs(gen->pdgid) == 13 && gen->status == 1 && (abs(gen->motherPdgID) == 23 || gen->motherPdgID == gen->pdgid)) {
          if (cmsana::deltaR(mu->eta,mu->phi,gen->eta,gen->phi) < 0.1) match = true;
        }
      }
      if (match) NMuons++;
    }
    if (!(NElectrons >= 2 || NMuons >= 2)) continue;


    //********************************************************
    // Loop over Gen Jets
    //********************************************************
    vector<float> eleUsedEta;
    vector<float> eleUsedPhi;

    for(Int_t i=0; i<genjetArr->GetEntries(); i++) {

      const cmsana::TGenJet *genjet = (cmsana::TGenJet*)((*genjetArr)[i]);
      
      //***************************************
      //Don't use real electrons
      //***************************************
      bool matchRealEle = false;
      bool matchRealMu = false;
      bool matchRealTau = false;
      for(Int_t j=0; j<genparticleArr->GetEntries(); j++) {
        const cmsana::TGenParticle *gen = (cmsana::TGenParticle*)((*genparticleArr)[j]);
        if (abs(gen->pdgid) == 11 && gen->status == 1 && (abs(gen->motherPdgID) == 23 || gen->motherPdgID == gen->pdgid)) {
          if (cmsana::deltaR(genjet->eta,genjet->phi,gen->eta,gen->phi) < 0.3) matchRealEle = true;
        }
        if (abs(gen->pdgid) == 13 && gen->status == 1 && (abs(gen->motherPdgID) == 23 || gen->motherPdgID == gen->pdgid)) {
          if (cmsana::deltaR(genjet->eta,genjet->phi,gen->eta,gen->phi) < 0.3) matchRealMu = true;
        }
        if (abs(gen->pdgid) == 15 && (abs(gen->motherPdgID) == 23 || gen->motherPdgID == gen->pdgid)) {
          if (cmsana::deltaR(genjet->eta,genjet->phi,gen->eta,gen->phi) < 0.3) matchRealTau = true;
        }
      }
      if (matchRealEle || matchRealMu || matchRealTau) continue;

      //***************************************
      //Don't use duplicate electrons
      //***************************************
      bool isDuplicate = false;
      for(Int_t j=0; j<eleUsedEta.size(); j++) {
        if (cmsana::deltaR(genjet->eta,genjet->phi,eleUsedEta[j],eleUsedPhi[j]) < 0.1) {
          isDuplicate = true;
          break;
        }
      }
      if (isDuplicate) continue;
      
      //***************************************
      //Find corresponding reco electron
      //***************************************
      double elept = -999;
      double eleeta = -999;
      double elephi = -999;
      
      for(Int_t j=0; j<electronArr->GetEntries(); j++) {
        const cmsana::TElectron *ele = (cmsana::TElectron*)((*electronArr)[j]);
        if (cmsana::deltaR(genjet->eta,genjet->phi,ele->eta,ele->phi) < 0.1) {
          elept = 0;
          if (ele->pt > 0) elept = ele->pt; //protect against possible weird negative pt values
          eleeta = ele->eta;
          elephi = ele->phi;
        }
      }
      
      //***************************************
      //Fill gen electron
      //***************************************
      eleTree->fWeight = eventweight;
      eleTree->fRunNumber = info->runNum;
      eleTree->fLumiSectionNumber = info->lumiSec;
      eleTree->fEventNumber = info->evtNum;
      eleTree->fEleEventNumberParity = (info->evtNum % 2 == 0);
      eleTree->fRho = rho;
      eleTree->fEleGenPt = genjet->pt;
      eleTree->fEleGenEta = genjet->eta; 
      eleTree->fEleGenPhi = genjet->phi;
      eleTree->fElePt = elept;
      eleTree->fEleEta = eleeta;
      eleTree->fElePhi = elephi;
      eleTree->fPdgId = genjet->matchedPdgId;
      eleTree->tree_->Fill();
      NElectronsFilled++;     

    } //loop over gen particles

  }//loop over events

  delete info;
  delete electronArr;
  delete genparticleArr;

  cout << "Total Electrons: " << NElectronsFilled << endl;

  outputFile->Write();
  outputFile->Close();
  
  delete outputFile;
  if (eleTree) delete eleTree;

  gBenchmark->Show("WWTemplate");       
} 



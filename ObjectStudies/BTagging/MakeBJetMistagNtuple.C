//================================================================================================
// root -l CMSAna/ObjectStudies/BTagging/MakeBJetMistagNtuple.C+'("root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-02/qcd_pt30To50-START53_V7A/BACONNtuple_80_1_80U.root","BJetMistagNtuple.test.root",0)'
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
#include "CMSAna/DataTree/interface/TPhoton.hh"
#include "CMSAna/DataTree/interface/TJet.hh"
#include "CMSAna/DataTree/interface/TPFCandidate.hh"
#include "CMSAna/DataTree/interface/TGenParticle.hh"
#include "CMSAna/DataTree/interface/TGenJet.hh"
#include "CMSAna/ObjectStudies/interface/EfficiencyTree.h"
#include "CMSAna/Utils/PhotonTools.hh"



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
    TDirectory *dir = (TDirectory*)inf->FindObjectAny("Events");
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


void MakeBJetMistagNtuple(const string inputFilename, const string outputFilename, int denominatorType = 0)
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
  cmsana::EfficiencyTree *effTree = new cmsana::EfficiencyTree;
  effTree->CreateTree();
  effTree->tree_->SetAutoFlush(0);

  UInt_t NDenominatorsFilled = 0;
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  cmsana::TEventInfo *info    = new cmsana::TEventInfo();
  TClonesArray *jetArr = new TClonesArray("cmsana::TJet");
  TClonesArray *pfcandidateArr = new TClonesArray("cmsana::TPFCandidate");
  TClonesArray *genparticleArr = new TClonesArray("cmsana::TGenParticle");
  TClonesArray *genjetArr = new TClonesArray("cmsana::TGenJet");

  Int_t NEvents = 0;



  //********************************************************
  // Get Tree
  //********************************************************
  eventTree = getTreeFromFile(inputFilename.c_str(),"Events"); 
  TBranch *infoBr;
  TBranch *jetBr;
  TBranch *pfcandidateBr;
  TBranch *genparticleBr;
  TBranch *genjetBr;
  
  
  //*****************************************************************************************
  //Loop over Data Tree
  //*****************************************************************************************
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);  infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("PFJet", &jetArr); jetBr  = eventTree->GetBranch("PFJet");
  eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr); pfcandidateBr = eventTree->GetBranch("PFCandidate");
  eventTree->SetBranchAddress("GenParticle", &genparticleArr); genparticleBr = eventTree->GetBranch("GenParticle");
  eventTree->SetBranchAddress("GenJet", &genjetArr); genjetBr = eventTree->GetBranch("GenJet");

  cout << "InputFile " << inputFilename << " --- Total Events : " << eventTree->GetEntries() << endl;
  for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
    infoBr->GetEntry(ientry);
		
    if (ientry % 1 == 0) cout << "Event " << ientry << endl;

    Double_t eventweight = info->eventweight;

    //********************************************************
    // Load the branches
    //********************************************************
    jetArr->Clear(); 
    pfcandidateArr->Clear(); 
    genparticleArr->Clear(); 
    genjetArr->Clear(); 
    jetBr->GetEntry(ientry);
    pfcandidateBr->GetEntry(ientry);
    genparticleBr->GetEntry(ientry);
    genjetBr->GetEntry(ientry);


    //********************************************************
    // Pileup Energy Density
    //********************************************************
    Double_t rho = info->RhoKt6PFJets;
    UInt_t EAEra = kDataEra_2012_Data;


    //********************************************************
    // genjets denominator
    //********************************************************
    if ( denominatorType == 0) {
      for(Int_t k=0; k<genjetArr->GetEntries(); k++) {
        const cmsana::TGenJet *genjet = (cmsana::TGenJet*)((*genjetArr)[k]);
               
        //some kinematic cuts to save ntupling time
        if (!(genjet->pt > 20)) continue;
        if (!(fabs(genjet->eta) < 2.5)) continue;

        //make sure it doesn't match to a real b-jet
        bool matchBJet = false;
        for(Int_t i=0; i<genparticleArr->GetEntriesFast(); i++) {
          const cmsana::TGenParticle *p = (cmsana::TGenParticle*)((*genparticleArr)[i]);
          if ( p->pt > 20 && abs(p->pdgid) == 5 && p->status == 2 ) {
            if (cmsana::deltaR(p->eta,p->phi,genjet->eta,genjet->phi) < 0.5) {
              matchBJet = true;
              break;
            }
          }
        }
        if (matchBJet) continue;


        bool pass = false;
        //Find matching jet
        for(Int_t i=0; i<jetArr->GetEntriesFast(); i++) {
          const cmsana::TJet *jet = (cmsana::TJet*)((*jetArr)[i]);
          if (!(jet->pt > 20)) continue;
          if (!(fabs(jet->eta) < 2.5)) continue;
          
          //match in dR?
          double DR = cmsana::deltaR(jet->eta,jet->phi,genjet->eta,genjet->phi);
          if (!(DR < 0.5)) continue;

          if (!(jet->CombinedSecondaryVertexBJetTagsDisc > 0.679)) continue;
          
          pass = true;
          break;
        }
        
        effTree->weight = eventweight;
        effTree->mass = 0;
        effTree->pt = genjet->pt;
        effTree->eta = genjet->eta;
        effTree->phi = genjet->phi;
        effTree->rho = info->RhoKt6PFJets;
        effTree->q = 0;
        effTree->npv = info->nGoodPV;
        effTree->npu = info->nPU;
        effTree->run = info->runNum;
        effTree->lumi = info->lumiSec;
        effTree->event = info->evtNum;
        effTree->pass = pass;

        //***********************
        //Fill Denominator
        //***********************
        NDenominatorsFilled++;
        effTree->tree_->Fill();

      }//loop over genjet denominators
    } // if denominator is genjet


  } //loop over events
  
  cout << "Total Denominators: " << NDenominatorsFilled << endl;

  delete info;
  delete jetArr;
  delete pfcandidateArr;
  delete genparticleArr;
  delete genjetArr;

  outputFile->Write();
  outputFile->Close();
  
  delete outputFile;
  if (effTree) delete effTree;
} 



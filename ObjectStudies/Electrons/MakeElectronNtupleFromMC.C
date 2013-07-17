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
// #include "CMSAna/Utils/LeptonID.hh"
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
void MakeElectronNtupleFromMC(Int_t Sample = 0) {
 
  //Phase 1 Samples: Age0
  if (Sample == 1) {
    MakeNtuple("root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-03_SLHC4/DYToEEAge0START_STAR17_61_V1A/BACONNtuple_1_1_Pmm.root","ElectronNtuple.DYToEEAge0START.1.root");
  }


} 


void MakeElectronNtupleFromMC(const string inputFilename, const string outputFilename) {
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
  eleTree->CreateTree();
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
  TClonesArray *pfcandidateArr = new TClonesArray("cmsana::TPFCandidate");
  TClonesArray *genparticleArr = new TClonesArray("cmsana::TGenParticle");

  Int_t NEvents = 0;



    //********************************************************
    // Get Tree
    //********************************************************
    eventTree = getTreeFromFile(inputFilename.c_str(),"Events"); 
    TBranch *infoBr;
    TBranch *electronBr;
    TBranch *muonBr;
    TBranch *pfcandidateBr;
    TBranch *genparticleBr;


    //*****************************************************************************************
    //Loop over Data Tree
    //*****************************************************************************************
    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr); pfcandidateBr = eventTree->GetBranch("PFCandidate");
    eventTree->SetBranchAddress("GenParticle", &genparticleArr); genparticleBr = eventTree->GetBranch("GenParticle");

    cout << "InputFile " << inputFilename << " --- Total Events : " << eventTree->GetEntries() << endl;
//     for(UInt_t ientry=0; ientry < 10000; ientry++) {       	
    for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
      if (ientry % 100 == 0) cout << "Event " << ientry << endl;

      Double_t eventweight = info->eventweight;

      //********************************************************
      // Load the branches
      //********************************************************
      electronArr->Clear(); 
      pfcandidateArr->Clear(); 
      genparticleArr->Clear(); 
      electronBr->GetEntry(ientry);
      pfcandidateBr->GetEntry(ientry);
      genparticleBr->GetEntry(ientry);


      //********************************************************
      // Pileup Energy Density
      //********************************************************
      Double_t rhoEleIso = info->RhoKt6PFJets;
      UInt_t EleEAEra = kDataEra_2012_Data;

      //********************************************************
      // Loop Over Electrons
      //********************************************************
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const cmsana::TElectron *ele = (cmsana::TElectron*)((*electronArr)[i]);

        //********************************************************
        // Select MC Truth Electrons
        //********************************************************              
        double minDR = 9999;
        const cmsana::TGenParticle *matchedGenElectron =0;
        for(Int_t k=0; k<genparticleArr->GetEntries(); k++) {
          const cmsana::TGenParticle *gen = (cmsana::TGenParticle*)((*genparticleArr)[k]);
          
          //status 1 match
          if (abs(gen->pdgid) == 11 && (gen->status == 1)) {
            double tmpDR = cmsana::deltaR( ele->eta, ele->phi, gen->eta , gen->phi);
            if (tmpDR < minDR && tmpDR < 0.3) {
              minDR = tmpDR;
              matchedGenElectron = gen;
            }
          }
        }
        //Don't use electrons from tau decays
        if (!matchedGenElectron || abs(matchedGenElectron->motherPdgID) == 15) {
          if (ele->pt > 10) {
//             cout << "fake ele: " << ele->pt << " " << ele->eta << " " << ele->phi << "\n";
//             for(Int_t k=0; k<genparticleArr->GetEntries(); k++) {
//               const cmsana::TGenParticle *gen = (cmsana::TGenParticle*)((*genparticleArr)[k]);
//               cout << "gen: " << k << " : " << gen->pdgid << " " << gen->status << " " << gen->pt << " " << gen->eta << " " << gen->phi << "\n";
//             }
          }

          continue;
        } else {
//           cout << "real ele: " << ele->pt << " " << ele->eta << " " << ele->phi << "\n";
        }

        //Fill These Electrons
        NElectronsFilled++;
        FillElectronTree( eleTree, ele, i, pfcandidateArr, rhoEleIso, EleEAEra, 
                          info->nGoodPV, info->runNum, info->lumiSec, info->evtNum, 1.0);
        
      } //loop over electrons

    }

    cout << "Total Electrons: " << NElectronsFilled << endl;


    delete info;
    delete electronArr;
    delete pfcandidateArr;
    delete genparticleArr;

    cout << "Total Electrons: " << NElectronsFilled << endl;

    outputFile->Write();
    outputFile->Close();
  
    delete outputFile;
    if (eleTree) delete eleTree;

    gBenchmark->Show("WWTemplate");       
} 



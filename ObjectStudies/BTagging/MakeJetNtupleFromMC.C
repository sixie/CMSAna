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
#include "CMSAna/ObjectStudies/interface/JetTree.h"
#include "CMSAna/Utils/JetTools.hh"
#include "CMSAna/Utils/CommonDefs.hh"
#include "CMSAna/Utils/CommonTools.hh"

// For Jet energy corrections
#include "CMSAna/JetEnergyCorrections/interface/FactorizedJetCorrector.h"
#include "CMSAna/JetEnergyCorrections/interface/JetCorrectorParameters.h"
#include "CMSAna/Utils/JetEnergyCorrections.hh"


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
void MakeNtuple(const string inputFilename,  const string outputFilename, int firstEvent, int processNEvents);

void MakeNtuple(const string inputFilename,  const string outputFilename) {
  MakeNtuple(inputFilename, outputFilename, -1,-1);
}

//=== MAIN MACRO =================================================================================================
void MakeJetNtupleFromMC(Int_t Sample = 0) {
 
  //Phase 1 Samples: Age0
  if (Sample == 1) {
    MakeNtuple("root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-04_SLHC4/WZHgg-125-Age0_STAR17_61_V1A/BACONNtuple_10_1_azS.root", "jetNtuple.summer12.real.root",  -1,-1);
  } else if (Sample==2) {
    MakeNtuple("root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-04/ttHbb-125-START53_V7A/BACONNtuple_15_1_9MJ.root", "jetNtuple.summer12.fake.root",  -1, -1);
  } else if (Sample == 11) {
    MakeNtuple("root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-04_SLHC4/WZHgg-125-Age0_STAR17_61_V1A/BACONNtuple_10_1_azS.root", "jetNtuple.highpileup.wzhgg.root", 100, 0);
  } else if (Sample == 12) {
    MakeNtuple("root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-04_SLHC4/ttHbb-125-Age0_STAR17_61_V1A/BACONNtuple_17_1_a9l.root", "jetNtuple.highpileup.fake2.root", 100, 0);
  }


} 

void MakeJetNtupleFromMC(const string inputFilename, const string outputFilename, int splitIntoNJobs, int jobNumber) {
  MakeNtuple(inputFilename,outputFilename, splitIntoNJobs, jobNumber);
}

void MakeJetNtupleFromMC(const string inputFilename, const string outputFilename) {
  MakeNtuple(inputFilename,outputFilename);
}



void MakeNtuple(const string inputFilename, const string outputFilename, int splitIntoNJobs, int jobNumber)
{  
  gBenchmark->Start("WWTemplate");

  //*****************************************************************************************
  //Define Data Era
  //*****************************************************************************************
  UInt_t DataEra = kDataEra_NONE;
  DataEra = kDataEra_2012_MC;

  cout << "using DataEra = " << DataEra << endl;

  //*****************************************************************************************
  //Setup Jet Energy Corrections
  //*****************************************************************************************
  std::vector<cmsana::JetCorrectorParameters> correctionParameters;
  correctionParameters.push_back(cmsana::JetCorrectorParameters( ( getenv("CMSSW_BASE") + string("/src/CMSAna/JetEnergyCorrections/data/GR_R_52_V9_L1FastJet_AK5PF.txt")).c_str()));
  correctionParameters.push_back(cmsana::JetCorrectorParameters( ( getenv("CMSSW_BASE") + string("/src/CMSAna/JetEnergyCorrections/data/GR_R_52_V9_L2Relative_AK5PF.txt")).c_str()));
  correctionParameters.push_back(cmsana::JetCorrectorParameters( ( getenv("CMSSW_BASE") + string("/src/CMSAna/JetEnergyCorrections/data/GR_R_52_V9_L3Absolute_AK5PF.txt")).c_str()));
  cmsana::FactorizedJetCorrector *JetCorrector = new cmsana::FactorizedJetCorrector(correctionParameters);


  //*****************************************************************************************
  //Setup Output Tree
  //*****************************************************************************************
  TFile *outputFile = new TFile(outputFilename.c_str(), "RECREATE");
  cmsana::JetTree *jetTree = new cmsana::JetTree;
  jetTree->CreateTree();
  jetTree->tree_->SetAutoFlush(0);

  UInt_t NJetsFilled = 0;
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  cmsana::TEventInfo *info    = new cmsana::TEventInfo();
  TClonesArray *jetArr = new TClonesArray("cmsana::TJet");
  TClonesArray *genparticleArr = new TClonesArray("cmsana::TGenParticle");
  TClonesArray *genjetArr = new TClonesArray("cmsana::TGenJet");

  Int_t NEvents = 0;



  //********************************************************
  // Get Tree
  //********************************************************
  eventTree = getTreeFromFile(inputFilename.c_str(),"Events"); 
  TBranch *infoBr;
  TBranch *jetBr;
  TBranch *genparticleBr;
  TBranch *genjetBr;
  
  
  //*****************************************************************************************
  //Loop over Data Tree
  //*****************************************************************************************
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("PFJet", &jetArr); jetBr = eventTree->GetBranch("PFJet");
  eventTree->SetBranchAddress("GenParticle", &genparticleArr); genparticleBr = eventTree->GetBranch("GenParticle");
  eventTree->SetBranchAddress("GenJet", &genjetArr); genjetBr = eventTree->GetBranch("GenJet");

  cout << "InputFile " << inputFilename << " --- Total Events : " << eventTree->GetEntries() << endl;

  int firstEvent = 0;
  int lastEvent = eventTree->GetEntries();
  if ( splitIntoNJobs > 0) {
    int NEventsPerJob = eventTree->GetEntries() / splitIntoNJobs;
    firstEvent = NEventsPerJob*jobNumber;
    lastEvent = NEventsPerJob*(jobNumber+1);
  }

  for(UInt_t ientry=firstEvent; ientry < lastEvent; ientry++) {       	
    infoBr->GetEntry(ientry);
		
    if (ientry % 1 == 0) cout << "Event " << ientry << endl;

    Double_t eventweight = info->eventweight;

    //********************************************************
    // Load the branches
    //********************************************************
    jetArr->Clear(); 
    genparticleArr->Clear(); 
    genjetArr->Clear(); 
    jetBr->GetEntry(ientry);
    genparticleBr->GetEntry(ientry);
    genjetBr->GetEntry(ientry);


    //********************************************************
    // Pileup Energy Density
    //********************************************************
    Double_t rho = info->RhoKt6PFJets;
    UInt_t EAEra = kDataEra_2012_Data;


    //********************************************************
    // Loop over jets
    //********************************************************
    for(Int_t i=0; i<jetArr->GetEntries(); i++) {
      
      const cmsana::TJet *rawjet = (cmsana::TJet*)((*jetArr)[i]);
      double JEC = JetEnergyCorrectionFactor(rawjet, JetCorrector, rho, false);      
      cmsana::TJet* jet = (cmsana::TJet*)rawjet->Clone(); //new object created
      jet->pt = JEC * rawjet->rawPt;
      jet->mass = JEC * rawjet->mass;

      //********************************************************
      // See if photon matches to real photon
      //********************************************************
      double minDR = 9999;
      const cmsana::TGenParticle *matchedGenPhoton = 0;
      const cmsana::TGenParticle *matchedGenLepton = 0;
      for(Int_t k=0; k<genparticleArr->GetEntries(); k++) {
        const cmsana::TGenParticle *gen = (cmsana::TGenParticle*)((*genparticleArr)[k]);
          
        //status 1 match
        if (abs(gen->pdgid) == 22 && (gen->status == 1)
            && ( gen->motherPdgID == 25 || abs(gen->motherPdgID) <= 6 || 
                 (abs(gen->motherPdgID) >= 11 && abs(gen->motherPdgID) <= 16) ||
                 abs(gen->motherPdgID) == 23 || abs(gen->motherPdgID) == 24 || abs(gen->motherPdgID) == 21
              )
          ) {
          double tmpDR = cmsana::deltaR( jet->eta, jet->phi, gen->eta , gen->phi);
          if (tmpDR < minDR && tmpDR < 0.1) {
            minDR = tmpDR;
            matchedGenPhoton = gen;
          }
        }

        if ( abs(gen->pdgid) >= 11 && abs(gen->pdgid) <= 16 && gen->status == 1) {
          if (cmsana::deltaR(jet->eta,jet->phi,gen->eta,gen->phi) < 0.1) {
            matchedGenLepton = gen;
          }
        }
      }

      const cmsana::TGenJet *matchedGenJet = 0;
      for(Int_t k=0; k<genjetArr->GetEntries(); k++) {
        const cmsana::TGenJet *genjet = (cmsana::TGenJet*)((*genjetArr)[k]);
        if (cmsana::deltaR(genjet->eta,genjet->phi,jet->eta,jet->phi) < 0.5) {
          matchedGenJet = genjet;
        }
        if (matchedGenJet) break;
      }

      //veto denominators matching to gen photons or leptons
      if (matchedGenPhoton || matchedGenLepton) {
        continue;
      }
      
      //Fill Jet
      NJetsFilled++;
      FillJetTree( jetTree, jet, matchedGenJet, rho, EAEra, info, 1.0);
      
      delete jet; //memory clean up
    } //loop over jets
    
  } //loop over events


  cout << "Total Jets: " << NJetsFilled << endl;


  delete info;
  delete jetArr;
  delete genparticleArr;


  outputFile->Write();
  outputFile->Close();
  
  delete outputFile;
  if (jetTree) delete jetTree;

  gBenchmark->Show("WWTemplate");       
} 



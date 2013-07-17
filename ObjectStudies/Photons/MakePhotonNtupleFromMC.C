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
#include "CMSAna/ObjectStudies/interface/PhotonTree.h"
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
void MakeNtuple(const string inputFilename,  const string outputFilename, bool selectRealPhotons, int firstEvent, int processNEvents);

void MakeNtuple(const string inputFilename,  const string outputFilename, bool selectRealPhotons) {
  MakeNtuple(inputFilename, outputFilename, selectRealPhotons, -1,-1);
}

//=== MAIN MACRO =================================================================================================
void MakePhotonNtupleFromMC(Int_t Sample = 0) {
 
  //Phase 1 Samples: Age0
  if (Sample == 1) {
    MakeNtuple("root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-04/HHtoBBGG-14tev-START53_V7A/BACONNtuple.dihiggs-bbgg-14tev.1.root", "photonNtuple.summer12.real.root", true, -1,-1);
  } else if (Sample==2) {
    MakeNtuple("root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-04/ttHbb-125-START53_V7A/BACONNtuple_15_1_9MJ.root", "photonNtuple.summer12.fake.root", false, -1, -1);
  } else if (Sample == 11) {
    MakeNtuple("root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-04_SLHC4/WZHgg-125-Age0_STAR17_61_V1A/BACONNtuple_14_1_FoQ.root", "photonNtuple.highpileup.real.root", true, 100, 0);
  } else if (Sample == 12) {
    MakeNtuple("root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-04_SLHC4/ttHbb-125-Age0_STAR17_61_V1A/BACONNtuple_17_1_a9l.root", "photonNtuple.highpileup.fake2.root", false, 100, 0);
  }


} 

void MakePhotonNtupleFromMC(const string inputFilename, const string outputFilename, bool selectRealPhotons, int splitIntoNJobs, int jobNumber) {
  MakeNtuple(inputFilename,outputFilename,selectRealPhotons, splitIntoNJobs, jobNumber);
}

void MakePhotonNtupleFromMC(const string inputFilename, const string outputFilename, bool selectRealPhotons) {
  MakeNtuple(inputFilename,outputFilename,selectRealPhotons);
}



void MakeNtuple(const string inputFilename, const string outputFilename, bool selectRealPhotons, int splitIntoNJobs, int jobNumber)
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
  cmsana::PhotonTree *phoTree = new cmsana::PhotonTree;
  phoTree->CreateTree();
  phoTree->tree_->SetAutoFlush(0);

  UInt_t NPhotonsFilled = 0;
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  cmsana::TEventInfo *info    = new cmsana::TEventInfo();
  TClonesArray *photonArr = new TClonesArray("cmsana::TPhoton");
  TClonesArray *pfcandidateArr = new TClonesArray("cmsana::TPFCandidate");
  TClonesArray *genparticleArr = new TClonesArray("cmsana::TGenParticle");
  TClonesArray *genjetArr = new TClonesArray("cmsana::TGenJet");

  Int_t NEvents = 0;



  //********************************************************
  // Get Tree
  //********************************************************
  eventTree = getTreeFromFile(inputFilename.c_str(),"Events"); 
  TBranch *infoBr;
  TBranch *photonBr;
  TBranch *pfcandidateBr;
  TBranch *genparticleBr;
  TBranch *genjetBr;
  
  
  //*****************************************************************************************
  //Loop over Data Tree
  //*****************************************************************************************
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Photon", &photonArr); photonBr = eventTree->GetBranch("Photon");
  eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr); pfcandidateBr = eventTree->GetBranch("PFCandidate");
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
    photonArr->Clear(); 
    pfcandidateArr->Clear(); 
    genparticleArr->Clear(); 
    genjetArr->Clear(); 
    photonBr->GetEntry(ientry);
    pfcandidateBr->GetEntry(ientry);
    genparticleBr->GetEntry(ientry);
    genjetBr->GetEntry(ientry);


    //********************************************************
    // Pileup Energy Density
    //********************************************************
    Double_t rho = info->RhoKt6PFJets;
    UInt_t EAEra = kDataEra_2012_Data;

    //********************************************************
    // For Real Photons Start with Gen Photons
    //********************************************************
    if (selectRealPhotons) {
      for(Int_t k=0; k<genparticleArr->GetEntries(); k++) {
        const cmsana::TGenParticle *gen = (cmsana::TGenParticle*)((*genparticleArr)[k]);
        
        //status 1 match
        if (abs(gen->pdgid) == 22 && (gen->status == 1)
            && ( gen->motherPdgID == 25 || abs(gen->motherPdgID) <= 6 || 
                 (abs(gen->motherPdgID) >= 11 && abs(gen->motherPdgID) <= 14) ||
                 abs(gen->motherPdgID) == 23 || abs(gen->motherPdgID) == 24 || abs(gen->motherPdgID) == 21
              )
          ) {
          
          double minDR = 9999;
          const cmsana::TPhoton *pho = 0;
          //Find Matching Reco Photon
          for(Int_t i=0; i<photonArr->GetEntries(); i++) {
            const cmsana::TPhoton *tmppho = (cmsana::TPhoton*)((*photonArr)[i]);
            double tmpDR = cmsana::deltaR( tmppho->eta, tmppho->phi, gen->eta , gen->phi);
            if (tmpDR < minDR && tmpDR < 0.1) {
              minDR = tmpDR;
              pho = tmppho;
            }                        
          }

          //Fill Photon
          NPhotonsFilled++;
          FillPhotonTree( phoTree, gen, pho, pfcandidateArr, genparticleArr, rho, EAEra, 
                          info, 1.0);
          
        } //found proper gen photon
      } //loop over gen particles
    } else {

      //********************************************************
      // For Fake Photons Start with Reco Photons
      //********************************************************
      for(Int_t i=0; i<photonArr->GetEntries(); i++) {
        const cmsana::TPhoton *pho = (cmsana::TPhoton*)((*photonArr)[i]);
        
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
            double tmpDR = cmsana::deltaR( pho->eta, pho->phi, gen->eta , gen->phi);
            if (tmpDR < minDR && tmpDR < 0.1) {
              minDR = tmpDR;
              matchedGenPhoton = gen;
            }
          }

          if ( abs(gen->pdgid) >= 11 && abs(gen->pdgid) <= 16 && gen->status == 1) {
            if (cmsana::deltaR(pho->eta,pho->phi,gen->eta,gen->phi) < 0.1) {
              matchedGenLepton = gen;
            }
          }
        }
        
        //genjets
        const cmsana::TGenJet *matchedGenJet = 0;
        for(Int_t k=0; k<genjetArr->GetEntries(); k++) {
          const cmsana::TGenJet *genjet = (cmsana::TGenJet*)((*genjetArr)[k]);
          if (cmsana::deltaR(genjet->eta,genjet->phi,pho->eta,pho->phi) < 0.5) {
            matchedGenJet = genjet;
          }
        }

        //veto denominators matching to gen photons or leptons
        if (matchedGenPhoton || matchedGenLepton) {
          continue;
        }

        //veto denominators not matching to any genjets in order to avoid photons from pileup
        if (!matchedGenJet) {
          continue;
        }               
        
        //Fill Photon
        NPhotonsFilled++;
        FillPhotonTree( phoTree, 0, pho, pfcandidateArr, genparticleArr, rho, EAEra, 
                        info, 1.0);
        
      } //loop over photons      
    } // for fake photons
    


  }

  cout << "Total Photons: " << NPhotonsFilled << endl;


  delete info;
  delete photonArr;
  delete pfcandidateArr;
  delete genparticleArr;


  outputFile->Write();
  outputFile->Close();
  
  delete outputFile;
  if (phoTree) delete phoTree;

  gBenchmark->Show("WWTemplate");       
} 



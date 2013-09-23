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

#include "CMSAna/ObjectStudies/interface/MuonTree.h"
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
void MakeHZZMuonNtuple(Int_t Sample = 100) {

  //default
  if (Sample == 100) {
    MakeNtuple("/data1/sixie/ntuples/BACON/HZZ4L/V1/BACONNtuple_100.root","MuonNtuple.HZZ4L.100.root");
  }
} 


void MakeHZZMuonNtuple(const string inputFilename, const string outputFilename) {
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
  cmsana::MuonTree *muTree = new cmsana::MuonTree;
  muTree->CreateTree(cmsana::MuonTree::kMuTreeLight);
  muTree->tree_->SetAutoFlush(0);

  UInt_t NMuonsFilled = 0;
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  cmsana::TEventInfo *info    = new cmsana::TEventInfo();
  TClonesArray *muonArr = new TClonesArray("cmsana::TMuon");
  TClonesArray *genparticleArr = new TClonesArray("cmsana::TGenParticle");

  Int_t NEvents = 0;

  //********************************************************
  // Get Tree
  //********************************************************
  eventTree = getTreeFromFile(inputFilename.c_str(),"Events"); 
  TBranch *infoBr;
  TBranch *muonBr;
  TBranch *genparticleBr;


  //*****************************************************************************************
  //Loop over Data Tree
  //*****************************************************************************************
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("GenParticle", &genparticleArr); genparticleBr = eventTree->GetBranch("GenParticle");

  cout << "InputFile " << inputFilename << " --- Total Events : " << eventTree->GetEntries() << endl;
  for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
    infoBr->GetEntry(ientry);
		
    if (ientry % 10000 == 0) cout << "Event " << ientry << endl;

    Double_t eventweight = info->eventweight;

    //********************************************************
    // Load the branches
    //********************************************************
    muonArr->Clear(); 
    genparticleArr->Clear(); 
    muonBr->GetEntry(ientry);
    genparticleBr->GetEntry(ientry);

    //********************************************************
    // Pileup Energy Density
    //********************************************************
    Double_t rho = info->RhoKt6PFJets;


    //********************************************************
    // Loop over Gen muons
    //********************************************************
    vector<float> muUsedEta;
    vector<float> muUsedPhi;

    for(Int_t i=0; i<genparticleArr->GetEntries(); i++) {
      const cmsana::TGenParticle *gen = (cmsana::TGenParticle*)((*genparticleArr)[i]);
      
      if (abs(gen->pdgid) == 13 && gen->status == 1 && (abs(gen->motherPdgID) == 23 || gen->motherPdgID == gen->pdgid)) {

        //***************************************
        //Don't use duplicate muons
        //***************************************
        bool isDuplicate = false;
        for(Int_t j=0; j<muUsedEta.size(); j++) {
          if (cmsana::deltaR(gen->eta,gen->phi,muUsedEta[j],muUsedPhi[j]) < 0.1) {
            isDuplicate = true;
            break;
          }
        }
        if (isDuplicate) continue;

        //***************************************
        //Find corresponding reco muon
        //***************************************
        double mupt = -999;
        double mueta = -999;
        double muphi = -999;

        for(Int_t j=0; j<muonArr->GetEntries(); j++) {
          const cmsana::TMuon *mu = (cmsana::TMuon*)((*muonArr)[j]);
          if (cmsana::deltaR(gen->eta,gen->phi,mu->eta,mu->phi) < 0.1) {
            mupt = 0;
            if (mu->pt > 0) mupt = mu->pt; //protect against possible weird negative pt values
            mueta = mu->eta;
            muphi = mu->phi;
          }
        }

        //***************************************
        //Fill gen muon
        //***************************************
        muTree->fWeight = eventweight;
        muTree->fRunNumber = info->runNum;
        muTree->fLumiSectionNumber = info->lumiSec;
        muTree->fEventNumber = info->evtNum;
        muTree->fMuEventNumberParity = (info->evtNum % 2 == 0);
        muTree->fRho = rho;
        muTree->fMuGenPt = gen->pt; 
        muTree->fMuGenEta = gen->eta; 
        muTree->fMuGenPhi = gen->phi; 
        muTree->fMuPt = mupt;
        muTree->fMuEta = mueta;
        muTree->fMuPhi = muphi;
        muTree->tree_->Fill();
        NMuonsFilled++;

      } // if gen mu is found

    } //loop over gen particles

  }

  delete info;
  delete muonArr;
  delete genparticleArr;

  cout << "Total Muons: " << NMuonsFilled << endl;

  outputFile->Write();
  outputFile->Close();
  
  delete outputFile;
  if (muTree) delete muTree;

  gBenchmark->Show("WWTemplate");       
} 



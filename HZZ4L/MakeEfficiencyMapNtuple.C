
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

// output data structs
#include "CMSAna/HZZ4L/interface/HZZEfficiencyMap.hh"
#include "CMSAna/ObjectStudies/interface/ElectronTree.h"
#include "CMSAna/ObjectStudies/interface/MuonTree.h"

#endif


//=== FUNCTION DECLARATIONS ======================================================================================

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




//=== MAIN MACRO =================================================================================================

void MakeEfficiencyMapNtuple(string EleNtupleFile, string MuonNtupleFile, string Label = "") 
{  
  string label = Label;
  if (Label != "") label = "_" + Label;


  //--------------------------------------------------------------------------------------------------------------
  // output ntuple structure
  //==============================================================================================================  
  HZZEfficiencyMap efficiencyMap;

  TFile *fOutputFile = new TFile(("HZZEfficiencyMap"+label+".root").c_str(), "RECREATE");
  TTree *fOutputTree = new TTree("EfficiencyMap","EfficiencyMap");
  fOutputTree->Branch("efficiencyMap",&efficiencyMap);


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  cmsana::ElectronTree eleTree;
  eleTree.LoadTree(EleNtupleFile.c_str());
  eleTree.InitTree(cmsana::ElectronTree::kEleTreeLight);

  cmsana::MuonTree muTree;
  muTree.LoadTree(MuonNtupleFile.c_str());
  muTree.InitTree(cmsana::MuonTree::kMuTreeLight);

  cout << "Total Electron Entries : " << eleTree.tree_->GetEntries() << "\n";
  for(UInt_t i=0; i<eleTree.tree_->GetEntries(); i++) {       	
    eleTree.tree_->GetEntry(i);
    if (i % 10000 == 0) cout << "Ele " << i << endl;
	
    //Fill Information
    efficiencyMap.type = 11;
    efficiencyMap.genpt = eleTree.fEleGenPt;
    efficiencyMap.geneta = eleTree.fEleGenEta;
    efficiencyMap.genphi = eleTree.fEleGenPhi;
    efficiencyMap.recopt = eleTree.fElePt;
    efficiencyMap.recoeta = eleTree.fEleEta;
    efficiencyMap.recophi = eleTree.fElePhi;
    efficiencyMap.recocharge = 1;
    efficiencyMap.pass = bool( eleTree.fElePt > 0 ) ;
    efficiencyMap.nvtx = eleTree.fNVertices;
    efficiencyMap.weight = eleTree.fWeight;

    //*********************************************************
    //Fill 
    //*********************************************************
    fOutputTree->Fill();    
  }
  
  cout << "Total Muon Entries : " << muTree.tree_->GetEntries() << "\n";
  for(UInt_t i=0; i<muTree.tree_->GetEntries(); i++) {       	
    muTree.tree_->GetEntry(i);
    if (i % 10000 == 0) cout << "Mu " << i << endl;
	
    //Fill Information
    efficiencyMap.type = 13;
    efficiencyMap.genpt = muTree.fMuGenPt;
    efficiencyMap.geneta = muTree.fMuGenEta;
    efficiencyMap.genphi = muTree.fMuGenPhi;
    efficiencyMap.recopt = muTree.fMuPt;
    efficiencyMap.recoeta = muTree.fMuEta;
    efficiencyMap.recophi = muTree.fMuPhi;
    efficiencyMap.recocharge = 1;
    efficiencyMap.pass = bool( muTree.fMuPt > 0 ) ;
    efficiencyMap.nvtx = muTree.fNVertices;
    efficiencyMap.weight = muTree.fWeight;

    //*********************************************************
    //Fill 
    //*********************************************************
    fOutputTree->Fill();    
  }
  

  fOutputFile->cd();
  fOutputFile->Write();
  fOutputFile->Close();
  delete fOutputFile;

   
} 




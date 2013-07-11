// *************************************
// From EOS
//
// root -l -b -q CMSAna/Skimming/SkimDoubleBTagPerFile.C+\(\"root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-00/diphjets-START53_V7A/BACONNtuple_42_1_D1l.root\",\"testSkim.root\",true\)
// *************************************


#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TVector3.h>               
#include <TChain.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h> 
#include <TBranch.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TBenchmark.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <TMath.h>
#endif

// define structures to read in ntuple
#include "CMSAna/DataTree/interface/DataTreeDefs.hh"
#include "CMSAna/DataTree/interface/TEventInfo.hh"
#include "CMSAna/DataTree/interface/TGenParticle.hh"
#include "CMSAna/DataTree/interface/TGenJet.hh"
#include "CMSAna/DataTree/interface/TMuon.hh"
#include "CMSAna/DataTree/interface/TElectron.hh"
#include "CMSAna/DataTree/interface/TPhoton.hh"
#include "CMSAna/DataTree/interface/TJet.hh"
#include "CMSAna/DataTree/interface/TPFCandidate.hh"
#include "CMSAna/DataTree/interface/TVertex.hh"


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


// Main macro function
//--------------------------------------------------------------------------------------------------
void SkimDoubleBTagPerFile(string inputFilename, string outputFilename, bool isMC = true) 
{
  gBenchmark->Start("SkimNtuples");
    
  TTree::SetMaxTreeSize(kMaxLong64);
  
  // Don't write TObject part of the objects
  cmsana::TEventInfo::Class()->IgnoreTObjectStreamer();
  cmsana::TGenParticle::Class()->IgnoreTObjectStreamer();
  cmsana::TGenJet::Class()->IgnoreTObjectStreamer();
  cmsana::TElectron::Class()->IgnoreTObjectStreamer();
  cmsana::TMuon::Class()->IgnoreTObjectStreamer();
  cmsana::TJet::Class()->IgnoreTObjectStreamer();
  cmsana::TPhoton::Class()->IgnoreTObjectStreamer();
  cmsana::TPFCandidate::Class()->IgnoreTObjectStreamer();
  cmsana::TVertex::Class()->IgnoreTObjectStreamer();

  // Data structures to store info from TTrees
  cmsana::TEventInfo *info      = new cmsana::TEventInfo();
  TClonesArray *genparticleArr  = new TClonesArray("cmsana::TGenParticle");
  TClonesArray *genjetArr       = new TClonesArray("cmsana::TGenJet");
  TClonesArray *muonArr         = new TClonesArray("cmsana::TMuon");
  TClonesArray *electronArr     = new TClonesArray("cmsana::TElectron");
  TClonesArray *photonArr       = new TClonesArray("cmsana::TPhoton");
  TClonesArray *jetArr          = new TClonesArray("cmsana::TJet");
  TClonesArray *pfcandidateArr  = new TClonesArray("cmsana::TPFCandidate",20000);
  TClonesArray *vertexArr       = new TClonesArray("cmsana::TVertex");
 
  UInt_t nInputEvts = 0;
  UInt_t nPassEvts  = 0;
  UInt_t nEventsTotal = 0;

  TFile* outfile = new TFile(outputFilename.c_str(), "RECREATE");
  
  //
  // Initialize data trees and structs
  // 
  TTree *outEventTree = new TTree("Events","Events"); 
  outEventTree->Branch("Info",     &info);
  if (isMC) {
    outEventTree->Branch("GenParticle",  &genparticleArr);
    outEventTree->Branch("GenJet",       &genjetArr);
  }
  outEventTree->Branch("Muon",          &muonArr);
  outEventTree->Branch("Electron",      &electronArr);
  outEventTree->Branch("Photon",        &photonArr);
  outEventTree->Branch("PFJet",         &jetArr);
  outEventTree->Branch("PFCandidate",   &pfcandidateArr);
  outEventTree->Branch("Vertex",        &vertexArr);

  // list input ntuple files to be skimmed
  
  cout << "Skimming " << inputFilename << endl;
  TTree *eventTree = getTreeFromFile(inputFilename.c_str(),"Events"); 
  nEventsTotal = eventTree->GetEntries();
  assert(eventTree);
    
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",     &info);                TBranch *infoBr     = eventTree->GetBranch("Info");
  TBranch *genparticleBr = 0;
  TBranch *genjetBr      = 0;
  if (isMC) {
    eventTree->SetBranchAddress("GenParticle", &genparticleArr); genparticleBr       = eventTree->GetBranch("GenParticle");
    eventTree->SetBranchAddress("GenJet", &genjetArr);           genjetBr            = eventTree->GetBranch("GenJet");
  }
  eventTree->SetBranchAddress("Muon",     &muonArr);             TBranch *muonBr     = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("Electron", &electronArr);         TBranch *electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("Photon",    &photonArr);          TBranch *photonBr   = eventTree->GetBranch("Photon");
  eventTree->SetBranchAddress("PFJet",    &jetArr);              TBranch *jetBr      = eventTree->GetBranch("PFJet");
  eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr);   TBranch *pfcandidateBr   = eventTree->GetBranch("PFCandidate");
  eventTree->SetBranchAddress("Vertex", &vertexArr);             TBranch *vertexBr   = eventTree->GetBranch("Vertex");
     
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) { 
    infoBr->GetEntry(ientry);
    if (isMC) {
      genparticleArr->Clear(); genparticleBr->GetEntry(ientry);
      genjetArr->Clear(); genjetBr->GetEntry(ientry);
    }
    muonArr->Clear();     muonBr->GetEntry(ientry);
    electronArr->Clear(); electronBr->GetEntry(ientry);      
    photonArr->Clear(); photonBr->GetEntry(ientry);
    jetArr->Clear(); jetBr->GetEntry(ientry);
    pfcandidateArr->Clear(); pfcandidateBr->GetEntry(ientry);
    vertexArr->Clear(); vertexBr->GetEntry(ientry);

    if (ientry % 100000 == 0) cout << "Events: " << ientry << endl;

    nInputEvts++;
      
    bool keep = false;

    //2 b-tagged jets
    Int_t NBTags = 0;
    for(Int_t i=0; i<jetArr->GetEntriesFast(); i++) {

      const cmsana::TJet *jet = (cmsana::TJet*)((*jetArr)[i]);
      if (!(jet->pt > 25)) continue;
      if (!(fabs(jet->eta) < 2.4)) continue;
      if (!(jet->CombinedSecondaryVertexBJetTagsDisc > 0.679)) continue;
      NBTags++;
    }
    if (NBTags >= 2) keep = true;

    if(keep) {
      outEventTree->Fill();
      nPassEvts++;
    }
  }
  outfile->Write();
  outfile->Close();
      
  std::cout << outputFilename << " created!" << std::endl;
  std::cout << " >>> Total Number of Events: " << nEventsTotal << std::endl;
  std::cout << " >>> Events processed: " << nInputEvts << std::endl;
  std::cout << " >>>   Events passing: " << nPassEvts << std::endl;

  gBenchmark->Show("SkimNtuples");
}  


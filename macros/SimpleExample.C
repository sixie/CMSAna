//================================================================================================
//
// Simple Example
//
//root -l -b -q  CMSAna/macros/SimpleExample.C+'("root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-03_SLHC4/DYToEEAge0START_STAR17_61_V1A/BACONNtuple_1_1_Pmm.root")'
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <TH1F.h>                
#include <TCanvas.h>                

// define structures to read in ntuple
#include "CMSAna/DataTree/interface/Types.h"
#include "CMSAna/DataTree/interface/TEventInfo.hh"
#include "CMSAna/DataTree/interface/TGenParticle.hh"
#include "CMSAna/DataTree/interface/TGenJet.hh"
#include "CMSAna/DataTree/interface/TPhoton.hh"
#include "CMSAna/DataTree/interface/TMuon.hh"
#include "CMSAna/DataTree/interface/TElectron.hh"
#include "CMSAna/DataTree/interface/TJet.hh"
#include "CMSAna/DataTree/interface/TPFCandidate.hh"


// Helper functions for Photon selection
#include "CMSAna/Utils/PhotonID.hh"
#include "CMSAna/Utils/LeptonID.hh"
#include "CMSAna/Utils/CommonDefs.hh"
#include "CMSAna/Utils/CommonTools.hh"

// For Jet energy corrections
#include "CMSAna/JetEnergyCorrections/interface/FactorizedJetCorrector.h"
#include "CMSAna/JetEnergyCorrections/interface/JetCorrectorParameters.h"
#include "CMSAna/Utils/JetEnergyCorrections.hh"

#endif


//=== MAIN MACRO ================================================================================================= 

void SimpleExample(const string inputfile) {
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;


  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  TH1F *histElectronPt = new TH1F ("histElectronPt",";Electron p_{T} [GeV/c^{2}];Number of Events", 100, 0 , 100);


  //*****************************************************************************************
  // Set up input
  //*****************************************************************************************
  TFile *infile=0;
  TTree *eventTree=0;
  
  // Data structures to store info from TTrees
  cmsana::TEventInfo *info  = new cmsana::TEventInfo();
  TClonesArray *genparticleArr = new TClonesArray("cmsana::TGenParticle");
  TClonesArray *genjetArr = new TClonesArray("cmsana::TGenJet");
  TClonesArray *photonArr = new TClonesArray("cmsana::TPhoton");
  TClonesArray *muonArr = new TClonesArray("cmsana::TMuon");
  TClonesArray *electronArr = new TClonesArray("cmsana::TElectron");
  TClonesArray *jetArr = new TClonesArray("cmsana::TJet");
  TClonesArray *pfcandidateArr = new TClonesArray("cmsana::TPFCandidate");

  // Read input file and get the TTrees
  cout << "Processing " << inputfile << "..." << endl;
  infile = TFile::Open(inputfile.c_str(),"read");
  assert(infile);

    
  eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
  eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("GenParticle", &genparticleArr); TBranch *genparticleBr = eventTree->GetBranch("GenParticle");
  eventTree->SetBranchAddress("GenJet", &genjetArr); TBranch *genjetBr = eventTree->GetBranch("GenJet");
  eventTree->SetBranchAddress("Photon", &photonArr); TBranch *photonBr = eventTree->GetBranch("Photon");
  eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("PFJet", &jetArr); TBranch *jetBr = eventTree->GetBranch("PFJet");
  eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr);  TBranch *pfcandidateBr = eventTree->GetBranch("PFCandidate");
  cout << "NEvents = " << eventTree->GetEntries() << endl;

  
  Double_t weight = 1;

  // null vector for default four vectors
  cmsana::FourVector null(0.0,0.0,0.0,0.0);

  // loop over events
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    if (ientry % 1000 == 0) cout << "Processed Event " << ientry << endl;
    infoBr->GetEntry(ientry);

    //***********************************************************
    // Definition of Pileup Energy density
    //***********************************************************
    Double_t rho = info->RhoKt6PFJets;

    genparticleArr->Clear(); 
    genjetArr->Clear(); 
    photonArr->Clear();
    jetArr->Clear(); 
    muonArr->Clear();
    electronArr->Clear();
    pfcandidateArr->Clear(); 

    genparticleBr->GetEntry(ientry);
    genjetBr->GetEntry(ientry);
    photonBr->GetEntry(ientry);
    jetBr->GetEntry(ientry);
    muonBr->GetEntry(ientry);
    electronBr->GetEntry(ientry);
    pfcandidateBr->GetEntry(ientry);

    //***********************************************************
    // Gen-Level particles
    //***********************************************************
    for(Int_t i=0; i<genparticleArr->GetEntriesFast(); i++) {
      const cmsana::TGenParticle *p = (cmsana::TGenParticle*)((*genparticleArr)[i]);

//       cout << p->pdgid << " " << p->status << " " << p->pt << " " << p->eta << " " << p->phi 
//            << " | " << p->motherPdgID << "\n";

    }


    //********************************************************
    //Loop Over Electrons
    //********************************************************
    for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
      const cmsana::TElectron *ele = (cmsana::TElectron*)((*electronArr)[i]);
      
      if (!(ele->pt > 20)) continue;
      if (!(fabs(ele->scEta) < 2.5)) continue;
      
      //pass loose veto cuts
      //if (!PassEleSimpleCutsVeto( ele, pfcandidateArr, rho,
      //                            kDataEra_2012_MC, false)) continue;

      //cout << "Electrons : " << ele->pt << " " << ele->eta << " " << ele->phi << "\n";
      histElectronPt->Fill(ele->pt);

    }


    //********************************************************
    //Loop Over Muons
    //********************************************************
    for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
      const cmsana::TMuon *mu = (cmsana::TMuon*)((*muonArr)[i]);

      if (!(mu->pt > 10)) continue;
      if (!(fabs(mu->eta) < 2.4)) continue;

      //pass loose veto cuts
      if (!PassMuonIDVeto(mu)) continue;
      if (!(ComputeMuonPFIsoRings(mu,pfcandidateArr, rho,
                                  kDataEra_2012_MC, kPFIso, 0.0, 0.3, false) / mu->pt < 0.2)) continue;

      //cout << "Muons: " << mu->pt << " " << mu->eta << " " << mu->phi << "\n";

    }

  }

  delete infile;
  infile=0, eventTree=0;      
  delete info;
  delete photonArr;
  delete jetArr;
  delete genparticleArr;
  delete genjetArr;
  
     
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  TCanvas *cv = new TCanvas("cv","cv", 800, 600);
  histElectronPt->Draw();
  cv->SaveAs("SimpleExample_ElectronPt.gif");
        

}

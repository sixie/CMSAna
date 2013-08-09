//================================================================================================
//
// Select HH->bbgammagamma
//
//root -l -b -q  CMSAna/HHToBBGG/HHToBBGGSelectionMistags.C+'("root://eoscms//eos/cms//store/group/phys_higgs/future/sixie/Madgraph/DiPhotonJJ_M60To200_14TeV-v2/BACON/BACONNtuple_GenOnly_DiphotonJJ_14TeV_16.root","HHToBBGGNtuple.DiPhotonJJ_M60To200_14TeV.16.root",7)'
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
#include <TH2F.h>                

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

#include "CMSAna/HHToBBGG/interface/HHToBBGGEventTree.hh"

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

void HHToBBGGSelectionMistags(const string inputfile,          // input file
                       const string outputfile,         // output directory
                       Int_t SampleType = 1) {
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  Double_t nEvents = 0;
  Double_t nEventsTwoRealPho = 0;
  
  //*****************************************************************************************
  // Set up output ntuple
  //*****************************************************************************************
  TFile *outFile = new TFile(outputfile.c_str(),"RECREATE"); 
  TH1F *NEvents =  new TH1F("NEvents",";;",1,-0.5,0.5);
  //edit
  TH1F *NEventsTwoRealPho =  new TH1F("NEventsTwoRealPho",";;",1,-0.5,0.5);
  TH1F *NPUMean = new TH1F("NPUMean",";NPUMean;Number of Events", 100, -0.5, 99.5);
  cmsana::HHToBBGGEventTree *outputEventTree = new cmsana::HHToBBGGEventTree;
  outputEventTree->CreateTree();


  //*****************************************************************************************
  // Set up input
  //*****************************************************************************************
  TFile *infile=0;
  TTree *eventTree=0;
  
  // Data structures to store info from TTrees
  cmsana::TEventInfo *info  = new cmsana::TEventInfo();
  TClonesArray *genparticleArr = new TClonesArray("cmsana::TGenParticle");
  TClonesArray *genjetArr = new TClonesArray("cmsana::TGenJet");

  TFile *BJetMistagRateLightJetsFile = TFile::Open("/afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/BTaggingEfficiency_LightJetsMistagRate.root");
  TFile *BJetMistagRateCharmJetsFile = TFile::Open("/afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/BTaggingEfficiency_CharmJetsMistagRate.root");
  TFile *PhotonEfficiencyFile = TFile::Open("/afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/PhotonEfficiency_PromptPhoton.root");
  TH2F *BJetMistagRateLightJetsHist = (TH2F*)BJetMistagRateLightJetsFile->Get("MistagRate_CSVMedium_Pt_Eta");
  TH2F *BJetMistagRateCharmJetsHist = (TH2F*)BJetMistagRateCharmJetsFile->Get("MistagRate_CSVMedium_Pt_Eta");
  TH2F *PhotonEfficiencyHist = (TH2F*)PhotonEfficiencyFile->Get("Efficiency_PtEta");
  TH2F *bjet1FakeRateHist = 0;
  TH2F *bjet2FakeRateHist = 0;
  double bjet1Eff = 0;
  double bjet2Eff = 0;

  // Read input file and get the TTrees
  cout << "Processing " << inputfile << "..." << endl;
  infile = TFile::Open(inputfile.c_str(),"read");
  assert(infile);

    
  eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
  eventTree->SetBranchAddress("Info", &info); TBranch *infoBr = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("GenParticle", &genparticleArr); TBranch *genparticleBr = eventTree->GetBranch("GenParticle");
  eventTree->SetBranchAddress("GenJet", &genjetArr); TBranch *genjetBr = eventTree->GetBranch("GenJet");
  cout << "NEvents = " << eventTree->GetEntries() << endl;

  
  Double_t weight = 1;

  // null vector for default four vectors
  cmsana::FourVector null(0.0,0.0,0.0,0.0);

  // loop over events
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    if (ientry % 1000 == 0) cout << "Processed Event " << ientry << endl;
    infoBr->GetEntry(ientry);

    NEvents->Fill(0);
    NEventsTwoRealPho->Fill(0);
    NPUMean->Fill(info->nPUMean);

    //***********************************************************
    // Definition of Pileup Energy density
    //***********************************************************
    Double_t rho = info->RhoKt6PFJets;

    genparticleArr->Clear();
    genjetArr->Clear(); 

    genparticleBr->GetEntry(ientry);
    genjetBr->GetEntry(ientry);

    double tmpHT = 0;

    //***********************************************************
    // Find Gen-Level particles
    //***********************************************************
    const cmsana::TGenParticle *genPhoton1 = 0;
    const cmsana::TGenParticle *genPhoton2 = 0;
    for(Int_t i=0; i<genparticleArr->GetEntriesFast(); i++) {
      const cmsana::TGenParticle *p = (cmsana::TGenParticle*)((*genparticleArr)[i]);

      //***********************************************************
      // Find Gen Photons
      //***********************************************************
      if ( SampleType == cmsana::HHToBBGGEventTree::GGPlusTwoMistag) {
        if (p->pdgid == 22 && ( p->motherPdgID == 21 || (abs(p->motherPdgID) >= 1 && abs(p->motherPdgID) <= 6) ) ) {
          if (!genPhoton1) {
            genPhoton1 = p;
          } else if (!genPhoton2) {
            genPhoton2 = p;
          }
        }
      }

    }
    
    //sampleType
    cmsana::HHToBBGGEventTree::SampleType stype = cmsana::HHToBBGGEventTree::none;
    if (SampleType == 7) stype = cmsana::HHToBBGGEventTree::GGPlusTwoMistag;

    outputEventTree->sampletype = stype;
    outputEventTree->run = info->runNum;
    outputEventTree->lumi = info->lumiSec;
    outputEventTree->event = info->evtNum;
    outputEventTree->npu = info->nPU;
    outputEventTree->rho = info->RhoKt6PFJets;
    outputEventTree->nvtx = info->nGoodPV;

    cmsana::FourVectorM genpho1v;
    cmsana::FourVectorM genpho2v;
    outputEventTree->genpho1 = null;
    outputEventTree->genpho2 = null;
    if (genPhoton1) {
      genpho1v.SetPt(genPhoton1->pt);
      genpho1v.SetEta(genPhoton1->eta);
      genpho1v.SetPhi(genPhoton1->phi);
      genpho1v.SetM(0);
      outputEventTree->genpho1 = genpho1v;
    }
    if (genPhoton2) {
      genpho2v.SetPt(genPhoton2->pt);
      genpho2v.SetEta(genPhoton2->eta);
      genpho2v.SetPhi(genPhoton2->phi);
      genpho2v.SetM(0);   
      outputEventTree->genpho2 = genpho2v;    
    }

    cmsana::FourVectorM genb1v;
    cmsana::FourVectorM genb2v;
    outputEventTree->genb1 = null;
    outputEventTree->genb2 = null;

    //****************************************************************
    // Loop over genjets to make collection of eligible denominators
    //****************************************************************   
    UInt_t njets = 0;
    UInt_t ncentraljets = 0;
    vector<const cmsana::TGenJet*> goodBJets;
    for(Int_t i=0; i<genjetArr->GetEntriesFast(); i++) {
      const cmsana::TGenJet *genjet = (cmsana::TGenJet*)((*genjetArr)[i]);

      bool passBJetCuts = true;
//       if ( !(genPhoton1 && genPhoton2) ) passBJetCuts = false; //Si: not sure this is useful

      //***************************************************
      //bjet cuts
      //***************************************************
      if (!(genjet->pt > 30)) passBJetCuts = false;
      if (!(fabs(genjet->eta) < 2.4)) passBJetCuts = false;
      if (passBJetCuts) {      
        //save good bjets
        goodBJets.push_back((cmsana::TGenJet*)genjet->Clone()); //new object created
      }

      //***************************************************
      //jet counting
      //***************************************************
      if (genjet->pt > 30) {
        tmpHT += genjet->pt;
        njets++;
        if (fabs(genjet->eta) < 2.5) ncentraljets++;
      }
    }


    //***************************************************
    //No real bjets in this sample
    //***************************************************
    cmsana::FourVectorM genbjet1v;
    cmsana::FourVectorM genbjet2v;
    outputEventTree->genbjet1 = null;    
    outputEventTree->genbjet2 = null;    
    
    
    //***************************************************
    //selected photons
    //***************************************************
    const cmsana::TGenParticle *photon1 = 0;
    const cmsana::TGenParticle *photon2 = 0;

    if (genPhoton1) photon1 = genPhoton1;
    if (genPhoton2) photon2 = genPhoton2;

    cmsana::FourVectorM photon1v;
    cmsana::FourVectorM photon2v;
    cmsana::FourVectorM diphotonv;
    double pho1eff = 0;
    double pho2eff = 0;
    outputEventTree->pho1 = null;      //default 4-vector
    outputEventTree->pho2 = null;
    outputEventTree->diphoton = null;    

    if (photon1) {
      photon1v.SetPt(photon1->pt);
      photon1v.SetEta(photon1->eta);
      photon1v.SetPhi(photon1->phi);
      photon1v.SetM(0);
      outputEventTree->pho1 = photon1v;
      pho1eff = PhotonEfficiencyHist->GetBinContent(PhotonEfficiencyHist->FindFixBin(fmax(fmin(photon1->pt,99.9),0.01),
                                                                    fmax(fmin(photon1->eta,2.49),-2.49)));
    } 

    if (photon2) {
      photon2v.SetPt(photon2->pt);
      photon2v.SetEta(photon2->eta);
      photon2v.SetPhi(photon2->phi);
      photon2v.SetM(0);
      outputEventTree->pho2 = photon2v;
      pho2eff = PhotonEfficiencyHist->GetBinContent(PhotonEfficiencyHist->FindFixBin(fmax(fmin(photon2->pt,99.9),0.01),
                                                                    fmax(fmin(photon2->eta,2.49),-2.49)));
    }
    
    if (photon1 && photon2) {
      diphotonv = photon1v + photon2v;
      outputEventTree->diphoton = diphotonv;
    }
    
    if (photon1) tmpHT += photon1->pt;
    if (photon2) tmpHT += photon2->pt;

    //cout << " Number of good BJets: " << goodBJets.size() << endl;

    //***************************************************************************************************************
    // now that I have all the possible goodBJets, loop through them and choose two pairwise unique goodBJets to
    // promote to bjet1 and bjet2 and store this event into the tree (don't forget to weight the event)
    //
    //***************************************************************************************************************
    //Select all pairwise goodBJets (i.e. possible bjets)
    //***************************************************************************************************************
    if (goodBJets.size() > 1) {
      for (UInt_t i=0; i<goodBJets.size()-1; i++) {
        for (UInt_t j=i+1; j<goodBJets.size(); j++) {
          const cmsana::TGenJet* bjet1 = 0;
          const cmsana::TGenJet* bjet2 = 0;
          bjet1 = goodBJets[i];
          bjet2 = goodBJets[j];
          if (bjet2->pt > bjet1->pt) {
            const cmsana::TGenJet* tmp = bjet2;
            bjet2 = bjet1;
            bjet1 = tmp;
          }

          //**************************************************
          //skip any jet pairs that have real b jets
          //**************************************************
          if ( abs(bjet1->matchedPdgId) == 5 || abs(bjet2->matchedPdgId) == 5) continue;

          //**************************************************
          //assign mistag rates for light jets and charm
          //**************************************************          
          if ( abs(bjet1->matchedPdgId) == 4) bjet1FakeRateHist = BJetMistagRateCharmJetsHist;
          else bjet1FakeRateHist = BJetMistagRateLightJetsHist;
          if ( abs(bjet2->matchedPdgId) == 4) bjet2FakeRateHist = BJetMistagRateCharmJetsHist;
          else bjet2FakeRateHist = BJetMistagRateLightJetsHist;
          bjet1Eff = bjet1FakeRateHist->GetBinContent(bjet1FakeRateHist->FindBin(fmax(fmin(bjet1->pt,119.9),0.01),fmax(fmin(bjet1->eta,2.49),-2.49)));
          bjet2Eff = bjet2FakeRateHist->GetBinContent(bjet2FakeRateHist->FindBin(fmax(fmin(bjet1->pt,119.9),0.01),fmax(fmin(bjet1->eta,2.49),-2.49)));
      
          cmsana::FourVectorM bjet1v;
          cmsana::FourVectorM bjet2v;
          cmsana::FourVectorM dibjetv;
          outputEventTree->bjet1 = null;      //default 4-vector
          outputEventTree->bjet2 = null;
          outputEventTree->dibjet = null;    
      
          if (bjet1) {
            bjet1v.SetPt(bjet1->pt);
            bjet1v.SetEta(bjet1->eta);
            bjet1v.SetPhi(bjet1->phi);
            bjet1v.SetM(bjet1->mass);
      
            outputEventTree->bjet1 = bjet1v;
          }
      
          if (bjet2) {
            bjet2v.SetPt(bjet2->pt);
            bjet2v.SetEta(bjet2->eta);
            bjet2v.SetPhi(bjet2->phi);
            bjet2v.SetM(bjet2->mass);
            outputEventTree->bjet2 = bjet2v;
          }

          if (bjet1 && bjet2) {
            dibjetv = bjet1v + bjet2v;
            outputEventTree->dibjet = dibjetv;
          }

          //********************************************************
          //bbgg system
          //********************************************************    
          cmsana::FourVectorM bbggSystemv;
          outputEventTree->bbgg = null;
          if (bjet1 && bjet2 && photon1 && photon2) {
            bbggSystemv = (photon1v + photon2v + bjet1v + bjet2v);
            outputEventTree->bbgg = bbggSystemv;
          }
      
      
      
          //********************************************************
          //NJets
          //********************************************************    
          outputEventTree->njets = njets;
          outputEventTree->ncentraljets = ncentraljets;      
          outputEventTree->nlep = 0;

          //********************************************************
          //Some kinematic variables
          //********************************************************
          outputEventTree->DRgg = -1;
          outputEventTree->DRbb = -1;
          outputEventTree->minDRgb = -1;
          if (photon1 && photon2 && bjet1 && bjet2 ) {
            outputEventTree->DRgg = cmsana::deltaR(photon1->eta, photon1->phi, photon2->eta, photon2->phi);
            outputEventTree->DRbb = cmsana::deltaR(bjet1->eta, bjet1->phi, bjet2->eta, bjet2->phi);
            outputEventTree->minDRgb = fmin(fmin(fmin( cmsana::deltaR(photon1->eta, photon1->phi, bjet1->eta, bjet1->phi), 
                                                       cmsana::deltaR(photon1->eta, photon1->phi, bjet2->eta, bjet2->phi)),
                                                 cmsana::deltaR(photon2->eta, photon2->phi, bjet1->eta, bjet1->phi)),
                                            cmsana::deltaR(photon2->eta, photon2->phi, bjet2->eta, bjet2->phi));      
          }
      
          outputEventTree->pfmet = 0;
          outputEventTree->pfTrackMET = 0;          
          outputEventTree->HT = tmpHT;
          
          //Note: currently not computing HT. we may add it later.

          //********************************************************
          //Fill Output Tree
          //********************************************************
          outputEventTree->weight = bjet1Eff * bjet2Eff * pho1eff * pho2eff;
          outputEventTree->tree_->Fill();
          nEventsTwoRealPho++;
        }
      }
    }


    //********************************************************
    //Clean up Memory
    //********************************************************
    for (uint k=0; k<goodBJets.size(); ++k) {
      if (goodBJets[k]) delete goodBJets[k];
    }
    nEvents++;

  } //loop over all events


  delete infile;
  infile=0, eventTree=0;      
  delete info;
  delete genparticleArr;
  delete genjetArr;
  
     
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  
  outFile->Write();
  outFile->Close();
  delete outFile;
  
  cout << " Number of events selected: " << nEvents << endl;
  cout << endl;
  cout << " Number of events with two real photons: " << nEventsTwoRealPho << endl;
  cout << endl;
  cout << "  <> Output saved in " << outputfile << endl;    
  cout << endl;  
      

}

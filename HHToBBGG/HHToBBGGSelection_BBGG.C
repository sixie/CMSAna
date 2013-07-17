//================================================================================================
//
// Select HH->bbgammagamma for Madgraph DiphotonBB sample 
//
//root -l -b -q  CMSAna/HHToBBGG/HHToBBGGSelection_BBGG.C+'("root://eoscms//eos/cms/store/group/phys_higgs/future/sixie/Madgraph/DiPhotonBB_M60To200_14TeV-v2/BACON/BACONNtuple_GenOnly_DiphotonBB_14TeV_89.root","HHToBBGG.BBGG.root",6)'
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

void HHToBBGGSelection_BBGG(const string inputfile,          // input file
                            const string outputfile,         // output directory
                            Int_t SampleType = 6
  ) {
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;


  //*****************************************************************************************
  //Setup
  //*****************************************************************************************
  TFile *photonEffFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/PhotonEfficiency_PromptPhoton.root","READ"); 
  TH2F* photonEffHist = (TH2F*)photonEffFile->Get("Efficiency_PtEta");
  photonEffHist->SetDirectory(0);
  photonEffFile->Close(); delete photonEffFile;

  TFile *btagEffFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/BTaggingEfficiency_BJetEfficiency.root","READ"); 
  TH2F* btagEffHist = (TH2F*)btagEffFile->Get("Efficiency_PtEta");
  btagEffHist->SetDirectory(0);
  btagEffFile->Close(); delete btagEffFile;


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  Double_t nEvents = 0;

  
  //*****************************************************************************************
  // Set up output ntuple
  //*****************************************************************************************
  TFile *outFile = new TFile(outputfile.c_str(),"RECREATE"); 
  TH1F *NEvents =  new TH1F("NEvents",";;",1,-0.5,0.5);
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

  // Read input file and get the TTrees
  cout << "Processing " << inputfile << "..." << endl;
  infile = TFile::Open(inputfile.c_str(),"read");
  assert(infile);

    
  eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
  eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("GenParticle", &genparticleArr); TBranch *genparticleBr = eventTree->GetBranch("GenParticle");
  eventTree->SetBranchAddress("GenJet", &genjetArr); TBranch *genjetBr = eventTree->GetBranch("GenJet");
  cout << "NEvents = " << eventTree->GetEntries() << endl;

  
  //Double_t weight = 1;

  // null vector for default four vectors
  cmsana::FourVector null(0.0,0.0,0.0,0.0);

  // loop over events
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    printdebug = false;
    if (ientry % 1000 == 0) cout << "Processed Event " << ientry << endl;
    infoBr->GetEntry(ientry);

    NEvents->Fill(0);
    NPUMean->Fill(info->nPUMean);

    //***********************************************************
    // Definition of Pileup Energy density
    //***********************************************************
    Double_t rho = info->RhoKt6PFJets;

    genparticleArr->Clear(); 
    genjetArr->Clear(); 
    genparticleBr->GetEntry(ientry);
    genjetBr->GetEntry(ientry);

    //***********************************************************
    // Find Gen-Level particles
    //***********************************************************
    const cmsana::TGenParticle *genPhoton1 = 0;
    const cmsana::TGenParticle *genPhoton2 = 0;
    const cmsana::TGenParticle *genB1 = 0;
    const cmsana::TGenParticle *genB2 = 0;
    const cmsana::TGenJet *genBJet1 = 0;
    const cmsana::TGenJet *genBJet2 = 0;
    int NBGenJets = 0;
    int NNonBGenJets = 0;
    int NBJets = 0;

    
    //***********************************************************
    // Find B Partons
    //***********************************************************    
    for(Int_t i=0; i<genparticleArr->GetEntriesFast(); i++) {
      const cmsana::TGenParticle *p = (cmsana::TGenParticle*)((*genparticleArr)[i]);

      if (abs(p->pdgid) == 5 && p->status == 3 && ( ( abs(p->motherPdgID) >= 1 && abs(p->motherPdgID) <= 5 ) || abs(p->motherPdgID) == 21) ) {
        if (!genB1) {
          genB1 = p;
        } else {
          if (!genB2 && ( p->pdgid != genB1->pdgid || cmsana::deltaR(genB1->eta,genB1->phi,p->eta,p->phi) > 0.5) ) genB2 = p;
        }
      }
    }
    

    //***********************************************************
    // Find B Gen-Jets
    //***********************************************************    
    for(Int_t i=0; i<genjetArr->GetEntriesFast(); i++) {
      const cmsana::TGenJet *genjet = (cmsana::TGenJet*)((*genjetArr)[i]);  
      if (!(genjet->pt > 30)) continue;

      if (abs(genjet->matchedPdgId) == 5) {
        NBGenJets++;
        if (abs(genjet->eta) < 2.4) {
          if(!genBJet1) genBJet1 = genjet;
          else if (!genBJet2) genBJet2 = genjet;
        }
      } else {
        NNonBGenJets++;
      }
    }

    //***********************************************************
    //Filter for events with 2 b genjets
    //***********************************************************
    if ( NBGenJets < 2 ) {
      continue;
    }



    //***********************************************************
    // Find Gen Photons
    //***********************************************************    
    for(Int_t i=0; i<genparticleArr->GetEntriesFast(); i++) {
      const cmsana::TGenParticle *p = (cmsana::TGenParticle*)((*genparticleArr)[i]);
      if ( SampleType == cmsana::HHToBBGGEventTree::BBGG) {
        if (p->pdgid == 22 && ( p->motherPdgID == 21 || (abs(p->motherPdgID) >= 1 && abs(p->motherPdgID) <= 6) ) ) {
          if (!genPhoton1) {
            genPhoton1 = p;
          } else if (!genPhoton2) {
            genPhoton2 = p;
          }
        }
      }
    }
    
    //***********************************************************
    // Count Jets
    //***********************************************************    
    int NJets = 0;
    int NCentralJets = 0;
    double jetSumPt = 0;
    for(Int_t i=0; i<genjetArr->GetEntriesFast(); i++) {
      const cmsana::TGenJet *genjet = (cmsana::TGenJet*)((*genjetArr)[i]);  
      if (!(genjet->pt > 30)) continue;

      //Si Note: need to implement jet efficiency too
      if (abs(genjet->eta) < 5.0) {
        jetSumPt += genjet->pt;
        NJets++;
        if (abs(genjet->eta) < 2.5) {
          NCentralJets++;
        }
      }

    }



    //sampleType
    cmsana::HHToBBGGEventTree::SampleType stype = cmsana::HHToBBGGEventTree::none;
    if (SampleType == 6) stype = cmsana::HHToBBGGEventTree::BBGG;
    else {
      cout << "Warning: sample type not diphoton jets. Unintended use of this macro. \n";
    }

    outputEventTree->sampletype = stype;
    //outputEventTree->weight = info->eventweight;
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
    if (genB1) {
      genb1v.SetPt(genB1->pt);
      genb1v.SetEta(genB1->eta);
      genb1v.SetPhi(genB1->phi);
      genb1v.SetM(0);
      outputEventTree->genb1 = genb1v;
    }
    if (genB2) {
      genb2v.SetPt(genB2->pt);
      genb2v.SetEta(genB2->eta);
      genb2v.SetPhi(genB2->phi);
      genb2v.SetM(0);   
      outputEventTree->genb2 = genb2v;
    }

    cmsana::FourVectorM genbjet1v;
    cmsana::FourVectorM genbjet2v;
    outputEventTree->genbjet1 = null;    
    outputEventTree->genbjet2 = null;    
    if (genBJet1) {
      genbjet1v.SetPt(genBJet1->pt);
      genbjet1v.SetEta(genBJet1->eta);
      genbjet1v.SetPhi(genBJet1->phi);
      genbjet1v.SetM(genBJet1->mass);
      outputEventTree->genbjet1 = genbjet1v;
    }
    if (genBJet2) {
      genbjet2v.SetPt(genBJet2->pt);
      genbjet2v.SetEta(genBJet2->eta);
      genbjet2v.SetPhi(genBJet2->phi);
      genbjet2v.SetM(genBJet2->mass);   
      outputEventTree->genbjet2 = genbjet2v;    
    }


    //********************************************************
    //Select photons
    //********************************************************
    cmsana::FourVectorM photon1v;
    cmsana::FourVectorM photon2v;
    cmsana::FourVectorM diphotonv;
    double pho1eff = 0;
    double pho2eff = 0;
    outputEventTree->pho1 = null;      //default 4-vector
    outputEventTree->pho2 = null;
    outputEventTree->diphoton = null;    

    if (genPhoton1) {
      photon1v.SetPt(genPhoton1->pt);
      photon1v.SetEta(genPhoton1->eta);
      photon1v.SetPhi(genPhoton1->phi);
      photon1v.SetM(0);
      outputEventTree->pho1 = photon1v;
      pho1eff = photonEffHist->GetBinContent( 
        photonEffHist->GetXaxis()->FindFixBin( fmax( fmin( genPhoton1->pt, 99.9), 0.01) ),
        photonEffHist->GetYaxis()->FindFixBin( fmax( fmin (genPhoton1->eta, 2.99), -2.99))
        );
      
      if (pho1eff > 1) {
        cout << "ERROR: " << genPhoton1->pt << " " << genPhoton1->eta << " | " << photonEffHist->GetXaxis()->FindFixBin( genPhoton1->pt) << " " << photonEffHist->GetYaxis()->FindFixBin( genPhoton1->eta) << " | " << pho1eff << "\n";
      }

    } 

    if (genPhoton2) {
      photon2v.SetPt(genPhoton2->pt);
      photon2v.SetEta(genPhoton2->eta);
      photon2v.SetPhi(genPhoton2->phi);
      photon2v.SetM(0);
      outputEventTree->pho2 = photon2v;
      pho2eff = photonEffHist->GetBinContent( 
        photonEffHist->GetXaxis()->FindFixBin( fmax( fmin( genPhoton2->pt, 99.9), 0.01) ),
        photonEffHist->GetYaxis()->FindFixBin( fmax( fmin (genPhoton2->eta, 2.99), -2.99))
        );
    }
    
    if (genPhoton1 && genPhoton2) {
      diphotonv = photon1v + photon2v;
      outputEventTree->diphoton = diphotonv;
    }
   

    //********************************************************
    //Select b-jets
    //********************************************************
    cmsana::FourVectorM bjet1v;
    cmsana::FourVectorM bjet2v;
    cmsana::FourVectorM dibjetv;
    double bjet1eff = 0;
    double bjet2eff = 0;
    outputEventTree->bjet1 = null;      //default 4-vector
    outputEventTree->bjet2 = null;
    outputEventTree->dibjet = null;    

    if (genBJet1) {
      bjet1v.SetPt(genBJet1->pt);
      bjet1v.SetEta(genBJet1->eta);
      bjet1v.SetPhi(genBJet1->phi);
      bjet1v.SetM(genBJet1->mass);
      bjet1eff = btagEffHist->GetBinContent( 
        btagEffHist->GetXaxis()->FindFixBin( fmax( fmin( genBJet1->pt, 199.9), 0.01) ),
        btagEffHist->GetYaxis()->FindFixBin( fmax( fmin (genBJet1->eta, 2.49), -2.49))
        );
      outputEventTree->bjet1 = bjet1v;
    } else {

    }

    if (genBJet2) {
      bjet2v.SetPt(genBJet2->pt);
      bjet2v.SetEta(genBJet2->eta);
      bjet2v.SetPhi(genBJet2->phi);
      bjet2v.SetM(genBJet2->mass);
      bjet2eff = btagEffHist->GetBinContent( 
        btagEffHist->GetXaxis()->FindFixBin( fmax( fmin( genBJet2->pt, 199.9), 0.01) ),
        btagEffHist->GetYaxis()->FindFixBin( fmax( fmin (genBJet2->eta, 2.49), -2.49))
        );
      outputEventTree->bjet2 = bjet2v;
    }

    if (genBJet1 && genBJet2) {
      dibjetv = bjet1v + bjet2v;
      outputEventTree->dibjet = dibjetv;
    }


    //********************************************************
    //implement efficiency weights
    //********************************************************    
    outputEventTree->weight = pho1eff*pho2eff*bjet1eff*bjet2eff;

    //********************************************************
    //bbgg system
    //********************************************************    
    cmsana::FourVectorM bbggSystemv;
    outputEventTree->bbgg = null;
    if (genBJet1 && genBJet2 && genPhoton1 && genPhoton2) {
      bbggSystemv = (photon1v + photon2v + bjet1v + bjet2v);
      outputEventTree->bbgg = bbggSystemv;
    }



    //********************************************************
    //NJets
    //********************************************************    
    outputEventTree->njets = NJets;
    outputEventTree->ncentraljets = NCentralJets;

    outputEventTree->nlep = 0;

    //********************************************************
    //Some kinematic variables
    //********************************************************
    outputEventTree->DRgg = -1;
    outputEventTree->DRbb = -1;
    outputEventTree->minDRgb = -1;
    if (genBJet1 && genBJet2 && genPhoton1 && genPhoton2) {
      outputEventTree->DRgg = cmsana::deltaR(genPhoton1->eta, genPhoton1->phi, genPhoton2->eta, genPhoton2->phi);
      outputEventTree->DRbb = cmsana::deltaR(genBJet1->eta, genBJet1->phi, genBJet2->eta, genBJet2->phi);
      outputEventTree->minDRgb = fmin(fmin(fmin( cmsana::deltaR(genPhoton1->eta, genPhoton1->phi, genBJet1->eta, genBJet1->phi), 
                                                 cmsana::deltaR(genPhoton1->eta, genPhoton1->phi, genBJet2->eta, genBJet2->phi)),
                                           cmsana::deltaR(genPhoton2->eta, genPhoton2->phi, genBJet1->eta, genBJet1->phi)),
                                      cmsana::deltaR(genPhoton2->eta, genPhoton2->phi, genBJet2->eta, genBJet2->phi));      
    }

    outputEventTree->pfmet = sqrt( info->pfMEx*info->pfMEx + info->pfMEy*info->pfMEy);
    outputEventTree->pfTrackMET = sqrt( info->pfTrackMEx*info->pfTrackMEx + info->pfTrackMEy*info->pfTrackMEy);
    
    outputEventTree->HT = jetSumPt;
    if (genPhoton1) outputEventTree->HT += genPhoton1->pt;
    if (genPhoton2) outputEventTree->HT += genPhoton2->pt;
    

    //***********************************************************
    //Filter for events with 2 bjets, 2photons, and mgg in [110, 140]
    //***********************************************************
    if ( !(genPhoton1 && genPhoton2 && genBJet1 && genBJet2 && genPhoton1->pt > 25 && genPhoton2->pt > 25 && genBJet1->pt > 30 && genBJet2->pt > 30
           && diphotonv.M() > 110 && diphotonv.M() < 140)
      ) continue;


    //********************************************************
    //Fill Output Tree
    //********************************************************
    outputEventTree->tree_->Fill();
    nEvents++;


    //********************************************************
    //Debug
    //********************************************************    
    if (printdebug) {
      cout << "\n\nDebug Event: " << info->runNum << " " << info->lumiSec << " " << info->evtNum << "\n";
      if (genB1) cout << "GenB1: " << genB1->pt << " " << genB1->eta << " " << genB1->phi << "\n";
      if (genB2) cout << "GenB2: " << genB2->pt << " " << genB2->eta << " " << genB2->phi << "\n";
      if (genPhoton1) cout << "GenPho1: " << genPhoton1->pt << " " << genPhoton1->eta << " " << genPhoton1->phi << "\n";
      if (genPhoton2) cout << "GenPho2: " << genPhoton2->pt << " " << genPhoton2->eta << " " << genPhoton2->phi << "\n";
      cout << "NBGenJets = " << NBGenJets << "\n";

      for(Int_t i=0; i<genparticleArr->GetEntriesFast(); i++) {
        const cmsana::TGenParticle *p = (cmsana::TGenParticle*)((*genparticleArr)[i]);
        if ( (abs(p->pdgid) == 5 && p->pt > 0)
            || 
            ( p->pdgid == 22 && p->status==3 && ( p->motherPdgID == 21 || (abs(p->motherPdgID) >= 1 && abs(p->motherPdgID) <= 6) ) )
          ) {
          cout << p->pdgid << " " << p->status << " " << p->pt << " " << p->eta << " " << p->phi 
               << " | " << p->motherPdgID << "\n";
        }
      }

      for(Int_t i=0; i<genjetArr->GetEntriesFast(); i++) {
        const cmsana::TGenJet *genjet = (cmsana::TGenJet*)((*genjetArr)[i]);
        
        bool isStablePhoton = false;
        bool isTruthB = false;
        for(Int_t j=0; j<genparticleArr->GetEntriesFast(); j++) {
          const cmsana::TGenParticle *p = (cmsana::TGenParticle*)((*genparticleArr)[j]);

          if (p->pdgid == 22 && ( p->motherPdgID == 21 || (abs(p->motherPdgID) >= 1 && abs(p->motherPdgID) <= 6) ) 
              && cmsana::deltaR(p->eta,p->phi,genjet->eta,genjet->phi) < 0.5
            ) {
            //genjet matches to a stable photon. don't consider these
            isStablePhoton = true;
            break;
          }
        }

        if (abs(genjet->matchedPdgId) == 5) isTruthB = true;

        cout << "genjet " << i << " : " << genjet->pt << " " << genjet->eta << " " << genjet->phi << " " << genjet->mass << " | ";
        if (isStablePhoton) cout << " matched photon";
        if (isTruthB) cout << " matched b";
        cout << "\n";

      }
      
 
      cout << "All Gen Particles\n";
      for(Int_t i=0; i<genparticleArr->GetEntriesFast(); i++) {
        const cmsana::TGenParticle *p = (cmsana::TGenParticle*)((*genparticleArr)[i]);
        cout << p->pdgid << " " << p->status << " " << p->pt << " " << p->eta << " " << p->phi 
             << " | " << p->motherPdgID << "\n";
      }
      cout << "\n\n";

    }
    //********************************************************
    //End Debug
    //********************************************************



    //********************************************************
    //Clean up Memory
    //********************************************************

  }

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
  cout << "  <> Output saved in " << outputfile << endl;    
  cout << endl;  
      

}

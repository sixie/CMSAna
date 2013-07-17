//================================================================================================
//
// Select HH->bbgammagamma
//
//root -l -b -q  CMSAna/HHToBBGG/HHToBBGGSelection.C+'("root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-02/HHtoBBGG-14tev-START53_V7A/BACONNtuple.dihiggs-bbgg-14tev.1.root","HHToBBGG.root",1)'
//root -l -b -q  CMSAna/HHToBBGG/HHToBBGGSelection.C+'("root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-00/ttHgg-125-START53_V7A/BACONNtuple_1_1_CaY.root","ttH.root",2)'
//root -l -b -q  CMSAna/HHToBBGG/HHToBBGGSelection.C+'("BACONNtuple.root","test.root",0)'
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


  //*****************************************************************************************
  //Setup Jet Energy Corrections
  //*****************************************************************************************
  std::vector<cmsana::JetCorrectorParameters> correctionParameters;
  correctionParameters.push_back(cmsana::JetCorrectorParameters( ( getenv("CMSSW_BASE") + string("/src/CMSAna/JetEnergyCorrections/data/GR_R_52_V9_L1FastJet_AK5PF.txt")).c_str()));
  correctionParameters.push_back(cmsana::JetCorrectorParameters( ( getenv("CMSSW_BASE") + string("/src/CMSAna/JetEnergyCorrections/data/GR_R_52_V9_L2Relative_AK5PF.txt")).c_str()));
  correctionParameters.push_back(cmsana::JetCorrectorParameters( ( getenv("CMSSW_BASE") + string("/src/CMSAna/JetEnergyCorrections/data/GR_R_52_V9_L3Absolute_AK5PF.txt")).c_str()));
  cmsana::FactorizedJetCorrector *JetCorrector = new cmsana::FactorizedJetCorrector(correctionParameters);


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
  //TClonesArray *photonArr = new TClonesArray("cmsana::TPhoton");
  //TClonesArray *muonArr = new TClonesArray("cmsana::TMuon");
  //TClonesArray *electronArr = new TClonesArray("cmsana::TElectron");
  //TClonesArray *jetArr = new TClonesArray("cmsana::TJet");
  //TClonesArray *pfcandidateArr = new TClonesArray("cmsana::TPFCandidate");
  //edit
  TFile *fileZeroWeight = TFile::Open("/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/Efficiency/Trees/BJetMistagRate_type0_nocuts.root");
  TFile *fileFourWeight = TFile::Open("/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/Efficiency/Trees/BJetMistagRate_type4_nocuts.root");
  TFile *fileRealPhoWeight = TFile::Open("/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/PhotonEfficiency_PromptPhoton.root");
  TH2F *typeZeroWeight = (TH2F*)fileZeroWeight->Get("MistagRate_CSVMedium_Pt_Eta");
  TH2F *typeFourWeight = (TH2F*)fileFourWeight->Get("MistagRate_CSVMedium_Pt_Eta");
  TH2F *RealPhoWeight = (TH2F*)fileRealPhoWeight->Get("Efficiency_PtEta");
  TH2F *bjet1Weight = 0;
  TH2F *bjet2Weight = 0;

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
    //edit
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

    //***********************************************************
    // Find Gen-Level particles
    //***********************************************************
    const cmsana::TGenParticle *genPhoton1 = 0;
    const cmsana::TGenParticle *genPhoton2 = 0;
    const cmsana::TGenParticle *genB1 = 0;
    const cmsana::TGenParticle *genB2 = 0;
    const cmsana::TGenJet *genBJet1 = 0;
    const cmsana::TGenJet *genBJet2 = 0;
    for(Int_t i=0; i<genparticleArr->GetEntriesFast(); i++) {
      const cmsana::TGenParticle *p = (cmsana::TGenParticle*)((*genparticleArr)[i]);

      if ( SampleType == cmsana::HHToBBGGEventTree::HHToBBGG) {
        if (abs(p->pdgid) == 5 && p->motherPdgID == 25 && p->status == 2) {
          if (!genB1) {
            genB1 = p;
          } else {
            if (!genB2) genB2 = p;
          }
        }
      }
      if ( SampleType == cmsana::HHToBBGGEventTree::ttHgg 
           || SampleType == cmsana::HHToBBGGEventTree::ttbar
        ) {
        if (abs(p->pdgid) == 5 && abs(p->motherPdgID) == 6 && p->status == 2) {
          if (!genB1) {
            genB1 = p;
          } else {
            if (!genB2) genB2 = p;
          }
        }
      }
      if ( SampleType == cmsana::HHToBBGGEventTree::ZHgg 
        ) {
        if (abs(p->pdgid) == 5 && p->motherPdgID == 23 && p->status == 2) {
          if (!genB1) {
            genB1 = p;
          } else {
            if (!genB2) genB2 = p;
          }
        }
      }

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
    if (SampleType == 0) stype = cmsana::HHToBBGGEventTree::data;
    if (SampleType == 1) stype = cmsana::HHToBBGGEventTree::HHToBBGG;
    if (SampleType == 2) stype = cmsana::HHToBBGGEventTree::ttHgg;
    if (SampleType == 3) stype = cmsana::HHToBBGGEventTree::ZHgg;
    if (SampleType == 4) stype = cmsana::HHToBBGGEventTree::ggHgg;
    if (SampleType == 5) stype = cmsana::HHToBBGGEventTree::ttbar;
    if (SampleType == 6) stype = cmsana::HHToBBGGEventTree::BBGG;
    if (SampleType == 7) stype = cmsana::HHToBBGGEventTree::GGPlusTwoMistag;
    if (SampleType == 8) stype = cmsana::HHToBBGGEventTree::BBPlusTwoFakePhotons;
    if (SampleType == 9) stype = cmsana::HHToBBGGEventTree::CCMistagPlusTwoFakePhotons;
    if (SampleType == 10) stype = cmsana::HHToBBGGEventTree::TwoLightJetsMistagPlusTwoFakePhotons;

    outputEventTree->sampletype = stype;
    outputEventTree->run = info->runNum;
    outputEventTree->lumi = info->lumiSec;
    outputEventTree->event = info->evtNum;
    outputEventTree->npu = info->nPU;
    outputEventTree->rho = info->RhoKt6PFJets;
    outputEventTree->nvtx = info->nGoodPV;

    cmsana::FourVectorM photon1v;
    cmsana::FourVectorM photon2v;
    cmsana::FourVectorM diphotonv;
    outputEventTree->pho1 = null;      //default 4-vector
    outputEventTree->pho2 = null;
    outputEventTree->diphoton = null;    

    if (genPhoton1) {
      photon1v.SetPt(genPhoton1->pt);
      photon1v.SetEta(genPhoton1->eta);
      photon1v.SetPhi(genPhoton1->phi);
      photon1v.SetM(0);
      outputEventTree->pho1 = photon1v;
    } 

    if (genPhoton2) {
      photon2v.SetPt(genPhoton2->pt);
      photon2v.SetEta(genPhoton2->eta);
      photon2v.SetPhi(genPhoton2->phi);
      photon2v.SetM(0);
      outputEventTree->pho2 = photon2v;
    }
    
    if (genPhoton1 && genPhoton2) {
      diphotonv = photon1v + photon2v;
      outputEventTree->diphoton = diphotonv;
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

    //***********************************************************
    // Find B Gen-Jets
    //***********************************************************    
    UInt_t njets = 0;
    UInt_t ncentraljets = 0;
    vector<const cmsana::TGenJet*> goodBJets;
    for(Int_t i=0; i<genjetArr->GetEntriesFast(); i++) {
      const cmsana::TGenJet *genjet = (cmsana::TGenJet*)((*genjetArr)[i]);

      bool isTruthB = false;
      if ( genB1 && cmsana::deltaR(genjet->eta,genjet->phi,genB1->eta,genB1->phi) < 0.5) isTruthB=true;
      if ( genB2 && cmsana::deltaR(genjet->eta,genjet->phi,genB2->eta,genB2->phi) < 0.5) isTruthB=true;

      if (isTruthB) {
        if (!genBJet1) {
          genBJet1 = genjet;
        } else {
          if (!genBJet2) genBJet2 = genjet;
        }      
      }

      bool passCuts = true;
      if ( !(genPhoton1 && genPhoton2) ) passCuts = false;
      //more cuts
      if (!(genjet->pt > 30)) passCuts = false;
      if (!(fabs(genjet->eta) < 2.4)) passCuts = false;

      if (passCuts && genjet->pt > 30) {
        njets++;
        if (fabs(genjet->eta) < 2.5) ncentraljets++;
      }

      if (passCuts) {      
        //save good bjets
        goodBJets.push_back((cmsana::TGenJet*)genjet->Clone()); //new object created
      }
    }

    // switch bjet1 and bjet2 to match gen bjets
    if (genBJet1) {
      if ( genB2 && cmsana::deltaR(genBJet1->eta, genBJet1->phi, genB2->eta, genB2->phi) < 0.5) {
        const cmsana::TGenJet* tmp = genBJet2;
        genBJet2 = genBJet1;
        genBJet1 = tmp;
      }
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

    //cout << " Number of good BJets: " << goodBJets.size() << endl;
    // now that I have all the possible goodBJets, loop through them and choose two pairwise unique goodBJets to
    // promote to bjet1 and bjet2 and store this event into the tree (don't forget to weight the event)

    //********************************************************
    //Select all pairwise goodBJets (i.e. possible bjets)
    //********************************************************
    if (goodBJets.size() > 1) {
      for (Int_t i=0; i<goodBJets.size()-1; i++) {
        for (Int_t j=i+1; j<goodBJets.size(); j++) {
          const cmsana::TGenJet* bjet1 = 0;
          const cmsana::TGenJet* bjet2 = 0;
          bjet1 = goodBJets[i];
          bjet2 = goodBJets[j];
          //cout << bjet1->pt << " | " << bjet2->pt << endl;
          if (bjet2->pt > bjet1->pt) {
            const cmsana::TGenJet* tmp = bjet2;
            bjet2 = bjet1;
            bjet1 = tmp;
          }


          cout << bjet1->matchedPdgId << " | " << bjet2->matchedPdgId << endl;
          if ( abs(bjet1->matchedPdgId) == 5 || abs(bjet2->matchedPdgId) == 5) continue;
          if ( abs(bjet1->matchedPdgId) == 4) bjet1Weight = typeFourWeight;
          else bjet1Weight = typeZeroWeight;
          if ( abs(bjet2->matchedPdgId) == 4) bjet2Weight = typeFourWeight;
          else bjet2Weight = typeZeroWeight;

      
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
          if (bjet1 && bjet2 && genPhoton1 && genPhoton2) {
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
          if (genPhoton1 && genPhoton2 && bjet1 && bjet2 ) {
            outputEventTree->DRgg = cmsana::deltaR(genPhoton1->eta, genPhoton1->phi, genPhoton2->eta, genPhoton2->phi);
            outputEventTree->DRbb = cmsana::deltaR(bjet1->eta, bjet1->phi, bjet2->eta, bjet2->phi);
            outputEventTree->minDRgb = fmin(fmin(fmin( cmsana::deltaR(genPhoton1->eta, genPhoton1->phi, bjet1->eta, bjet1->phi), 
                                                       cmsana::deltaR(genPhoton1->eta, genPhoton1->phi, bjet2->eta, bjet2->phi)),
                                                 cmsana::deltaR(genPhoton2->eta, genPhoton2->phi, bjet1->eta, bjet1->phi)),
                                            cmsana::deltaR(genPhoton2->eta, genPhoton2->phi, bjet2->eta, bjet2->phi));      
          }
      
          outputEventTree->pfmet = sqrt( info->pfMEx*info->pfMEx + info->pfMEy*info->pfMEy);
          outputEventTree->pfTrackMET = sqrt( info->pfTrackMEx*info->pfTrackMEx + info->pfTrackMEy*info->pfTrackMEy);
          
          //Note: currently not computing HT. we may add it later.

          //********************************************************
          //Fill Output Tree
          //********************************************************
          outputEventTree->weight = bjet1Weight->GetBinContent(bjet1Weight->FindBin(fmax(fmin(bjet1->pt,99.9),0.01),fmax(fmin(bjet1->eta,2.49),-2.49))) * bjet2Weight->GetBinContent(bjet2Weight->FindBin(fmax(fmin(bjet2->pt,99.9),0.01),fmax(fmin(bjet2->eta,2.49),-2.49))) * RealPhoWeight->GetBinContent(RealPhoWeight->FindBin(fmax(fmin(genPhoton1->pt,99.9),0.01),fmax(fmin(genPhoton1->eta,2.49),-2.49))) * RealPhoWeight->GetBinContent(RealPhoWeight->FindBin(fmax(fmin(genPhoton2->pt,99.9),0.01),fmax(fmin(genPhoton2->eta,2.49),-2.49)));
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
  cout << " Number of events with two real photons: " << nEventsTwoRealPho << endl;
  cout << endl;
  cout << "  <> Output saved in " << outputfile << endl;    
  cout << endl;  
      

}

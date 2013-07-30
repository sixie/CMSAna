//================================================================================================
//
// Select HH->bbgammagamma
//
//root -l -b -q  CMSAna/HHToBBGG/HHToBBGGSelection.C+'("root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-04/HHtoBBGG-14tev-START53_V7A/BACONNtuple.dihiggs-bbgg-14tev.1.root","HHToBBGG.root",1)'
//root -l -b -q  CMSAna/HHToBBGG/HHToBBGGSelection.C+'("root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-04/ttHgg-125-START53_V7A/BACONNtuple_1_1_LJ8.root","ttH.root",2)'
//root -l -b -q  CMSAna/HHToBBGG/HHToBBGGSelection.C+'("root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-04/ZHgg-125-START53_V7A/BACONNtuple_1_1_pgY.root","ZH.root",3)'
//root -l -b -q  CMSAna/HHToBBGG/HHToBBGGSelection.C+'("BACONNtuple.root","test.root",0)'
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

void HHToBBGGSelection(const string inputfile,          // input file
                       const string outputfile,         // output directory
                       Int_t SampleType = 1
  ) {
  
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

    NEvents->Fill(0);
    NPUMean->Fill(info->nPUMean);

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
//       cout << p->pdgid << " " << p->status << " " << p->pt << " " << p->eta << " " << p->phi 
//            << " | " << p->motherPdgID << "\n";

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

      if ( SampleType == cmsana::HHToBBGGEventTree::HHToBBGG
           || SampleType == cmsana::HHToBBGGEventTree::ttHgg
           || SampleType == cmsana::HHToBBGGEventTree::ZHgg
           || SampleType == cmsana::HHToBBGGEventTree::ggHgg
        ) {
        if (p->pdgid == 22 && p->motherPdgID == 25) {
          if (!genPhoton1) {
            genPhoton1 = p;
          } else {
            if (!genPhoton2) genPhoton2 = p;
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
    outputEventTree->weight = info->eventweight;
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

    //***********************************************************
    // Find B Gen-Jets
    //***********************************************************    
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
    }

    // switch bjet1 and bjet2 to match gen bjets
    if (genBJet1) {
      if ( genB2 && cmsana::deltaR(genBJet1->eta, genBJet1->phi, genB2->eta, genB2->phi) < 0.5) {
        const cmsana::TGenJet* tmp = genBJet2;
        genBJet2 = genBJet1;
        genBJet1 = tmp;
      }
    }

    //*********************************************************************************************************
    //for debug
    //Note: if the 2 b quarks are too close (ie. within dR of 0.6 or so) the genjets can merge and not be reconstructed
    //*********************************************************************************************************
//     if ((!genBJet1 || !genBJet2) && genB1->pt > 25 && genB2->pt > 25 
//         && cmsana::deltaR(genB1->eta,genB1->phi,genB2->eta,genB2->phi) > 0.5
//           ) {
//     if ( genB1 && genBJet1 && genB1->pt < 25 && genBJet1->pt > 25  ) {
//       cout << "\n\nb quarks\n";

//       cout << "b1,b2 : " << genB1->pt << " , " << genB2->pt << " \n";
//       cout << "gen bjet1 : " << genBJet1->pt << "\n";
//       cout << "DRbb = " << cmsana::deltaR(genB1->eta,genB1->phi,genB2->eta,genB2->phi) << "\n";
//       for(Int_t i=0; i<genparticleArr->GetEntriesFast(); i++) {
//         const cmsana::TGenParticle *p = (cmsana::TGenParticle*)((*genparticleArr)[i]);
//         if (abs(p->pdgid) == 5) {
//           cout << p->pdgid << " " << p->status << " " << p->pt << " " << p->eta << " " << p->phi 
//                << " | " << p->motherPdgID << "\n";
//         }
//       }
//       for(Int_t i=0; i<genjetArr->GetEntriesFast(); i++) {
//         const cmsana::TGenJet *genjet = (cmsana::TGenJet*)((*genjetArr)[i]);
        
//         bool isTruthB = false;
//         if ( cmsana::deltaR(genjet->eta,genjet->phi,genB1->eta,genB1->phi) < 0.5) isTruthB=true;
//         if ( cmsana::deltaR(genjet->eta,genjet->phi,genB2->eta,genB2->phi) < 0.5) isTruthB=true;        
//         cout << "genjet " << i << " : " << genjet->pt << " " << genjet->eta << " " << genjet->phi << " " << genjet->mass << " | " << isTruthB << "\n";                
//       }

//       for(Int_t i=0; i<genparticleArr->GetEntriesFast(); i++) {
//         const cmsana::TGenParticle *p = (cmsana::TGenParticle*)((*genparticleArr)[i]);
//         cout << p->pdgid << " " << p->status << " " << p->pt << " " << p->eta << " " << p->phi 
//              << " | " << p->motherPdgID << "\n";
//       }

//     }
    
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

    const cmsana::TPhoton* photon1 = 0;
    const cmsana::TPhoton* photon2 = 0;
    vector<const cmsana::TPhoton*> goodPhotons;
    for(Int_t i=0; i<photonArr->GetEntriesFast(); i++) {

      const cmsana::TPhoton *pho = (cmsana::TPhoton*)((*photonArr)[i]);

      if (!(pho->pt > 25)) continue;
      if (!(fabs(pho->scEta) < 2.5)) continue;
      if (fabs(pho->scEta) > 1.4442 && fabs(pho->scEta) < 1.566) continue;

      if ( !(passPhotonIDSimpleLoose( pho, pfcandidateArr, info->RhoKt6PFJets,
                                       kDataEra_2012_MC, false))) continue;

      //Do tighter electron veto
      bool passEleVeto = true;
      for(Int_t e=0; e<electronArr->GetEntriesFast(); e++) {
        const cmsana::TElectron *tmpele =
          (cmsana::TElectron*)((*electronArr)[e]);
        if (cmsana::deltaR(tmpele->eta,tmpele->phi,pho->eta,pho->phi) < 0.1) {
          passEleVeto = false;
          break;
        }
      }
      if (!passEleVeto) continue;
      

      //save good photons
      goodPhotons.push_back(pho);     
      if (!photon1) {
        photon1 = pho;
      } else if (pho->pt > photon1->pt) {
        photon2 = photon1;
        photon1 = pho;
      } else if (!photon2) { 
        photon2 = pho;
      } else if (pho->pt > photon2->pt) {
        photon2 = pho;
      } else {
        //it's a third photon, and we don't care about it
      }
    }

    // switch photon1 and photon2 to match gen photons
    if (photon1) {
      if ( genPhoton2 && cmsana::deltaR(photon1->eta, photon1->phi, genPhoton2->eta, genPhoton2->phi) < 0.3) {
        const cmsana::TPhoton* tmp = photon2;
        photon2 = photon1;
        photon1 = tmp;
      }
    }
     


  
    cmsana::FourVectorM photon1v;
    cmsana::FourVectorM photon2v;
    cmsana::FourVectorM diphotonv;
    outputEventTree->pho1 = null;      //default 4-vector
    outputEventTree->pho2 = null;
    outputEventTree->diphoton = null;    

    if (photon1) {
      photon1v.SetPt(photon1->pt);
      photon1v.SetEta(photon1->eta);
      photon1v.SetPhi(photon1->phi);
      photon1v.SetM(0);
      outputEventTree->pho1 = photon1v;
    } 

    if (photon2) {
      photon2v.SetPt(photon2->pt);
      photon2v.SetEta(photon2->eta);
      photon2v.SetPhi(photon2->phi);
      photon2v.SetM(0);
      outputEventTree->pho2 = photon2v;
    }
    
    if (photon1 && photon2) {
      diphotonv = photon1v + photon2v;
      outputEventTree->diphoton = diphotonv;
    }
   

    //********************************************************
    //Select b-jets
    //********************************************************
    double jetSumPt = 0;
    UInt_t njets = 0;
    UInt_t ncentraljets = 0;
    const cmsana::TJet* bjet1 = 0;
    const cmsana::TJet* bjet2 = 0;
    vector<const cmsana::TJet*> goodBJets;
    for(Int_t i=0; i<jetArr->GetEntriesFast(); i++) {

      const cmsana::TJet *rawjet = (cmsana::TJet*)((*jetArr)[i]);
            
      //******************************************************
      //apply jet energy corrections
      //******************************************************
      double JEC = JetEnergyCorrectionFactor(rawjet, JetCorrector, rho, false);      
      cmsana::TJet* jet = (cmsana::TJet*)rawjet->Clone(); //new object created
      jet->pt = JEC * rawjet->rawPt;
      jet->mass = JEC * rawjet->mass;

      bool passCuts = true;

      //remove jets matching to selected photons
      if (photon1 && cmsana::deltaR(photon1->eta, photon1->phi, jet->eta, jet->phi) < 0.3) passCuts = false;
      if (photon2 && cmsana::deltaR(photon2->eta, photon2->phi, jet->eta, jet->phi) < 0.3) passCuts = false;

      //jet counting      
      if (passCuts && jet->pt > 30) {
        njets++;
        jetSumPt += jet->pt;
        if (fabs(jet->eta) < 2.5) ncentraljets++;
      }

      //more cuts
      if (!(jet->pt > 30)) passCuts = false;
      if (!(fabs(jet->eta) < 2.4)) passCuts = false;
      if (!(jet->CombinedSecondaryVertexBJetTagsDisc > 0.679)) passCuts = false;

      if (passCuts) {      
        //save good bjets
        goodBJets.push_back((cmsana::TJet*)jet->Clone()); //new object created
        const cmsana::TJet* goodBJet = goodBJets.back();
        if (!bjet1) {
          bjet1 = goodBJet;
        } else if (goodBJet->pt > bjet1->pt) {
          bjet2 = bjet1;
          bjet1 = goodBJet;
        } else if (!bjet2) { 
          bjet2 = goodBJet;
        } else if (goodBJet->pt > bjet2->pt) {
          bjet2 = goodBJet;
        } else {
          //it's a third bjet, and we don't care about it
        }
      }

      delete jet; //memory clean up
    }

    // switch bjet1 and bjet2 to match gen bjets
    if (bjet1) {
      if ( genBJet2 && cmsana::deltaR(bjet1->eta, bjet1->phi, genBJet2->eta, genBJet2->phi) < 0.5) {
        const cmsana::TJet* tmp = bjet2;
        bjet2 = bjet1;
        bjet1 = tmp;
      }
    }


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
    } else {

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


    //********************************************************
    //Count Additional Leptons
    //********************************************************    
    double leptonSumPt = 0;
    Int_t NLeptons = 0;
    for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
      const cmsana::TMuon *mu = (cmsana::TMuon*)((*muonArr)[i]);

      if (!(mu->pt > 10)) continue;
      if (!(fabs(mu->eta) < 2.4)) continue;

      //pass loose veto cuts
      if (!PassMuonIDVeto(mu)) continue;
      if (!(ComputeMuonPFIsoRings(mu,pfcandidateArr, rho,
                                  kDataEra_2012_MC, kPFIso, 0.0, 0.3, false) / mu->pt < 0.2)) continue;

      //make sure the lepton doesn't overlap with bjets or photons
      if (bjet1 && cmsana::deltaR(mu->eta, mu->phi, bjet1->eta, bjet1->phi) < 0.5) continue;
      if (bjet2 && cmsana::deltaR(mu->eta, mu->phi, bjet2->eta, bjet2->phi) < 0.5) continue;
      if (photon1 && cmsana::deltaR(mu->eta, mu->phi, photon1->eta, photon1->phi) < 0.3) continue;
      if (photon2 && cmsana::deltaR(mu->eta, mu->phi, photon2->eta, photon2->phi) < 0.3) continue;

      leptonSumPt += mu->pt;
      NLeptons++;
    }
    for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
      const cmsana::TElectron *ele = (cmsana::TElectron*)((*electronArr)[i]);

      if (!(ele->pt > 20)) continue;
      if (!(fabs(ele->scEta) < 2.5)) continue;

      //pass loose veto cuts
      if (!PassEleSimpleCutsVeto( ele, pfcandidateArr, rho,
                                  kDataEra_2012_MC, false)) continue;

      //make sure the lepton doesn't overlap with bjets or photons
      if (bjet1 && cmsana::deltaR(ele->eta, ele->phi, bjet1->eta, bjet1->phi) < 0.5) continue;
      if (bjet2 && cmsana::deltaR(ele->eta, ele->phi, bjet2->eta, bjet2->phi) < 0.5) continue;
      if (photon1 && cmsana::deltaR(ele->eta, ele->phi, photon1->eta, photon1->phi) < 0.3) continue;
      if (photon2 && cmsana::deltaR(ele->eta, ele->phi, photon2->eta, photon2->phi) < 0.3) continue;

      leptonSumPt += ele->pt;
      NLeptons++;
    }

    outputEventTree->nlep = NLeptons;

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

    outputEventTree->pfmet = sqrt( info->pfMEx*info->pfMEx + info->pfMEy*info->pfMEy);
    outputEventTree->pfTrackMET = sqrt( info->pfTrackMEx*info->pfTrackMEx + info->pfTrackMEy*info->pfTrackMEy);
    
    outputEventTree->HT = jetSumPt + leptonSumPt;
    if (photon1) outputEventTree->HT += photon1->pt;
    if (photon2) outputEventTree->HT += photon2->pt;
    

    //********************************************************
    //Fill Output Tree
    //********************************************************
    outputEventTree->tree_->Fill();
    nEvents++;


    //********************************************************
    //Debug
    //********************************************************
    printdebug = false;

    if (printdebug) {
      cout << "\n\nDebug Event: " << info->runNum << " " << info->lumiSec << " " << info->evtNum << "\n";
      for(Int_t i=0; i<genparticleArr->GetEntriesFast(); i++) {
        const cmsana::TGenParticle *p = (cmsana::TGenParticle*)((*genparticleArr)[i]);
        if (abs(p->pdgid) == 5 || ( p->pdgid == 22 && p->motherPdgID == 25 )) {
          cout << p->pdgid << " " << p->status << " " << p->pt << " " << p->eta << " " << p->phi 
               << " | " << p->motherPdgID << "\n";
        }
      }

      for(Int_t i=0; i<genjetArr->GetEntriesFast(); i++) {
        const cmsana::TGenJet *genjet = (cmsana::TGenJet*)((*genjetArr)[i]);
        
        bool isTruthB = false;
        if ( genB1 && cmsana::deltaR(genjet->eta,genjet->phi,genB1->eta,genB1->phi) < 0.5) isTruthB=true;
        if ( genB2 && cmsana::deltaR(genjet->eta,genjet->phi,genB2->eta,genB2->phi) < 0.5) isTruthB=true;        
        for(Int_t j=0; j<genparticleArr->GetEntriesFast(); j++) {
          const cmsana::TGenParticle *p = (cmsana::TGenParticle*)((*genparticleArr)[j]);
          if (abs(p->pdgid) == 5  && p->status == 3 && p->pt > 1 && cmsana::deltaR(p->eta,p->phi,genjet->eta,genjet->phi) < 0.5) {
            isTruthB = true;
          }
        }
        cout << "genjet " << i << " : " << genjet->pt << " " << genjet->eta << " " << genjet->phi << " " << genjet->mass << " | " << (isTruthB ? "matched b":"not b") << "\n";

      }
      
      for(Int_t i=0; i<photonArr->GetEntriesFast(); i++) {
        const cmsana::TPhoton *pho = (cmsana::TPhoton*)((*photonArr)[i]);
        if (!(pho->pt > 25)) continue;
        if (!(fabs(pho->scEta) < 2.5)) continue;
        if (fabs(pho->scEta) > 1.4442 && fabs(pho->scEta) < 1.566) continue;
        cout << "Photon " << i << " : " << pho->pt << " " << pho->eta << " " << pho->phi << " : " 
             << passPhotonIDSimpleLoose( pho, pfcandidateArr, info->RhoKt6PFJets, kDataEra_2012_MC, false)
             << "\n";
      }
      for(Int_t i=0; i<jetArr->GetEntriesFast(); i++) {
        const cmsana::TJet *rawjet = (cmsana::TJet*)((*jetArr)[i]);
        //******************************************************
        //apply jet energy corrections
        //******************************************************
        double JEC = JetEnergyCorrectionFactor(rawjet, JetCorrector, rho, false);      
        cmsana::TJet* jet = (cmsana::TJet*)rawjet->Clone(); //new object created
        jet->pt = JEC * rawjet->rawPt;
        jet->mass = JEC * rawjet->mass;

        bool passCuts = true;
        //remove jets matching to selected photons
        if (photon1 && cmsana::deltaR(photon1->eta, photon1->phi, jet->eta, jet->phi) < 0.3) passCuts = false;
        if (photon2 && cmsana::deltaR(photon2->eta, photon2->phi, jet->eta, jet->phi) < 0.3) passCuts = false;
        if (!(jet->pt > 30)) passCuts = false;
        if (passCuts) {
          cout << "Jet " << i << " : " << jet->pt << " " << jet->eta << " " << jet->phi << " : "
               << jet->CombinedSecondaryVertexBJetTagsDisc << " " 
               << ( (jet->CombinedSecondaryVertexBJetTagsDisc > 0.679) ? "BTAG" : " Not BTAG") << " "
               << "\n";
        }
        delete jet;
      }
      for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
        const cmsana::TMuon *mu = (cmsana::TMuon*)((*muonArr)[i]);

        if (!(mu->pt > 10)) continue;
        if (!(fabs(mu->eta) < 2.4)) continue;

        bool overlap = false;
        if (bjet1 && cmsana::deltaR(mu->eta, mu->phi, bjet1->eta, bjet1->phi) < 0.5) overlap = true;
        if (bjet2 && cmsana::deltaR(mu->eta, mu->phi, bjet2->eta, bjet2->phi) < 0.5) overlap = true;
        if (photon1 && cmsana::deltaR(mu->eta, mu->phi, photon1->eta, photon1->phi) < 0.3) overlap = true;
        if (photon2 && cmsana::deltaR(mu->eta, mu->phi, photon2->eta, photon2->phi) < 0.3) overlap = true;

        cout << "muon : " << mu->pt << " " << mu->eta << " " << mu->phi << " : " 
             << PassMuonIDVeto(mu) << " " 
             << (ComputeMuonPFIsoRings(mu,pfcandidateArr, rho, kDataEra_2012_MC, kPFIso, 0.0, 0.3, false) / mu->pt < 0.2)
             << " : " << overlap << " "
             << "\n";
      }
      for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
        const cmsana::TElectron *ele = (cmsana::TElectron*)((*electronArr)[i]);
        
        if (!(ele->pt > 20)) continue;
        if (!(fabs(ele->scEta) < 2.5)) continue;

        bool overlap = false;
        if (bjet1 && cmsana::deltaR(ele->eta, ele->phi, bjet1->eta, bjet1->phi) < 0.5) overlap = true;
        if (bjet2 && cmsana::deltaR(ele->eta, ele->phi, bjet2->eta, bjet2->phi) < 0.5) overlap = true;
        if (photon1 && cmsana::deltaR(ele->eta, ele->phi, photon1->eta, photon1->phi) < 0.3) overlap = true;
        if (photon2 && cmsana::deltaR(ele->eta, ele->phi, photon2->eta, photon2->phi) < 0.3) overlap = true;

        cout << "ele : " << ele->pt << " " << ele->eta << " " << ele->phi << " : " 
             << PassEleSimpleCutsVeto( ele, pfcandidateArr, rho, kDataEra_2012_MC, false) << " : "
             << overlap << " "
             << "\n";

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
    for (uint k=0; k<goodBJets.size(); ++k) {
      if (goodBJets[k]) delete goodBJets[k];
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
   
  
  outFile->Write();
  outFile->Close();
  delete outFile;
  
  cout << " Number of events selected: " << nEvents << endl;
  cout << endl;
  cout << "  <> Output saved in " << outputfile << endl;    
  cout << endl;  
      

}

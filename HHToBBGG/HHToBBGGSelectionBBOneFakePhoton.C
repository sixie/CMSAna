//================================================================================================
//
// Select HH->bbgammagamma
//
//root -l -b -q  CMSAna/HHToBBGG/HHToBBGGSelectionBBOneFakePhoton.C+'("root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-02/HHtoBBGG-14tev-START53_V7A/BACONNtuple.dihiggs-bbgg-14tev.1.root","HHToBBGG.root",1)'
//root -l -b -q  CMSAna/HHToBBGG/HHToBBGGSelectionBBOneFakePhoton.C+'("root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-00/ttHgg-125-START53_V7A/BACONNtuple_1_1_CaY.root","ttH.root",2)'
//root -l -b -q  CMSAna/HHToBBGG/HHToBBGGSelectionBBOneFakePhoton.C+'("BACONNtuple.root","test.root",0)'
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

void HHToBBGGSelectionBBOneFakePhoton(const string inputfile,          // input file
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
  TH1F *N2bjets = new TH1F("N2bjets","Number of Events w/ 2 bjets",1,-0.5,0.5); // Count the number of events that have two real bjets
  TH1F *N4genjets = new TH1F("N4genjets", "Number of Events w/ at 2 bjets and >=2 genjets", 1, -0.5, 0.5); // Count the number of events that have two real bjets + at least two other jets
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
  cout << "NEvents = " << eventTree->GetEntries() << endl;

  // Read efficiency file for Fake Photon Rate
  TFile *glu = TFile::Open("/afs/cern.ch/work/v/vlambert/public/releases/CMSSW_5_3_9_patch3/src/PhotonEfficiencies/Efficiencies/PhotonMistagRate_gluon.root");
  TH2F *gluon_efficiencies = (TH2F*)glu->Get("MistagRate_CSVMedium_PtEta");  // access 2D histogram with efficiency weights for gluons
  TFile *qk = TFile::Open("/afs/cern.ch/work/v/vlambert/public/releases/CMSSW_5_3_9_patch3/src/PhotonEfficiencies/Efficiencies/PhotonMistagRate_quark.root");
  TH2F *quark_efficiencies = (TH2F*)qk->Get("MistagRate_CSVMedium_PtEta");  // access 2D histogram with efficiency weights for quarks

  // Efficiencies for Btagging
  TFile *Btagging = TFile::Open("/afs/cern.ch/work/v/vlambert/public/releases/CMSSW_5_3_9_patch3/src/PhotonEfficiencies/Efficiencies/BTaggingEfficiency_BJetEfficiency.root");
  TH2F *BJet_efficiency = (TH2F*)Btagging->Get("Efficiency_PtEta");
  // Efficiencies for Photons
  TFile *Photon = TFile::Open("/afs/cern.ch/work/v/vlambert/public/releases/CMSSW_5_3_9_patch3/src/PhotonEfficiencies/Efficiencies/PhotonEfficiency_PromptPhoton.root");
  TH2F *Photon_efficiency = (TH2F*)Photon->Get("Efficiency_PtEta");

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

    genparticleBr->GetEntry(ientry);
    genjetBr->GetEntry(ientry);

    //***********************************************************
    // Find Gen-Level particles
    //***********************************************************
    const cmsana::TGenParticle *genPhoton1 = 0;
    const cmsana::TGenParticle *genPhoton2 = 0;
    const cmsana::TGenParticle *photon1 = 0;
    const cmsana::TGenParticle *genB1 = 0;
    const cmsana::TGenParticle *genB2 = 0;
    const cmsana::TGenJet* genBJet1 = 0;
    const cmsana::TGenJet* genBJet2 = 0;
    const cmsana::TGenJet* bjet1 = 0;
    const cmsana::TGenJet* bjet2 = 0;
    Float_t FakeRate1 = 0;
    Float_t FakeRate3 = 0;
    Float_t FakeRate4 = 0;
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
      //********************************************************
      //Select Photons
      //********************************************************
      if (abs(p->pdgid) == 22
	  && ( p->motherPdgID == 25 || abs(p->motherPdgID) <= 6 ||
	       (abs(p->motherPdgID) >= 11 && abs(p->motherPdgID) <= 14) ||
	       abs(p->motherPdgID) == 23 || abs(p->motherPdgID) == 24 || 
	       abs(p->motherPdgID) == 21)) {
	if (p->pt < 20) continue;
	if (fabs(p->eta) > 2.5) continue;
	if (!genPhoton1) {
	  genPhoton1 = p;
	  photon1 = p;
	} else {
	  if (!genPhoton2) {
	    if (cmsana::deltaR(p->eta, p->phi, photon1->eta, photon1->phi) < 0.1) continue;
	    genPhoton2 = p;
	  }
	}
      }
    }
    //********************************************************
    //Select b-jets
    //********************************************************
    for(Int_t gen=0; gen<genjetArr->GetEntriesFast(); gen++) {
      const cmsana::TGenJet *genJ = (cmsana::TGenJet*)((*genjetArr)[gen]);
      if (fabs(genJ->matchedPdgId) == 5) { 
	if (genJ->pt < 30) continue;
	if (fabs(genJ->eta) > 2.4) continue;
	if (!genBJet1) {
	  genBJet1 = genJ;
	  bjet1 = genJ;
	} else {
	  if (!genBJet2) {
	    genBJet2 = genJ;
	    bjet2 = genJ;
	  }
	}
      }
    }
    // Discount all events without any photons and with more than one photon
    if (!genPhoton1) continue;
    if (genPhoton2) continue;
    
    // Cut out events without 2 bjets
    if (!(bjet1 && bjet2)) continue;
    
    // Count the number of events with 2 bjets
    N2bjets->Fill(0);

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
    if (SampleType == 11) stype = cmsana::HHToBBGGEventTree::bbjg;
    if (SampleType == 12) stype = cmsana::HHToBBGGEventTree::ccjg;
    if (SampleType == 13) stype = cmsana::HHToBBGGEventTree::jjjg;

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
    
    //Fill the event bjets and genbjets
    cmsana::FourVectorM bjet1v;
    cmsana::FourVectorM bjet2v;
    cmsana::FourVectorM dibjetv;
    outputEventTree->bjet1 = null;      //default 4-vector
    outputEventTree->bjet2 = null;
    outputEventTree->dibjet = null;    
    
    bjet1v.SetPt(bjet1->pt);
    bjet1v.SetEta(bjet1->eta);
    bjet1v.SetPhi(bjet1->phi);
    bjet1v.SetM(bjet1->mass);
    outputEventTree->bjet1 = bjet1v;
    
    bjet2v.SetPt(bjet2->pt);
    bjet2v.SetEta(bjet2->eta);
    bjet2v.SetPhi(bjet2->phi);
    bjet2v.SetM(bjet2->mass);
    outputEventTree->bjet2 = bjet2v;
    
    dibjetv = bjet1v + bjet2v;
    outputEventTree->dibjet = dibjetv;
    
    cmsana::FourVectorM genbjet1v;
    cmsana::FourVectorM genbjet2v;
    outputEventTree->genbjet1 = null;    
    outputEventTree->genbjet2 = null;    
    
    genbjet1v.SetPt(genBJet1->pt);
    genbjet1v.SetEta(genBJet1->eta);
    genbjet1v.SetPhi(genBJet1->phi);
    genbjet1v.SetM(genBJet1->mass);
    outputEventTree->genbjet1 = genbjet1v;
    
    genbjet2v.SetPt(genBJet2->pt);
    genbjet2v.SetEta(genBJet2->eta);
    genbjet2v.SetPhi(genBJet2->phi);
    genbjet2v.SetM(genBJet2->mass);   
    outputEventTree->genbjet2 = genbjet2v; 
    
    //Count other genjets
    Int_t addgenjets = 0;
    for(Int_t c=0; c<genjetArr->GetEntriesFast(); c++) {
      const cmsana::TGenJet *counterjet = (cmsana::TGenJet*)((*genjetArr)[c]);
      if (counterjet == genBJet1 || counterjet == genBJet2) continue;
      addgenjets++;
    }
    if (addgenjets > 1) {
      N4genjets->Fill(0);
    }

    FakeRate3=0;
    FakeRate4=0;
    //Assign efficiencies for bjets
    FakeRate3 = BJet_efficiency->GetBinContent(BJet_efficiency->FindBin(fmax(fmin(bjet1->pt,199.9),30.1),fmax(fmin(bjet1->eta, 2.39),-2.39)));
    FakeRate4 = BJet_efficiency->GetBinContent(BJet_efficiency->FindBin(fmax(fmin(bjet2->pt,199.9),30.1),fmax(fmin(bjet2->eta, 2.39),-2.39)));
       


    //********************************************************
    // Photons
    //********************************************************
    
    // Set the real photon as photon1
    cmsana::FourVectorM photon1v;
    Float_t photonsPt = 0;
    FakeRate1=0;
    //Assign efficiency for photon
    FakeRate1 = Photon_efficiency->GetBinContent(Photon_efficiency->FindBin(fmax(fmin(bjet1->pt,99.9),20.1),fmax(fmin(bjet1->eta, 2.49),-2.49)));

    outputEventTree->pho1 = null;
    photon1v.SetPt(photon1->pt);
    photon1v.SetEta(photon1->eta);
    photon1v.SetPhi(photon1->phi);
    photon1v.SetM(0);
    photonsPt+=photon1->pt;
    outputEventTree->pho1 = photon1v;

    // Define the second photon
    Float_t FakeRate2 = 0;
    const cmsana::TGenJet* photon2 = 0;
    // Store as photon
    cmsana::FourVectorM photon2v;
    cmsana::FourVectorM diphotonv;
    UInt_t njets = 0;
    UInt_t ncentraljets = 0;
    Float_t goodjetsPt = 0;

    // select genjet
    for(Int_t i=0; i<genjetArr->GetEntriesFast(); i++) {
      const cmsana::TGenJet *genjet2 = (cmsana::TGenJet*)((*genjetArr)[i]);
      
      // look only at non-bjets
      if (genjet2 == genBJet1 || genjet2 == genBJet2) continue;
      if (cmsana::deltaR(photon1->eta, photon1->phi, genjet2->eta, genjet2->phi) < 0.5) continue;

      if (genjet2->pt < 20) continue;
      if (fabs(genjet2->eta) > 2.5) continue;
	
	//Assign genjet as fake photons
	outputEventTree->pho2 = null;
	outputEventTree->diphoton = null; 
	FakeRate2 = 0;
       
	//Assign fake rate for photon 2
	if (fabs(genjet2->matchedPdgId) == 21) { // jet corresponds to a gluon
	  FakeRate2 = gluon_efficiencies->GetBinContent(gluon_efficiencies->FindBin(fmax(fmin(genjet2->pt,99.9),20.1),fmax(fmin(genjet2->eta, 2.49),-2.49)));
	}
	else { // jet corresponds to a quark
	  FakeRate2 = quark_efficiencies->GetBinContent(quark_efficiencies->FindBin(fmax(fmin(genjet2->pt,99.9),20.1),fmax(fmin(genjet2->eta, 2.49),-2.49)));
	}
	photon2 = genjet2;	
	photon2v.SetPt(genjet2->pt);
	photon2v.SetEta(genjet2->eta);
	photon2v.SetPhi(genjet2->phi);
	photon2v.SetM(0);
	photonsPt+=genjet2->pt;

	outputEventTree->pho2 = photon2v;
	diphotonv = photon1v + photon2v;
	outputEventTree->diphoton = diphotonv;

	// count additional jets
	njets = 0;
	ncentraljets = 0;
	goodjetsPt = 0;
	
	for(Int_t k=0; k<genjetArr->GetEntriesFast(); k++) {
	  const cmsana::TGenJet *othergenjet = (cmsana::TGenJet*)((*genjetArr)[k]);
	  
	  // look only at the genjets that are not bjets or photons
	  if ( othergenjet == photon2 || othergenjet == genBJet1 ||othergenjet == genBJet2) continue;
	  if (cmsana::deltaR(photon1->eta, photon1->phi, othergenjet->eta, othergenjet->phi) < 0.5) continue;
	  if (!(othergenjet->pt >30)) continue;
	  njets++;
	  goodjetsPt += othergenjet->pt;
	  if (fabs(othergenjet->eta) < 2.5) ncentraljets++;  
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
	
	if (bjet1 && bjet2 && photon1 && photon2 && njets >= 4) printdebug = true;

	//********************************************************
	//Count Additional Leptons
	//********************************************************    
	Int_t NLeptons = 0;
	Float_t ElectronsPt = 0;
	Float_t MuonsPt = 0;
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
	  MuonsPt+= mu->pt;
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
	  ElectronsPt+= ele->pt;
	  NLeptons++;
	}

	outputEventTree->nlep = NLeptons;
	
	//********************************************************
	//Some kinematic variables
	//********************************************************
	outputEventTree->DRgg = -1;
	outputEventTree->DRbb = -1;
	outputEventTree->minDRgb = -1;
	outputEventTree->HT = 0;
	outputEventTree->MET = 0;
	outputEventTree->pfTrackMET = 0;
	if (photon1 && photon2 && bjet1 && bjet2 ) {
	  outputEventTree->DRgg = cmsana::deltaR(photon1->eta, photon1->phi, photon2->eta, photon2->phi);
	  outputEventTree->DRbb = cmsana::deltaR(bjet1->eta, bjet1->phi, bjet2->eta, bjet2->phi);
	  outputEventTree->minDRgb = fmin(fmin(fmin( cmsana::deltaR(photon1->eta, photon1->phi, bjet1->eta, bjet1->phi), 
						     cmsana::deltaR(photon1->eta, photon1->phi, bjet2->eta, bjet2->phi)),
					       cmsana::deltaR(photon2->eta, photon2->phi, bjet1->eta, bjet1->phi)),
					  cmsana::deltaR(photon2->eta, photon2->phi, bjet2->eta, bjet2->phi));
	};     
	outputEventTree->HT = photonsPt + goodjetsPt + ElectronsPt + MuonsPt;
	outputEventTree->MET = sqrt( info->pfMEx*info->pfMEx + info->pfMEy*info->pfMEy);
	outputEventTree->pfTrackMET = sqrt( info->pfTrackMEx*info->pfTrackMEx + info->pfTrackMEy*info->pfTrackMEy);
	
	
	outputEventTree->pfmet = sqrt( info->pfMEx*info->pfMEx + info->pfMEy*info->pfMEy);

	//********************************************************
	//Fill Output Tree
	//********************************************************

	outputEventTree->weight = FakeRate1 * FakeRate2 * FakeRate3 *FakeRate4;
	outputEventTree->tree_->Fill();
	nEvents++;

	// Debug Fake Rate
	if (FakeRate1*FakeRate2*FakeRate3*FakeRate4 ==0) {
	  cout<< "FakeRate 1 = "<< FakeRate1 <<endl;
	  cout<< "FakeRate 2 = "<< FakeRate2 <<endl;
	  cout<< "FakeRate 3 = "<< FakeRate3 <<endl;
	  cout<< "FakeRate 4 = "<< FakeRate4 <<endl;
	}
	//cout<< "weight = " << FakeRate1*FakeRate2*FakeRate3*FakeRate4<<endl;
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



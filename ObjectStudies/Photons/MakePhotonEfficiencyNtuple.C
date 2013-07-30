//================================================================================================
// root -l CMSAna/ObjectStudies/Photons/MakePhotonEfficiencyNtuple.C+'("root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-04/diphjets-START53_V7A/BACONNtuple_48_1_doS.root","PhotonEfficiencyNtuple.test.root",10)'
// root -l CMSAna/ObjectStudies/Photons/MakePhotonEfficiencyNtuple.C+'("root://eoscms//eos/cms/store/user/sixie/BACON/V00-00-04/ttjll-START53_V7A/BACONNtuple_73_1_GAG.root","PhotonEfficiencyNtuple.electrons.root",20)'
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
#include "CMSAna/DataTree/interface/TPhoton.hh"
#include "CMSAna/DataTree/interface/TJet.hh"
#include "CMSAna/DataTree/interface/TPFCandidate.hh"
#include "CMSAna/DataTree/interface/TGenParticle.hh"
#include "CMSAna/DataTree/interface/TGenJet.hh"
#include "CMSAna/DataTree/interface/TElectron.hh"
#include "CMSAna/ObjectStudies/interface/EfficiencyTree.h"
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
    TDirectory *dir = (TDirectory*)inf->FindObjectAny("Events");
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


void MakePhotonEfficiencyNtuple(const string inputFilename, const string outputFilename, int denominatorType = 10)
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
  cmsana::EfficiencyTree *effTree = new cmsana::EfficiencyTree;
  effTree->CreateTree();
  effTree->tree_->SetAutoFlush(0);

  UInt_t NDenominatorsFilled = 0;
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
  TClonesArray *electronArr = new TClonesArray("cmsana::TElectron");

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
  TBranch *electronBr;
 
  
  //*****************************************************************************************
  //Loop over Data Tree
  //*****************************************************************************************
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);  infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Photon", &photonArr); photonBr  = eventTree->GetBranch("Photon");
  eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr); pfcandidateBr = eventTree->GetBranch("PFCandidate");
  eventTree->SetBranchAddress("GenParticle", &genparticleArr); genparticleBr = eventTree->GetBranch("GenParticle");
  eventTree->SetBranchAddress("GenJet", &genjetArr); genjetBr = eventTree->GetBranch("GenJet");
  eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");

  cout << "InputFile " << inputFilename << " --- Total Events : " << eventTree->GetEntries() << endl;
  for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
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
    electronArr->Clear(); 
    photonBr->GetEntry(ientry);
    pfcandidateBr->GetEntry(ientry);
    genparticleBr->GetEntry(ientry);
    genjetBr->GetEntry(ientry);
    electronBr->GetEntry(ientry);


    //********************************************************
    // Pileup Energy Density
    //********************************************************
    Double_t rho = info->RhoKt6PFJets;
    UInt_t EAEra = kDataEra_2012_Data;


    //********************************************************
    // genjets denominator
    //********************************************************
    if ( denominatorType == 0) {
      for(Int_t k=0; k<genjetArr->GetEntries(); k++) {
        const cmsana::TGenJet *genjet = (cmsana::TGenJet*)((*genjetArr)[k]);
               
        //some kinematic cuts to save ntupling time
        if (!(genjet->pt > 20)) continue;
        if (!(fabs(genjet->eta) < 2.5)) continue;
	
	//make sure it doesn't match to a prompt photon
	bool matchPhoton = false;
	for(Int_t i=0; i<genparticleArr->GetEntriesFast(); i++) {
	  const cmsana::TGenParticle *gen =
	    (cmsana::TGenParticle*)((*genparticleArr)[i]);
	  
	  if (abs(gen->pdgid) == 22 
	      && ( gen->motherPdgID == 25 || abs(gen->motherPdgID) <= 6
		   || 
		   (abs(gen->motherPdgID) >= 11 && abs(gen->motherPdgID)
		    <= 16) ||
		   abs(gen->motherPdgID) == 23 || abs(gen->motherPdgID)
		   == 24 || abs(gen->motherPdgID) == 21
		   )
	      ) {
	    if
	      (cmsana::deltaR(gen->eta,gen->phi,genjet->eta,genjet->phi) < 0.1) {
	      matchPhoton = true;
	      break;
	    }
	  } 
	}
	if (matchPhoton) continue;
	
        bool pass = false;
        //pass Photon cuts
        for(Int_t i=0; i<photonArr->GetEntriesFast(); i++) {
          const cmsana::TPhoton *photon = (cmsana::TPhoton*)((*photonArr)[i]);
          //if (!(photon->pt > 25)) continue;
          //if (!(fabs(photon->scEta) < 2.5)) continue;
	  //if (fabs(photon->scEta ) > 1.4442 && fabs(photon->scEta) < 1.566) continue;
          //match in dR?
          double DR = cmsana::deltaR(photon->eta,photon->phi,genjet->eta,genjet->phi);
          if (!(DR < 0.5)) continue;

	  if ( !(passPhotonIDSimpleLoose( photon, pfcandidateArr, info->RhoKt6PFJets,kDataEra_2012_MC, false))) continue;          
          
          //Do tighter electron veto
          bool passEleVeto = true;
          for(Int_t e=0; e<electronArr->GetEntriesFast(); e++) {
            const cmsana::TElectron *tmpele =
              (cmsana::TElectron*)((*electronArr)[e]);
            if (cmsana::deltaR(tmpele->eta,tmpele->phi,photon->eta,photon->phi) < 0.1) {
              passEleVeto = false;
              break;
            }
          }
          if (!passEleVeto) continue;
          
          pass = true;
          break;
        }
        
        effTree->weight = eventweight;
        effTree->mass = 0;
        effTree->pt = genjet->pt;
        effTree->eta = genjet->eta;
        effTree->phi = genjet->phi;
        effTree->rho = info->RhoKt6PFJets;
        effTree->q = 0;
        effTree->npv = info->nGoodPV;
        effTree->npu = info->nPU;
        effTree->run = info->runNum;
        effTree->lumi = info->lumiSec;
        effTree->event = info->evtNum;
        effTree->pass = pass;
	effTree->matchedPdgId = genjet->matchedPdgId;
        //***********************
        //Fill Denominator
        //***********************
        NDenominatorsFilled++;
        effTree->tree_->Fill();

      }//loop over genjet denominators
    } // if denominator is genjet




    //********************************************************
    // gen photon denominator
    //********************************************************
    if ( denominatorType == 10) {
      vector<double> genphotonsEta;
      vector<double> genphotonsPhi;
      genphotonsEta.clear();
      genphotonsPhi.clear();

      for(Int_t i=0; i<genparticleArr->GetEntriesFast(); i++) {
        const cmsana::TGenParticle *gen =
          (cmsana::TGenParticle*)((*genparticleArr)[i]);
        
        //Look for real photons 
        if (!(abs(gen->pdgid) == 22 
            && ( gen->motherPdgID == 25 
                 || abs(gen->motherPdgID) <= 6
                 || (abs(gen->motherPdgID) >= 11 && abs(gen->motherPdgID) <= 16) 
                 || abs(gen->motherPdgID) == 23 
                 || abs(gen->motherPdgID) == 24 
                 || abs(gen->motherPdgID) == 21
              )
              )) {

          continue;
        }


        //Don't use photon that was already used
        bool usedAlready = false;
        for (int q=0;q<genphotonsEta.size();++q) {
          if (cmsana::deltaR(genphotonsEta[q],genphotonsPhi[q],gen->eta,gen->phi) < 0.1) 
            usedAlready = true;
        }
        if (usedAlready) continue;

        genphotonsEta.push_back(gen->eta);
        genphotonsPhi.push_back(gen->phi);


        bool pass = false;
        //pass Photon cuts
        for(Int_t j=0; j<photonArr->GetEntriesFast(); j++) {
          const cmsana::TPhoton *photon = (cmsana::TPhoton*)((*photonArr)[j]);

          //match in dR?
          double DR = cmsana::deltaR(gen->eta,gen->phi,photon->eta,photon->phi);
          if (!(DR < 0.1)) continue;
	  if ( !(passPhotonIDSimpleLoose( photon, pfcandidateArr, info->RhoKt6PFJets,kDataEra_2012_MC, false))) continue;
          
          //Do tighter electron veto
          bool passEleVeto = true;
          for(Int_t e=0; e<electronArr->GetEntriesFast(); e++) {
            const cmsana::TElectron *tmpele =
              (cmsana::TElectron*)((*electronArr)[e]);
            if (cmsana::deltaR(tmpele->eta,tmpele->phi,photon->eta,photon->phi) < 0.1) {
              passEleVeto = false;
              break;
            }
          }
          if (!passEleVeto) continue;
          
          pass = true;
          break;
        }
        
        effTree->weight = eventweight;
        effTree->mass = 0;
        effTree->pt = gen->pt;
        effTree->eta = gen->eta;
        effTree->phi = gen->phi;
        effTree->rho = info->RhoKt6PFJets;
        effTree->q = 0;
        effTree->npv = info->nGoodPV;
        effTree->npu = info->nPU;
        effTree->run = info->runNum;
        effTree->lumi = info->lumiSec;
        effTree->event = info->evtNum;
        effTree->pass = pass;
	effTree->matchedPdgId = 22;
        //***********************
        //Fill Denominator
        //***********************
        NDenominatorsFilled++;
        effTree->tree_->Fill();

      }
    }
  
    //********************************************************
    // electron denominator
    //********************************************************
    if ( denominatorType == 20) {

      vector<double> geneleEta;
      vector<double> genelePhi;
      geneleEta.clear();
      genelePhi.clear();

      for(Int_t k=0; k<genparticleArr->GetEntriesFast(); k++) {
        const cmsana::TGenParticle *gen =
          (cmsana::TGenParticle*)((*genparticleArr)[k]);
        
        //use only electrons that are daughters of W's or Z's
        if ( !(abs(gen->pdgid) == 11 &&  ( abs(gen->motherPdgID) == 24 || abs(gen->motherPdgID) == 23 )
               ) 
          ) {
          continue;
        }
        
        //Don't use photon that was already used
        bool usedAlready = false;
        for (int q=0;q<geneleEta.size();++q) {
          if (cmsana::deltaR(geneleEta[q],genelePhi[q],gen->eta,gen->phi) < 0.1) 
            usedAlready = true;
        }
        if (usedAlready) continue;
        
        geneleEta.push_back(gen->eta);
        genelePhi.push_back(gen->phi);
   
        
        
        bool pass = false;
        //pass Photon cuts
        for(Int_t j=0; j<photonArr->GetEntriesFast(); j++) {
          const cmsana::TPhoton *photon = (cmsana::TPhoton*)((*photonArr)[j]);
          if (!(photon->pt > 25)) continue;
          if (!(fabs(photon->scEta) < 2.5)) continue;
	  if (fabs(photon->scEta ) > 1.4442 && fabs(photon->scEta) < 1.566) continue;

          //match in dR?
          double DR = cmsana::deltaR(gen->eta,gen->phi,photon->eta,photon->phi);
          if (!(DR < 0.1)) continue;
	  if ( !(passPhotonIDSimpleLoose( photon, pfcandidateArr, info->RhoKt6PFJets,kDataEra_2012_MC, false))) continue;

          //Do tighter electron veto
          bool passEleVeto = true;
          for(Int_t e=0; e<electronArr->GetEntriesFast(); e++) {
            const cmsana::TElectron *tmpele =
              (cmsana::TElectron*)((*electronArr)[e]);
            if (cmsana::deltaR(tmpele->eta,tmpele->phi,photon->eta,photon->phi) < 0.1) {
              passEleVeto = false;
              break;
            }
          }
          if (!passEleVeto) continue;

          pass = true;
          break;
        }
        
        effTree->weight = eventweight;
        effTree->mass = 0;
        effTree->pt = gen->pt;
        effTree->eta = gen->eta;
        effTree->phi = gen->phi;
        effTree->rho = info->RhoKt6PFJets;
        effTree->q = 0;
        effTree->npv = info->nGoodPV;
        effTree->npu = info->nPU;
        effTree->run = info->runNum;
        effTree->lumi = info->lumiSec;
        effTree->event = info->evtNum;
        effTree->pass = pass;
	effTree->matchedPdgId = 11;
        //***********************
        //Fill Denominator
        //***********************
        NDenominatorsFilled++;
        effTree->tree_->Fill();

        }
    }
  


  } //loop over events
  
  cout << "Total Denominators: " << NDenominatorsFilled << endl;

  delete info;
  delete photonArr;
  delete pfcandidateArr;
  delete genparticleArr;
  delete genjetArr;

  outputFile->Write();
  outputFile->Close();
  
  delete outputFile;
  if (effTree) delete effTree;
} 



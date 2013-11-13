//================================================================================================
//
// HZZ4l selection macro
//
//  * plots distributions associated with selected events
//  * prints list of selected events from data
//  * outputs ROOT files of events passing selection for each sample, 
//    which can be processed by plotSelect.C
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
// #include "CMSAna/DataTree/interface/CMSAnaDefs.hh"
#include "CMSAna/DataTree/interface/Types.h"
#include "CMSAna/DataTree/interface/TEventInfo.hh"
#include "CMSAna/DataTree/interface/TElectron.hh"
#include "CMSAna/DataTree/interface/TMuon.hh"
#include "CMSAna/DataTree/interface/TGenParticle.hh"
#include "CMSAna/Utils/CommonTools.hh"
#include "CMSAna/HZZ4L/HLL/LeptonResponseMap.hh"
#include "CMSAna/HZZ4L/HLL/LeptonResolutionData_LegacyPaper/ElectronEfficiencyMap.h"
#include "CMSAna/HZZ4L/HLL/LeptonResolutionData_LegacyPaper/MuonEfficiencyMap.h"

// output data structs
#include "CMSAna/HZZ4L/interface/HZZEventTree.h"

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

void DoMakeHZZEventNtupleUsingFastsim(const string inputFilename, const string outputFilename, int splitIntoNJobs, int jobNumber) 
{  

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  Bool_t printDebug = kFALSE;
  TFile *LeptonResponseFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/syncForRun1LegacyPaper/CMSSW_5_3_9_patch3/src/CMSAna/HZZ4L/HLL/LeptonResolutionData_LegacyPaper/PtResolutionModel_combined.LegacyPaper.root","READ");

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  //--------------------------------------------------------------------------------------------------------------
  // output ntuple structure
  //==============================================================================================================  
  TFile *outputFile = new TFile(outputFilename.c_str(), "RECREATE");
  HZZEventTree *hzzEventTree = new HZZEventTree;
  hzzEventTree->CreateTree();
  hzzEventTree->tree_->SetAutoFlush(0);


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  //
  // Access samples and fill histograms
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  cmsana::TEventInfo *info    = new cmsana::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("cmsana::TElectron");
  TClonesArray *muonArr = new TClonesArray("cmsana::TMuon");
  TClonesArray *genparticleArr = new TClonesArray("cmsana::TGenParticle");
  

  printDebug = kFALSE;

  //********************************************************
  // Get Tree
  //********************************************************
  cout << "Reading File " << inputFilename << endl;
  eventTree = getTreeFromFile(inputFilename.c_str(),"Events"); 
  TBranch *infoBr;
  TBranch *electronBr;
  TBranch *muonBr;
  TBranch *genparticleBr;


  //*****************************************************************************************
  //Loop over muon Data Tree
  //*****************************************************************************************
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("GenParticle", &genparticleArr);         genparticleBr = eventTree->GetBranch("GenParticle");

  cout << "InputFile " << inputFilename << " --- Total Events : " << eventTree->GetEntries() << endl;
  int firstEvent = 0;
  int lastEvent = eventTree->GetEntries();
  if ( splitIntoNJobs > 0) {
    int NEventsPerJob = eventTree->GetEntries() / splitIntoNJobs;
    firstEvent = NEventsPerJob*jobNumber;
    lastEvent = NEventsPerJob*(jobNumber+1);
  }
  cout << "JobNumber " << jobNumber << " : Processing Event " << firstEvent << " to "  << lastEvent << "\n";

  for(UInt_t ientry=firstEvent; ientry < lastEvent; ientry++) {       	
    infoBr->GetEntry(ientry);
    if (ientry % 100 == 0) cout << "Event " << ientry << endl;
	
    //********************************************************
    // Load the branches
    //********************************************************
    electronArr->Clear(); 
    muonArr->Clear(); 
    genparticleArr->Clear(); 
    electronBr->GetEntry(ientry);
    muonBr->GetEntry(ientry);
    genparticleBr->GetEntry(ientry);

    printDebug = false;
    Double_t minGenLeptonDR = 9999;

    const cmsana::TGenParticle *genLep1 = 0;
    const cmsana::TGenParticle *genLep2 = 0;
    const cmsana::TGenParticle *genLep3 = 0;
    const cmsana::TGenParticle *genLep4 = 0;

    int NGenLep = 0;
    bool hasTauTauDecay = false;

    for(Int_t k=0; k<genparticleArr->GetEntries(); k++) {
      const cmsana::TGenParticle *gen = (cmsana::TGenParticle*)((*genparticleArr)[k]);

      if ( (abs(gen->pdgid) == 11 || abs(gen->pdgid) == 13) && gen->status == 3 && ( gen->motherPdgID == 23 || gen->motherPdgID == gen->pdgid )  ) {
        if (!genLep1) genLep1 = gen;
        else if (!genLep2) genLep2 = gen;
        else if (!genLep3) genLep3 = gen;
        else if (!genLep4) genLep4 = gen;
        NGenLep++;
      }

      if ( abs(gen->pdgid) == 15 && gen->status == 3 && gen->motherPdgID == 23 ) {
        hasTauTauDecay = true;
        break;
      }
    }
    	
    //********************************************************
    // Skip Z->tautau decay events
    //********************************************************
    if (hasTauTauDecay) continue;


    //Gen Leptons
    if (NGenLep > 4) {
      cout << "****************************" << endl;
      cout << "warning: more than 4 gen leptons\n";
      cout << "****************************" << endl;
      printDebug = true;
    } 
    
    if (!(genLep1 && genLep2 && genLep3 && genLep4)) {
      cout << "****************************" << endl;
      cout << "warning: did not find 4 gen leptons\n";
      cout << "****************************" << endl;
      printDebug = true;
    } else {

      assert(abs(genLep1->pdgid) == abs(genLep2->pdgid));
      assert(abs(genLep3->pdgid) == abs(genLep4->pdgid));
      assert(genLep1->pdgid*genLep2->pdgid < 0);
      assert(genLep3->pdgid*genLep4->pdgid < 0);

      cmsana::FourVectorM tmpgenLep1v;
      cmsana::FourVectorM tmpgenLep2v;
      cmsana::FourVectorM tmpgenLep3v;
      cmsana::FourVectorM tmpgenLep4v;
      tmpgenLep1v.SetPt(genLep1->pt);
      tmpgenLep1v.SetEta(genLep1->eta);
      tmpgenLep1v.SetPhi(genLep1->phi);
      tmpgenLep1v.SetM(genLep1->mass);
      tmpgenLep2v.SetPt(genLep2->pt);
      tmpgenLep2v.SetEta(genLep2->eta);
      tmpgenLep2v.SetPhi(genLep2->phi);
      tmpgenLep2v.SetM(genLep2->mass);
      tmpgenLep3v.SetPt(genLep3->pt);
      tmpgenLep3v.SetEta(genLep3->eta);
      tmpgenLep3v.SetPhi(genLep3->phi);
      tmpgenLep3v.SetM(genLep3->mass);
      tmpgenLep4v.SetPt(genLep4->pt);
      tmpgenLep4v.SetEta(genLep4->eta);
      tmpgenLep4v.SetPhi(genLep4->phi);
      tmpgenLep4v.SetM(genLep4->mass);

      if (abs((tmpgenLep1v+tmpgenLep2v).M() - 91.1876) >  abs((tmpgenLep3v+tmpgenLep4v).M() - 91.1876)) {
        const cmsana::TGenParticle *tmp1;
        const cmsana::TGenParticle *tmp2;
        tmp1 = genLep1;
        tmp2 = genLep2;
        genLep1 = genLep3;
        genLep2 = genLep4;
        genLep3 = tmp1;
        genLep4 = tmp2;
      }

      if (genLep1->pt < genLep2->pt) {
        const cmsana::TGenParticle *tmp = genLep1;
        genLep1 = genLep2;
        genLep2 = tmp;
      }
      if (genLep3->pt < genLep4->pt) {
        const cmsana::TGenParticle *tmp = genLep3;
        genLep3 = genLep4;
        genLep4 = tmp;
      }

    }

    //******************************************************************
    //print out gen leptons
    //******************************************************************
    if (printDebug) {
      for(Int_t k=0; k<genparticleArr->GetEntries(); k++) {
        const cmsana::TGenParticle *gen = (cmsana::TGenParticle*)((*genparticleArr)[k]);
        
        cout << "GenLepton " << k << " : " << gen->pdgid << " " << gen->status << " : " << gen->pt << " " << gen->eta << " " << gen->phi << " " << gen->motherPdgID  << endl;
      }
     
      continue;
    }

    cmsana::FourVectorM genLep1v;
    cmsana::FourVectorM genLep2v;
    cmsana::FourVectorM genLep3v;
    cmsana::FourVectorM genLep4v;
    genLep1v.SetPt(genLep1->pt);
    genLep1v.SetEta(genLep1->eta);
    genLep1v.SetPhi(genLep1->phi);
    genLep1v.SetM(genLep1->mass);
    genLep2v.SetPt(genLep2->pt);
    genLep2v.SetEta(genLep2->eta);
    genLep2v.SetPhi(genLep2->phi);
    genLep2v.SetM(genLep2->mass);
    genLep3v.SetPt(genLep3->pt);
    genLep3v.SetEta(genLep3->eta);
    genLep3v.SetPhi(genLep3->phi);
    genLep3v.SetM(genLep3->mass);
    genLep4v.SetPt(genLep4->pt);
    genLep4v.SetEta(genLep4->eta);
    genLep4v.SetPhi(genLep4->phi);
    genLep4v.SetM(genLep4->mass);
    cmsana::FourVectorM genZ1v = genLep1v+genLep2v;
    cmsana::FourVectorM genZ2v = genLep3v+genLep4v;
    cmsana::FourVectorM genZZv = genLep1v+genLep2v+genLep3v+genLep4v;

    int genchannel = -1;
    if (abs(genLep1->pdgid) == 11 && abs(genLep3->pdgid) == 11) genchannel = HZZEventTree::kFourE;
    if (abs(genLep1->pdgid) == 13 && abs(genLep3->pdgid) == 13) genchannel = HZZEventTree::kFourMu;
    if (abs(genLep1->pdgid) == 11 && abs(genLep3->pdgid) == 13) genchannel = HZZEventTree::kTwoETwoMu;
    if (abs(genLep1->pdgid) == 13 && abs(genLep3->pdgid) == 11) genchannel = HZZEventTree::kTwoMuTwoE;

    
    //********************************************************
    // Debug Printout
    //********************************************************
//     cout << "GenLep1: " << genLep1->pdgid << " " << genLep1->pt << " " << genLep1->eta << " " << genLep1->phi << "\n";
//     cout << "GenLep2: " << genLep2->pdgid << " " << genLep2->pt << " " << genLep2->eta << " " << genLep2->phi << "\n";
//     cout << "GenLep3: " << genLep3->pdgid << " " << genLep3->pt << " " << genLep3->eta << " " << genLep3->phi << "\n";
//     cout << "GenLep4: " << genLep4->pdgid << " " << genLep4->pt << " " << genLep4->eta << " " << genLep4->phi << "\n";
//     cout << "Z1: " << (genLep1v+genLep2v).M() << " " << (genLep1v+genLep2v).Pt() << " " << (genLep1v+genLep2v).Eta() << "\n";
//     cout << "Z2: " << (genLep3v+genLep4v).M() << " " << (genLep3v+genLep4v).Pt() << " " << (genLep3v+genLep4v).Eta() << "\n";
//     cout << "ZZ: " << (genLep1v+genLep2v+genLep3v+genLep4v).M() << " " << (genLep1v+genLep2v+genLep3v+genLep4v).Pt() << " " << (genLep1v+genLep2v+genLep3v+genLep4v).Eta() << "\n";
    
    



    //********************************************************
    // Selection 
    //********************************************************
    double eventweight = 1;

    cmsana::FourVectorM Lep1v;
    cmsana::FourVectorM Lep2v;
    cmsana::FourVectorM Lep3v;
    cmsana::FourVectorM Lep4v;
    Lep1v.SetPt(0);
    Lep1v.SetEta(genLep1->eta);
    Lep1v.SetPhi(genLep1->phi);
    Lep1v.SetM(genLep1->mass);
    Lep2v.SetPt(0);
    Lep2v.SetEta(genLep2->eta);
    Lep2v.SetPhi(genLep2->phi);
    Lep2v.SetM(genLep2->mass);
    Lep3v.SetPt(0);
    Lep3v.SetEta(genLep3->eta);
    Lep3v.SetPhi(genLep3->phi);
    Lep3v.SetM(genLep3->mass);
    Lep4v.SetPt(0);
    Lep4v.SetEta(genLep4->eta);
    Lep4v.SetPhi(genLep4->phi);
    Lep4v.SetM(genLep4->mass);

    //Efficiency
    double Lep1Eff = 1.0;
    double Lep2Eff = 1.0;
    double Lep3Eff = 1.0;
    double Lep4Eff = 1.0;
    if (abs(genLep1->pdgid) == 11) Lep1Eff = GetElectronEfficiencyPtEta(genLep1->pt, fabs(genLep1->eta));
    else Lep1Eff = GetMuonEfficiencyPtEta(genLep1->pt, fabs(genLep1->eta));
    if (abs(genLep2->pdgid) == 11) Lep2Eff = GetElectronEfficiencyPtEta(genLep2->pt, fabs(genLep2->eta));
    else Lep2Eff = GetMuonEfficiencyPtEta(genLep2->pt, fabs(genLep2->eta));
    if (abs(genLep3->pdgid) == 11) Lep3Eff = GetElectronEfficiencyPtEta(genLep3->pt, fabs(genLep3->eta));
    else Lep3Eff = GetMuonEfficiencyPtEta(genLep3->pt, fabs(genLep3->eta));
    if (abs(genLep4->pdgid) == 11) Lep4Eff = GetElectronEfficiencyPtEta(genLep4->pt, fabs(genLep4->eta));
    else Lep4Eff = GetMuonEfficiencyPtEta(genLep4->pt, fabs(genLep4->eta));

    //Momentum Smearing
    Lep1v.SetPt(GenerateLeptonPtFromModel(LeptonResponseFile, genLep1->pdgid, genLep1->pt, fabs(genLep1->eta)));
    Lep2v.SetPt(GenerateLeptonPtFromModel(LeptonResponseFile, genLep2->pdgid, genLep2->pt, fabs(genLep2->eta)));
    Lep3v.SetPt(GenerateLeptonPtFromModel(LeptonResponseFile, genLep3->pdgid, genLep3->pt, fabs(genLep3->eta)));
    Lep4v.SetPt(GenerateLeptonPtFromModel(LeptonResponseFile, genLep4->pdgid, genLep4->pt, fabs(genLep4->eta)));


    cmsana::FourVectorM Z1v = Lep1v+Lep2v;
    cmsana::FourVectorM Z2v = Lep3v+Lep4v;
    cmsana::FourVectorM ZZv = Lep1v+Lep2v+Lep3v+Lep4v;


    //*******************************************
    // Fill HZZ Event
    //*******************************************
    hzzEventTree->weight = Lep1Eff*Lep2Eff*Lep3Eff*Lep4Eff;
    hzzEventTree->run = info->runNum;
    hzzEventTree->lumi = info->lumiSec;
    hzzEventTree->event = info->evtNum;
    hzzEventTree->rho = info->RhoKt6PFJets;
    hzzEventTree->nvtx = info->nGoodPV;
    hzzEventTree->met = sqrt( info->pfMEx*info->pfMEx + info->pfMEy*info->pfMEy);

    hzzEventTree->genl1id = genLep1->pdgid;
    hzzEventTree->genl1pt = genLep1->pt;
    hzzEventTree->genl1eta = genLep1->eta;
    hzzEventTree->genl1phi = genLep1->phi;
    hzzEventTree->genl2id = genLep2->pdgid;
    hzzEventTree->genl2pt = genLep2->pt;
    hzzEventTree->genl2eta = genLep2->eta;
    hzzEventTree->genl2phi = genLep2->phi;
    hzzEventTree->genl3id = genLep3->pdgid;
    hzzEventTree->genl3pt = genLep3->pt;
    hzzEventTree->genl3eta = genLep3->eta;
    hzzEventTree->genl3phi = genLep3->phi;
    hzzEventTree->genl4id = genLep4->pdgid;
    hzzEventTree->genl4pt = genLep4->pt;
    hzzEventTree->genl4eta = genLep4->eta;
    hzzEventTree->genl4phi = genLep4->phi;

    hzzEventTree->genchannel = genchannel;
    hzzEventTree->genz1mass = genZ1v.M(); 
    hzzEventTree->genz1pt = genZ1v.Pt(); 
    hzzEventTree->genz1eta = genZ1v.Eta(); 
    hzzEventTree->genz1rapidity = genZ1v.Rapidity(); 
    hzzEventTree->genz2mass = genZ2v.M(); 
    hzzEventTree->genz2pt = genZ2v.Pt(); 
    hzzEventTree->genz2eta = genZ2v.Eta(); 
    hzzEventTree->genz2rapidity = genZ2v.Rapidity(); 
    hzzEventTree->genzzmass = genZZv.M(); 
    hzzEventTree->genzzpt = genZZv.Pt(); 
    hzzEventTree->genzzeta = genZZv.Eta(); 
    hzzEventTree->genzzrapidity = genZZv.Rapidity(); 

    hzzEventTree->l1id = genLep1->pdgid;
    hzzEventTree->l1pt = Lep1v.Pt();
    hzzEventTree->l1eta = Lep1v.Eta();
    hzzEventTree->l1phi = Lep1v.Phi();
    hzzEventTree->l2id = genLep2->pdgid;
    hzzEventTree->l2pt = Lep2v.Pt();
    hzzEventTree->l2eta = Lep2v.Eta();
    hzzEventTree->l2phi = Lep2v.Phi();
    hzzEventTree->l3id = genLep3->pdgid;
    hzzEventTree->l3pt = Lep3v.Pt();
    hzzEventTree->l3eta = Lep3v.Eta();
    hzzEventTree->l3phi = Lep3v.Phi();
    hzzEventTree->l4id = genLep4->pdgid;
    hzzEventTree->l4pt = Lep4v.Pt();
    hzzEventTree->l4eta = Lep4v.Eta();
    hzzEventTree->l4phi = Lep4v.Phi();

    hzzEventTree->channel = genchannel;
    hzzEventTree->z1mass = Z1v.M(); 
    hzzEventTree->z1pt = Z1v.Pt(); 
    hzzEventTree->z1eta = Z1v.Eta(); 
    hzzEventTree->z1rapidity = Z1v.Rapidity(); 
    hzzEventTree->z2mass = Z2v.M(); 
    hzzEventTree->z2pt = Z2v.Pt(); 
    hzzEventTree->z2eta = Z2v.Eta(); 
    hzzEventTree->z2rapidity = Z2v.Rapidity(); 
    hzzEventTree->zzmass = ZZv.M(); 
    hzzEventTree->zzpt = ZZv.Pt(); 
    hzzEventTree->zzeta = ZZv.Eta(); 
    hzzEventTree->zzrapidity = ZZv.Rapidity(); 


    hzzEventTree->tree_->Fill();

  } //end loop over data     
  


  delete info;
  delete electronArr;
  delete muonArr;
  delete genparticleArr;

  outputFile->cd();
  outputFile->Write();
  outputFile->Close();
  delete outputFile;
  if (hzzEventTree) delete hzzEventTree;
  
} 



//=== MAIN MACRO =================================================================================================

void MakeHZZEventNtupleUsingFastsim(const string inputFilename, const string outputFilename, int splitIntoNJobs, int jobNumber )
{  

  DoMakeHZZEventNtupleUsingFastsim(inputFilename, outputFilename, splitIntoNJobs, jobNumber);

}

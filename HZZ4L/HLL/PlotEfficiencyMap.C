//root -l CMSAna/HZZ4L/HLL/CreateEfficiencyMap.C+'("HZZEfficiencyMap_HZZ125.root","HZZ125",0)'


//================================================================================================
//
// Create Efficiency Map
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
#include <TH2F.h>                   // 1D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

// data structs
#include "CMSAna/HZZ4L/HLL/LeptonResolutionData_LegacyPaper/ElectronEfficiencyMap.h"
#include "CMSAna/HZZ4L/HLL/LeptonResolutionData_LegacyPaper/MuonEfficiencyMap.h"

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

void PlotEfficiencyMap(const string Label = "") 
{  
  gBenchmark->Start("HZZTemplate");
  string label = Label;
  if (Label != "") label = "_" + Label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================

  //********************************************************
  // Create Arrays to store the map
  //********************************************************
   const UInt_t NPtBins = 15; 
   const UInt_t NEtaBins = 16;
   const UInt_t NPhiBins = 12;
   double ptBins[NPtBins+1] = { 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 50};
   double etaBins[NEtaBins+1] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4442, 1.566, 1.8, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6};
   double phiBins[NPhiBins+1] = { -3.2, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5,  2, 2.5, 3.2 };

   TH2F *ElectronEfficiencyMap = new TH2F("ElectronEfficiencyMap", ";p_{T} [GeV/c];#eta;Efficiency", 50, 0, 100, 50, 0, 2.5);
   TH2F *MuonEfficiencyMap = new TH2F("MuonEfficiencyMap", ";p_{T} [GeV/c];|#eta|;Efficiency", 50, 0, 100, 50, 0, 2.5);

   for (int i=0; i<ElectronEfficiencyMap->GetXaxis()->GetNbins()+1; ++i) {
     for (int j=0; j<ElectronEfficiencyMap->GetYaxis()->GetNbins()+1; ++j) {
       double x = ElectronEfficiencyMap->GetXaxis()->GetBinCenter(i);
       double y = ElectronEfficiencyMap->GetYaxis()->GetBinCenter(j);

       double weightModifier = 1.0;
       if (x < 8) weightModifier = 1.25;
       if (x>= 8 && x < 10) weightModifier = 1.03;
       if (x >= 10 &&  x < 16) weightModifier = 1.02;
       if (x >= 16 &&  x < 20) weightModifier = 1.01;
       if (x >= 20 &&  x < 40) weightModifier = 0.975;
       if (x >= 40 ) weightModifier = 0.97;

       ElectronEfficiencyMap->SetBinContent(i,j, weightModifier*GetElectronEfficiencyPtEta(x,y));
       //cout << x << " " << y << " : " << weightModifier*GetElectronEfficiencyPtEta(x,y) << "\n";

     }
   }

   for (int i=0; i<MuonEfficiencyMap->GetXaxis()->GetNbins()+1; ++i) {
     for (int j=0; j<MuonEfficiencyMap->GetYaxis()->GetNbins()+1; ++j) {
       double x = MuonEfficiencyMap->GetXaxis()->GetBinCenter(i);
       double y = MuonEfficiencyMap->GetYaxis()->GetBinCenter(j);

       double weightModifier = 1.0;
       if (x < 6) weightModifier *= 1.15;
       if (x >= 6 && x < 7) weightModifier *= 1.05;
       if (x >= 7 && x < 10) weightModifier *= 1.015;
       if (x >= 10 &&  x < 20) weightModifier *= 1.015;
       if (x >= 20 &&  x < 40) weightModifier *= 1.00;
       if (x >= 40 ) weightModifier *= 0.985;

       MuonEfficiencyMap->SetBinContent(i,j, weightModifier*GetMuonEfficiencyPtEta(x,y));
       //cout << x << " " << y << " : " << weightModifier*GetMuonEfficiencyPtEta(x,y) << "\n";

     }
   }


   TCanvas *cv = new TCanvas("cv","cv", 800, 600);
   cv->SetRightMargin(0.15);
   ElectronEfficiencyMap->SetMinimum(0.5);
   ElectronEfficiencyMap->Draw("colz");
   cv->SaveAs("ElectronEfficiencyMap.gif");


   cv = new TCanvas("cv","cv", 800, 600);
   cv->SetRightMargin(0.15);
   MuonEfficiencyMap->SetMinimum(0.5);
   MuonEfficiencyMap->Draw("colz");
   cv->SaveAs("MuonEfficiencyMap.gif");

} 



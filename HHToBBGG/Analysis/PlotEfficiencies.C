//================================================================================================
//
// Plot Efficiency Histogram
//
//root -l CMSAna/HHToBBGG/Analysis/PlotEfficiencies.C+'("CMSAna/HHToBBGG/data/PhotonEfficiency_PromptPhoton.root","Efficiency_PtEta","PromptPhoton")'
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
#include <TGraphAsymmErrors.h>                

#include "CMSAna/Utils/EfficiencyUtils.hh"
#include "CMSAna/ObjectStudies/interface/EfficiencyTree.h"

#endif



//=== MAIN MACRO ================================================================================================= 

void PlotEfficiencies(const string inputfile, string histname, string xaxislabel, string yaxislabel, string label = "", double min = -999, double max = -999, bool setLogz = false,
                      double xmin=-999, double xmax=-999
  ) {
  
  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;

  TFile *file = new TFile(inputfile.c_str(),"READ");
  assert(file);
  TH2F *efficiencyHist = (TH2F*)file->Get(histname.c_str());
  assert(efficiencyHist);

  //--------------------------------------------------------------------------------------------------------------
  // Draw
  //==============================================================================================================
  TCanvas *cv =0;

  cv = new TCanvas("cv","cv",800,600);
  efficiencyHist->Draw("colz");
  efficiencyHist->SetTitle(";;;Efficiency");
  if (setLogz) cv->SetLogz();
  if (xaxislabel != "") efficiencyHist->GetXaxis()->SetTitle(xaxislabel.c_str());
  if (yaxislabel != "") efficiencyHist->GetYaxis()->SetTitle(yaxislabel.c_str());
  if (min != -999)  efficiencyHist->SetMinimum(min);
  if (max != -999)  efficiencyHist->SetMaximum(max);
  if (xmin != -999 && xmax != -999)  efficiencyHist->GetXaxis()->SetRangeUser(xmin,xmax);
  cv->SetRightMargin(0.15);
  cv->SaveAs(("EfficiencyPtEta" + Label + ".gif").c_str());
  cv->SaveAs(("EfficiencyPtEta" + Label + ".pdf").c_str());


}

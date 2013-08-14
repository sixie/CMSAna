//root -l CMSAna/ObjectStudies/Photons/plotPhotonIDEffPU.C+

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
#include <TF1.h>
#include "CMSAna/Plotting/MitStyleRemix.hh"
#include "CMSAna/Plotting/CPlot.hh"
#include <TLatex.h>
#include "CMSAna/Utils/EfficiencyUtils.hh"
#include "CMSAna/ObjectStudies/interface/EfficiencyTree.h"

#endif

void plotPhotonIDEffPU() { 
  
  TFile *file1 = new TFile("/afs/cern.ch/user/s/sixie/work/public/Phase2Upgrade/HHToBBGG/PlotsForECFANoteSummer2013/PhotonIDnoheEfficiency_HHtoBBGG-diphjets-START53_V7A.root","READ");
  TFile *file2 = new TFile("/afs/cern.ch/user/s/sixie/work/public/Phase2Upgrade/HHToBBGG/PlotsForECFANoteSummer2013/PhotonIDnoheEfficiency_WZHgg-125-Age0DES_DES17_61_V5.root","READ");
  TGraphAsymmErrors *efficiency_eta1_1 = (TGraphAsymmErrors*) file1->Get("PhotonEfficiency_Eta1");
  TGraphAsymmErrors *efficiency_eta2_1 = (TGraphAsymmErrors*) file1->Get("PhotonEfficiency_Eta2");
  TGraphAsymmErrors *efficiency_eta1_2 = (TGraphAsymmErrors*) file2->Get("PhotonEfficiency_Eta1");
  TGraphAsymmErrors *efficiency_eta2_2 = (TGraphAsymmErrors*) file2->Get("PhotonEfficiency_Eta2");
 
  assert(efficiency_eta1_1); 
  assert(efficiency_eta2_1); 
  assert(efficiency_eta1_2); 
  assert(efficiency_eta2_2); 

  efficiency_eta1_1->GetXaxis()->SetRangeUser( 0, 200); 
  efficiency_eta2_1->GetXaxis()->SetRangeUser( 0, 200); 
  efficiency_eta1_2->GetXaxis()->SetRangeUser( 0, 200); 
  efficiency_eta2_2->GetXaxis()->SetRangeUser( 0, 200); 

  TCanvas *cv1 = 0; 
  TCanvas *cv2 = 0; 
 
  TString outdirName = "CMSAna/ECFAPlots"; 
  SetStyle(); 
  
  TString name = ""; 
  TString title = ""; 
  TString format = ""; 
  TString xtitle = ""; 
  TString ytitle = ""; 
  
  cv1 = new TCanvas("cv1","cv1",800,600); 
  name = "PhotonIDnoheEfficiencyTrend_Eta1"; 
  title = "Photon ID Efficiency versus Pileup (#eta < 1.5)"; 
  format = "gif"; 
  xtitle = "number pileup"; 
  ytitle = "Efficiency";

  CPlot *plot1 = new CPlot(name,title,xtitle,ytitle);
  plot1->sOutDir = outdirName; 
  gSystem->mkdir(outdirName,kTRUE);
  plot1->SetXRange( 0, 200);
  plot1->AddGraph(efficiency_eta1_1,"Run 1 CMS Detector","",kBlue,8,1);
  plot1->AddGraph(efficiency_eta1_2,"Phase I CMS Detector","",kGreen,8,1);
  plot1->SetLegend(0.6,0.2,0.95,0.35); 
  plot1->Draw(cv1,kTRUE,format); 
  delete plot1; 

//-------------------------------------------------------------------------------

  cv2 = new TCanvas("cv2","cv2",800,600);
  name = "PhotonIDnoheEfficiencyTrend_Eta2";
  title = "Photon ID Efficiency versus Pileup (1.5 < #eta < 2.5)";
  format = "gif";
  xtitle = "number pileup";  
  ytitle = "Efficiency";

  CPlot *plot2 = new CPlot(name,title,xtitle,ytitle);
  plot2->sOutDir = outdirName;
  gSystem->mkdir(outdirName,kTRUE);
  plot2->SetXRange( 0, 200);
  plot2->AddGraph(efficiency_eta2_1,"Run 1 CMS Detector","",kBlue,8,1);
  plot2->AddGraph(efficiency_eta2_2,"Phase I CMS detector","",kGreen,8,1);
  plot2->SetLegend(0.6,0.2,0.95,0.35);
  plot2->Draw(cv2,kTRUE,format); 
  delete plot2;

}













// root -l CMSAna/ObjectStudies/BTagging/plotBJetMistagRateTrend.C+'()'

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
#include "CMSAna/ObjectStudies/interface/EfficiencyTree.h"

#endif

void plotBJetMistagRateTrend() { 

  TFile *file1 = new TFile("/afs/cern.ch/user/s/sixie/work/public/Phase2Upgrade/HHToBBGG/PlotsForECFANoteSummer2013/MisBTagRate_CSVtight2_charm_qcd_pt30To120-START53_V7A.root","READ");
  TFile *file2 = new TFile("/afs/cern.ch/user/s/sixie/work/public/Phase2Upgrade/HHToBBGG/PlotsForECFANoteSummer2013/MisBTagRate_CSVtight2_charm_phaseIAge0.root","READ");
  TGraphAsymmErrors *efficiency_eta1_1 = (TGraphAsymmErrors*) file1->Get("MisBTagEfficiency_Eta1");
  TGraphAsymmErrors *efficiency_eta2_1 = (TGraphAsymmErrors*) file1->Get("MisBTagEfficiency_Eta2");
  TGraphAsymmErrors *efficiency_eta1_2 = (TGraphAsymmErrors*) file2->Get("MisBTagEfficiency_Eta1");
  TGraphAsymmErrors *efficiency_eta2_2 = (TGraphAsymmErrors*) file2->Get("MisBTagEfficiency_Eta2");
//  TGraphAsymmErrors *efficiency_eta3_1 = (TGraphAsymmErrors*) file1->Get("MisBTagEfficiency_Eta3");
//  TGraphAsymmErrors *efficiency_eta4_1 = (TGraphAsymmErrors*) file1->Get("MisBTagEfficiency_Eta4");
//  TGraphAsymmErrors *efficiency_eta3_2 = (TGraphAsymmErrors*) file2->Get("MisBTagEfficiency_Eta3");
//  TGraphAsymmErrors *efficiency_eta4_2 = (TGraphAsymmErrors*) file2->Get("MisBTagEfficiency_Eta4");


  assert(efficiency_eta1_1);
  assert(efficiency_eta1_2); 
  assert(efficiency_eta2_1); 
  assert(efficiency_eta2_2); 
//  assert(efficiency_eta3_1);
//  assert(efficiency_eta3_2);
//  assert(efficiency_eta4_1);
//  assert(efficiency_eta4_2);

  efficiency_eta1_1->GetXaxis()->SetRangeUser(0,200); 
  efficiency_eta2_1->GetXaxis()->SetRangeUser(0,200); 
  efficiency_eta1_2->GetXaxis()->SetRangeUser(0,200); 
  efficiency_eta2_2->GetXaxis()->SetRangeUser(0,200);
//  efficiency_eta3_1->GetXaxis()->SetRangeUser(0,200);
//  efficiency_eta4_1->GetXaxis()->SetRangeUser(0,200);
//  efficiency_eta3_2->GetXaxis()->SetRangeUser(0,200);
//  efficiency_eta4_2->GetXaxis()->SetRangeUser(0,200);

  TCanvas *cv1 = 0; 
  TCanvas *cv2 = 0; 
  TCanvas *cv3 = 0; 
  TCanvas *cv4 = 0; 

  TString outdirName = "CMSAna/ECFAPlots"; 
  SetStyle(); 

  TString name = ""; 
  TString title = "";
  TString format = "";
  TString xtitle = "";
  TString ytitle = "";

  cv1 = new TCanvas("cv1","cv1",800,600); 
  name = "MisBTagRatePU_CSVtight2_charm_Eta1"; 
  title = "Mistagged Charm Jets versus Pileup (#eta < 1.2)"; 
  format = "gif"; 
  xtitle = "number of pileup"; 
  ytitle = "Rate"; 
  CPlot *plot1 = new CPlot(name,title,xtitle,ytitle); 
  plot1->SetXRange(0,200);
  plot1->sOutDir = outdirName; 
  gSystem->mkdir(outdirName,kTRUE); 
  plot1->AddGraph(efficiency_eta1_1,"Run 1 CMS Detector","",kBlue,8,1); 
  plot1->AddGraph(efficiency_eta1_2,"Phase I CMS Detector","",kGreen,8,1); 
  plot1->SetLegend(0.5,0.8,0.95,0.9); 
  plot1->Draw(cv1,kTRUE,format); 
  delete plot1; 

  cv2 = new TCanvas("cv2","cv2",800,600);
  name = "MisBTagRatePU_CSVtight2_charm_Eta2";
  title = "Mistagged Charm Jets versus Pileup (1.2 < #eta < 2.5)";
  format = "gif";
  xtitle = "number of pileup";
  ytitle = "Rate";
  CPlot *plot2 = new CPlot(name,title,xtitle,ytitle);
  plot2->SetXRange(0,200);
  plot2->sOutDir = outdirName;
  gSystem->mkdir(outdirName,kTRUE);
  plot2->AddGraph(efficiency_eta2_1,"Run 1 CMS Detector","",kBlue,8,1);
  plot2->AddGraph(efficiency_eta2_2,"Phase I CMS Detector","",kGreen,8,1);
  plot2->SetLegend(0.5,0.8,0.95,0.9);
  plot2->Draw(cv2,kTRUE,format);
  delete plot2;

/*
  cv3 = new TCanvas("cv3","cv3",800,600);
  name = "MisBTagRatePU_charm_Eta3";
  title = "MisBTag Rate (charm) versus Pileup, 1.6 < eta < 2.0";
  format = "gif";
  xtitle = "number of pileup";
  ytitle = "Rate";
  CPlot *plot3 = new CPlot(name,title,xtitle,ytitle);
  plot3->SetXRange(0,200);
  plot3->sOutDir = outdirName;
  gSystem->mkdir(outdirName,kTRUE);
  plot3->AddGraph(efficiency_eta3_1,"qcd_pt30To120-START53_V7A","",kBlue,1,1);
  plot3->AddGraph(efficiency_eta3_2,"WZHgg-DYToEE-DYToMM-Age0","",kGreen,1,1);
  plot3->SetLegend(0.5,0.8,0.95,0.9);
  plot3->Draw(cv3,kTRUE,format);
  delete plot3;

  cv4 = new TCanvas("cv4","cv4",800,600);
  name = "MisBTagRatePU_charm_Eta4";
  title = "MisBTag Rate (charm) versus Pileup, 2.0 < eta < 2.5";
  format = "gif";
  xtitle = "number of pileup";
  ytitle = "Rate";
  CPlot *plot4 = new CPlot(name,title,xtitle,ytitle);
  plot4->SetXRange(0,200);
  plot4->sOutDir = outdirName;
  gSystem->mkdir(outdirName,kTRUE);
  plot4->AddGraph(efficiency_eta4_1,"qcd_pt30To120-START53_V7A","",kBlue,1,1);
  plot4->AddGraph(efficiency_eta4_2,"WZHgg-DYToEE-DYToMM-Age0","",kGreen,1,1);
  plot4->SetLegend(0.5,0.8,0.95,0.9);
  plot4->Draw(cv4,kTRUE,format);
  delete plot4;
*/


}
  


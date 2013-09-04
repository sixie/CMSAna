//
// root -l CMSAna/HHToBBGG/fits/FitMassTwoD.C+'("/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_diphotonMass.root", "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_dibjetMass.root", "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_twoDMass.root", 1, 1)'
//
// for plotOption argument: 0 = don't plot, 1 = plot and save
//
// All fits in this macro are binned fits
//
//


#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TStyle.h>                 // class to handle ROOT plotting style
#include <TLatex.h>									// class for latex
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1D.h>                   // 1D histograms
#include <TH2D.h>                   // 2D histograms
#include <TF1.h>										// 1D fitting functions 
#include <TBenchmark.h>             // class to track macro running statistics
#include <TEfficiency.h>            // class to handle efficiency calculations
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#endif

// RooFit headers
#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include/RooRealVar.h"
#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include/RooCategory.h"
#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include/RooDataSet.h"
#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include/RooDataHist.h"
#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include/RooFormulaVar.h"
#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include/RooSimultaneous.h"
#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include/RooAbsPdf.h"
#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include/RooAddPdf.h"
#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include/RooProdPdf.h"
#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include/RooFitResult.h"
#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include/RooExtendPdf.h"
#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include/RooWorkspace.h"
#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include/RooExponential.h"
#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include/RooPolynomial.h"
#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include/RooGaussian.h"
#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include/RooLandau.h"
#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include/RooCBShape.h"
#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include/RooPlot.h"
#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include/RooArgList.h"

#include "TRandom3.h"
#include "TTree.h"
#include "CMSAna/HHToBBGG/interface/HHToBBGGEventTree.hh"

using namespace RooFit;
TCanvas *cv = 0;
TLatex *tex = 0;

void AddModels(RooWorkspace *ws, const string inputfilePho, Int_t plotOption, Int_t constBkg);

//-------------------------------------------------------------
//Main macro for generating data and fitting
//=============================================================  
void PrelimFits_DiphotonMassInBJetWindow(const string inputfilePho = "/afs/cern.ch/user/s/sixie/work/public/Phase2Upgrade/HHToBBGG/MassFit/InputMassHistograms/noPU/HHToBBGG_SignalBkgd_AfterCuts_diphotonMassBJetWin.root", Int_t plotOption = 1, Int_t constBkg = 0) {

  RooWorkspace *ws = new RooWorkspace("MassFitWorkspace");
  AddModels(ws, inputfilePho, plotOption, constBkg);
  ws->writeToFile("CMSAna/HHToBBGG/data/FitWorkspace_DiphotonMassInBJetWindow.root",kTRUE);
}

void AddModels(RooWorkspace *ws, const string inputfilePho, Int_t plotOption, Int_t constbkg) {

  gStyle->SetOptFit(0111);
  
  //Import the sample histograms
  TFile *filePho = new TFile(inputfilePho.c_str(), "READ");
  TH1F *sigPhoMass = (TH1F*)filePho->Get("sigPhoMassBJetWin");
  TH1F *nonresPhoMass = (TH1F*)filePho->Get("nonresPhoMassBJetWin");
  
  //-------------------------------------------------------------
  //Make Variables for Fit
  //=============================================================
    
  //Variables for diPhoton
  RooRealVar massPho("massPho","M_{#gamma#gamma}",100.,150.,"GeV/c^{2}");
  RooRealVar massPhoCut("massPhoCut","M_{#gamma#gamma} (Cut)",115.,130.,"GeV/c^{2}");
  //Signal diPhoton
  RooRealVar sigMeanPho("sigMeanPho", "#bar{M}_{#gamma#gamma}", 125., 120., 130., "GeV/c^{2}");
  RooRealVar sigSigmaPho("sigSigmaPho", "#sigma_{#gamma#gamma}", 1.5, 1., 2., "GeV/c^{2}");
  RooRealVar sigAlphaPho("sigAlphaPho","#alpha_{#gamma#gamma} (cut)", 1., 0.1, 2.);
  RooRealVar sigPowerPho("sigPowerPho", "n_{#gamma#gamma} (power)", 3., 1., 5.);
  //Non-Resonant Background diPhoton
  RooRealVar expRatePho("expRatePho", "#lambda_{#gamma#gamma}", -.027, -.1, 0.);
  
  //Weights for the fit
  RooRealVar nsig("N (Sig)", "# signal events", 18.5, -50,100); //(12.7 + 5.8)
  RooRealVar nnonres("N (NonResBkg)", "# non-resonant background events", 99.7, 0,1000);
  
  //-------------------------------------------------------------
  //Make PDFs for Signal + Background
  //=============================================================
  //PDFs for diPhoton
  //Preliminary Fit for Cut Signal
  RooCBShape sigPDFPhoCut("sigPDFPhoCut", "Signal PDF #gamma#gamma", massPhoCut, sigMeanPho, sigSigmaPho, sigAlphaPho, sigPowerPho);
  RooDataHist sigPhoData("sigPhoData", "sigPhoData", massPhoCut, sigPhoMass);
  RooFitResult *sigPhoCBResult = sigPDFPhoCut.fitTo(sigPhoData, RooFit::Strategy(2));
  //Signal
  RooCBShape sigPDFPho("sigPDFPho", "Signal PDF #gamma#gamma", massPho, sigMeanPho, sigSigmaPho, sigAlphaPho, sigPowerPho);

  //Non-Resonant Background
  RooDataHist nonresPhoData("nonresPhoData", "nonresPhoData", massPho, nonresPhoMass);
  RooExponential nonresPDFPho("nonresPDFPho", "Non-Resonant Background PDF #gamma#gamma", massPho, expRatePho);
  RooFitResult *nonresPhoResult = nonresPDFPho.fitTo(nonresPhoData, RooFit::Strategy(2));
  
  RooAddPdf modelpdf("modelpdf", "Signal+Background PDF", RooArgList(sigPDFPho, nonresPDFPho), RooArgList(nsig, nnonres));
  
  //Variables to set constant
  sigMeanPho.setConstant();
  sigSigmaPho.setConstant();
  sigAlphaPho.setConstant();
  sigPowerPho.setConstant();
  expRatePho.setConstant(); 
  
  //Add all stuff to workspace
  modelpdf.graphVizTree("fullModel.dot");
  ws->import(sigPhoData, Rename("sigPhoData"));
  ws->import(nonresPhoData, Rename("nonresPhoData"));
  ws->import(sigPDFPhoCut);
  
  ws->import(modelpdf);
  
  //-------------------------------------------------------------
  //Plot the fitting functions for signal+background
  //=============================================================
  if (plotOption == 1) {
    RooPlot *framex = 0;
    
    //Plot the Crystal Ball Fit for Signal Photon
    cv = new TCanvas("cv","cv",800,600);
    framex = massPhoCut.frame(Bins(50));
    sigPhoData.plotOn(framex, DrawOption("B"), DataError(RooAbsData::None), FillColor(kRed));
    sigPDFPhoCut.plotOn(framex);
    framex->SetTitle(""); framex->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");  framex->SetYTitle("Number of Events");
    framex->Draw();
    tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.048);
    tex->SetTextFont(42);
    tex->DrawLatex(0.64, 0.84, Form("#bar{M}_{#gamma#gamma} = %.2f", sigMeanPho.getVal()));
    tex->DrawLatex(0.64, 0.79, Form("#sigma_{#gamma#gamma} = %.2f", sigSigmaPho.getVal()));
    //tex->DrawLatex(0.64, 0.74, Form("#alpha_{#gamma#gamma} = %.2f", sigAlpha.getVal()));
    //tex->DrawLatex(0.64, 0.69, Form("n_{#gamma#gamma} = %.2f", sigPower.getVal()));
    tex->Draw();
    cv->Update();
    cv->SaveAs("Plots/AllSignalBkgd/Fits/InitialFits/sigPhoCutCB.gif");
    
   
    //Plot Exponential Fit to Non-Resonant Background Bjet 
    cv = new TCanvas("cv","cv",800,600);
    framex = massPho.frame(Bins(50));
    nonresPhoData.plotOn(framex, DrawOption("B"), DataError(RooAbsData::None), FillColor(kGreen));
    nonresPDFPho.plotOn(framex);
    framex->SetTitle(""); framex->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");  framex->SetYTitle("Number of Events");
    framex->Draw();
    tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.048);
    tex->SetTextFont(42);
    tex->DrawLatex(0.64, 0.84, Form("#lambda_{#gamma#gamma} = %.3f", expRatePho.getVal()));
    tex->Draw();
    cv->Update();
    cv->SaveAs("Plots/AllSignalBkgd/Fits/InitialFits/nonresPhoExpo.gif");
    
    
  }//end of plotting
}


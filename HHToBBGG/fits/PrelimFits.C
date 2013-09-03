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

void AddModels(RooWorkspace *ws, RooWorkspace *wsGenerate, const string inputfilePho, const string inputfileBjet, const string inputfileTwoD, Int_t plotOption, Int_t constBkg);

//-------------------------------------------------------------
//Main macro for generating data and fitting
//=============================================================  
void PrelimFits(const string inputfilePho = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_diphotonMass.root", const string inputfileBjet = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_dibjetMass.root", const string inputfileTwoD = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_twoDMass.root", Int_t plotOption = 1, Int_t constBkg = 0) {

  RooWorkspace *ws = new RooWorkspace("MassFitWorkspace");
  RooWorkspace *wsGenerate = new RooWorkspace("DataGeneratingWorkspace");
  AddModels(ws, wsGenerate, inputfilePho, inputfileBjet, inputfileTwoD, plotOption, constBkg);
  ws->writeToFile("CMSAna/HHToBBGG/data/FitWorkspace.root",kTRUE);
  wsGenerate->writeToFile("CMSAna/HHToBBGG/data/DataGenerateWorkspace.root");
}

void AddModels(RooWorkspace *ws, RooWorkspace *wsGenerate, const string inputfilePho, const string inputfileBjet, const string inputfileTwoD, Int_t plotOption, Int_t constbkg) {

  gStyle->SetOptFit(0111);
  
  //Import the sample histograms
  TFile *filePho = new TFile(inputfilePho.c_str(), "READ");
  TH1F *sigPhoMass = (TH1F*)filePho->Get("sigPhoMass");
  TH1F *resPhoMass = (TH1F*)filePho->Get("resPhoMass");
  TH1F *nonresPhoMass = (TH1F*)filePho->Get("nonresPhoMass");

  TFile *fileBjet = new TFile(inputfileBjet.c_str(), "READ");
  TH1F *sigBjetMass = (TH1F*)fileBjet->Get("sigBjetMass");
  TH1F *resBjetMass = (TH1F*)fileBjet->Get("resBjetMass");
  TH1F *resBjetMassExt = (TH1F*)fileBjet->Get("resBjetMassExt");
  TH1F *nonresBjetMass = (TH1F*)fileBjet->Get("nonresBjetMass");
  
  TFile *fileTwoD = new TFile(inputfileTwoD.c_str(), "READ");
  TH2F *sigMassTwoD = (TH2F*)fileTwoD->Get("sigMassPhoBjet");
  TH2F *resMassTwoD = (TH2F*)fileTwoD->Get("resMassPhoBjet");
  TH2F *nonresMassTwoD = (TH2F*)fileTwoD->Get("nonresMassPhoBjet");
  TH2F *allMassTwoD = (TH2F*)fileTwoD->Get("allMassPhoBjet");
  
  std::cout << "Correlation Factors {sig, res, nonres}: {" << sigMassTwoD->GetCorrelationFactor() << ", " << resMassTwoD->GetCorrelationFactor() << ", " << nonresMassTwoD->GetCorrelationFactor() << "}"<< std::endl;
  
  //-------------------------------------------------------------
  //Make Variables for 2D Fits
  //=============================================================
    
  //Variables for diBjet
  RooRealVar massBjet("massBjet","M_{bb}",70.,200.,"GeV/c^{2}");
  RooRealVar massBjetExt("massBjetExt","M_{bb}",40.,200.,"GeV/c^{2}");
  RooRealVar massBjetCut("massBjetCut","M_{bb}",70.,170.,"GeV/c^{2}");
  //Signal diBjet
  RooRealVar sigMeanBjet("sigMeanBjet", "#bar{M}_{bb}", 125., 119., 131., "GeV/c^{2}");
  RooRealVar sigSigmaBjet("sigSigmaBjet", "#sigma_{bb}", 14., 9., 19., "GeV/c^{2}");
  RooRealVar sigAlphaBjet("sigAlphaBjet","#alpha_{bb} (cut)", 1., 0.1, 2.);
  RooRealVar sigPowerBjet("sigPowerBjet", "n_{bb} (power)", 3., 1., 5.);
  //Resonant Background diBjet
  RooRealVar resMeanBjet("resMeanBjet", "#bar{M}_{bb} (ResBkg)", 90., 80., 100., "GeV/c^{2}");
  RooRealVar resSigmaBjet("resSigmaBjet", "#sigma_{bb} (ResBkg)", 10., 5., 15., "GeV/c^{2}");
  RooRealVar resAlphaBjet("resAlphaBjet","#alpha_{bb} (ResBkg) (cut)", -1.1, -2., -0.1);
  RooRealVar resPowerBjet("resPowerBjet", "n_{bb} (ResBkg) (power)", 3., 1., 5.);
  RooRealVar resExpo("resExpo", "#lambda_{bb} (ResBkg)", -.01, -.1, 0.);
  //Non-Resonant Background diBjet
  RooRealVar expRateBjet("expRateBjet", "#lambda_{bb}", -.018, -.1, 0.);
  RooRealVar p0expBjet("p0expBjet", "a_{bb}", -2.11e-02, -1., 0.);
  RooRealVar p1Bjet("p1Bjet", "b_{bb}", 0., -1., 1.);
  RooRealVar p2Bjet("p2Bjet", "c_{bb}", 0., -1., 1.);
  RooRealVar p3Bjet("p3Bjet", "d_{bb}", 0., -1., 1.);
  RooRealVar p4Bjet("p4Bjet", "e_{bb}", 0., -1., 1.);
  //Weights for Resonant Background
  RooRealVar nbbH("nbbH", "# bbH events", 3.2, 0,5);
  RooRealVar nOthers("nOthers", "# resonant background events - bbH", 24.2, 0,50);
  
  //Variables for diPhoton
  RooRealVar massPho("massPho","M_{#gamma#gamma}",100.,150.,"GeV/c^{2}");
  RooRealVar massPhoCut("massPhoCut","M_{#gamma#gamma} (Cut)",115.,130.,"GeV/c^{2}");
  //Signal diPhoton
  RooRealVar sigMeanPho("sigMeanPho", "#bar{M}_{#gamma#gamma}", 125., 120., 130., "GeV/c^{2}");
  RooRealVar sigSigmaPho("sigSigmaPho", "#sigma_{#gamma#gamma}", 1.5, 1., 2., "GeV/c^{2}");
  RooRealVar sigAlphaPho("sigAlphaPho","#alpha_{#gamma#gamma} (cut)", 1., 0.1, 2.);
  RooRealVar sigPowerPho("sigPowerPho", "n_{#gamma#gamma} (power)", 3., 1., 5.);
  //Resonant Background diPhoton
  RooRealVar resMeanPho("resMeanPho", "#bar{M}_{#gamma#gamma} (ResBkg)", 125., 120., 130., "GeV/c^{2}");
  RooRealVar resSigmaPho("resSigmaPho", "#sigma_{#gamma#gamma} (ResBkg)", 1.5, 1., 2., "GeV/c^{2}");
  RooRealVar resAlphaPho("resAlphaPho","#alpha_{#gamma#gamma} (ResBkg) (cut)", 1., 0.1, 2.);
  RooRealVar resPowerPho("resPowerPho", "n_{#gamma#gamma} (ResBkg) (power)", 3., 1., 5.);
  //Non-Resonant Background diPhoton
  RooRealVar expRatePho("expRatePho", "#lambda_{#gamma#gamma}", -.027, -.1, 0.);
  RooRealVar p0expPho("p0expPho", "a_{#gamma#gamma}", -2.11e-02, -1., 0.);
  RooRealVar p1Pho("p1Pho", "b_{#gamma#gamma}", 0., -1., 1.);
  RooRealVar p2Pho("p2Pho", "c_{#gamma#gamma}", 0., -1., 1.);
  RooRealVar p3Pho("p3Pho", "d_{#gamma#gamma}", 0., -1., 1.);
  RooRealVar p4Pho("p4Pho", "e_{#gamma#gamma}", 0., -1., 1.);
  
  //Weights for the 2D fit
  RooRealVar nsig("N (Sig)", "# signal events", 16.3, -50,100);
  RooRealVar nres("N (ResBkg)", "# resonant background events", 27.4, -50,100);
  RooRealVar nnonres("N (NonResBkg)", "# non-resonant background events", 281.3, 0,1000);
  
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
  //Preliminary Fit for Cut Resonant Background 
  RooCBShape resPDFPhoCut("resPDFPhoCut", "Resonant Background Cut PDF #gamma#gamma", massPhoCut, resMeanPho, resSigmaPho, resAlphaPho, resPowerPho);
  RooDataHist resPhoData("resPhoData", "resPhoData", massPhoCut, resPhoMass);
  RooFitResult *resPhoCBResult = resPDFPhoCut.fitTo(resPhoData, RooFit::Strategy(2));
  //Resonant Background
  RooCBShape resPDFPho("resPDFPho", "Resonant Background PDF #gamma#gamma", massPho, resMeanPho, resSigmaPho, resAlphaPho, resPowerPho);
  //Expo*Poly4 Fit for Non-Resonant Background
  RooExponential nonresPDFexpPho("nonresPDFexpPho", "Non-Resonant Exponential PDF #gamma#gamma", massPho, p0expPho);
  RooPolynomial nonresPDFpolyPho("nonresPDFpolyPho", "Non-Resonant Linear PDF #gamma#gamma", massPho, RooArgList(p1Pho, p2Pho, p3Pho, p4Pho));
  RooProdPdf nonresPDFprodPho("nonresPDFprodPho", "Non-Resonant Background PDF #gamma#gamma", RooArgList(nonresPDFpolyPho,nonresPDFexpPho));
  RooDataHist nonresPhoData("nonresPhoData", "nonresPhoData", massPho, nonresPhoMass);
  RooFitResult *nonresPhoResult = nonresPDFprodPho.fitTo(nonresPhoData, RooFit::Strategy(2));
  //Non-Resonant Background
  RooExponential nonresPDFPho("nonresPDFPho", "Non-Resonant Background PDF #gamma#gamma", massPho, expRatePho);
  nonresPhoResult = nonresPDFPho.fitTo(nonresPhoData, RooFit::Strategy(2));
  
  //PDFs for diBjet
  //Preliminary Fit for Cut Signal
  RooCBShape sigPDFBjetCut("sigPDFBjetCut", "Signal PDF bb", massBjetCut, sigMeanBjet, sigSigmaBjet, sigAlphaBjet, sigPowerBjet);
  RooDataHist sigBjetData("sigBjetData", "sigBjetData", massBjetCut, sigBjetMass);
  RooFitResult *sigBjetCBResult = sigPDFBjetCut.fitTo(sigBjetData, RooFit::Strategy(2));
  //Signal
  RooCBShape sigPDFBjet("sigPDFBjet", "Signal PDF bb", massBjet, sigMeanBjet, sigSigmaBjet, sigAlphaBjet, sigPowerBjet);
  //Preliminary Fit for Extended Resonant Background
  RooCBShape resPDFBjetCBExt("resPDFBjetCBExt", "Resonant CB Background PDF bb", massBjetExt, resMeanBjet, resSigmaBjet, resAlphaBjet, resPowerBjet);
  RooExponential resPDFExpoExt("resPDFExpoExt", "Resonant Exponential Background PDF bb", massBjetExt, resExpo);
  RooAddPdf resPDFBjetExt("resPDFBjetExt", "Extended Resonant Background PDF bb", RooArgList(resPDFBjetCBExt, resPDFExpoExt), RooArgList(nOthers, nbbH));
  RooDataHist resBjetDataExt("resBjetDataExt", "resBjetDataExt", massBjetExt, resBjetMassExt);
  RooFitResult *resBjetCBResult = resPDFBjetExt.fitTo(resBjetDataExt, RooFit::Strategy(2));
  //Resonant Background
  RooCBShape resPDFBjetCB("resPDFBjetCB", "Resonant CB Background PDF bb", massBjet, resMeanBjet, resSigmaBjet, resAlphaBjet, resPowerBjet);
  RooExponential resPDFExpo("resPDFExp", "Resonant Exponential Background PDF bb", massBjet, resExpo);
  RooAddPdf resPDFBjet("resPDFBjet", "Resonant Background PDF bb", RooArgList(resPDFBjetCB, resPDFExpo), RooArgList(nOthers, nbbH));
  //Expo*Poly4 Fit for Non-Resonant Background
  RooExponential nonresPDFexpBjet("nonresPDFexpBjet", "Non-Resonant Exponential PDF bb", massBjet, p0expBjet);
  RooPolynomial nonresPDFpolyBjet("nonresPDFpolyBjet", "Non-Resonant Polynomial PDF bb", massBjet, RooArgList(p1Bjet, p2Bjet, p3Bjet, p4Bjet));
  RooProdPdf nonresPDFprodBjet("nonresPDFprodBjet", "Non-Resonant Background PDF bb", RooArgList(nonresPDFexpBjet,nonresPDFpolyBjet));
  RooDataHist nonresBjetData("nonresBjetData", "nonresBjetData", massBjet, nonresBjetMass);
  RooFitResult *nonresBjetResult = nonresPDFprodBjet.fitTo(nonresBjetData, RooFit::Strategy(2));
  //Non-Resonant Background
  RooExponential nonresPDFBjet("nonresPDFBjet", "Non-Resonant Background PDF bb", massBjet, expRateBjet);
  nonresBjetResult = nonresPDFBjet.fitTo(nonresBjetData, RooFit::Strategy(2));
  
  //PDFs for diPhoton*diBjet (2D)
  RooProdPdf sig2Dpdf("sig2Dpdf", "2D Signal PDF", RooArgList(sigPDFPho, sigPDFBjet));
  RooProdPdf res2Dpdf("res2Dpdf", "2D Resonant Background PDF", RooArgList(resPDFPho, resPDFBjet));
  RooProdPdf nonres2Dpdf("nonres2Dpdf", "2D Non-Resonant Background PDF", RooArgList(nonresPDFPho, nonresPDFBjet));
  RooProdPdf nonres2Dgenerate("nonres2Dgenerate", "2D Non-Resonant Background PDF", RooArgList(nonresPDFprodPho, nonresPDFprodBjet));
  RooAddPdf model2Dpdf("model2Dpdf", "2D Signal+Background PDF", RooArgList(sig2Dpdf, res2Dpdf, nonres2Dpdf), RooArgList(nsig, nres, nnonres));
  RooAddPdf model2Dgenerate("model2Dgenerate", "2D Signal+Backgrond Generator", RooArgList(sig2Dpdf, res2Dpdf, nonres2Dgenerate), RooArgList(nsig, nres, nnonres));
  
  //Variables to set constant
  sigMeanBjet.setConstant();  sigMeanPho.setConstant();
  sigSigmaBjet.setConstant(); sigSigmaPho.setConstant();
  sigAlphaBjet.setConstant(); sigAlphaPho.setConstant();
  sigPowerBjet.setConstant(); sigPowerPho.setConstant();
  resMeanBjet.setConstant();  resMeanPho.setConstant();
  resSigmaBjet.setConstant(); resSigmaPho.setConstant();
  resAlphaBjet.setConstant(); resAlphaPho.setConstant();
  resPowerBjet.setConstant(); resPowerPho.setConstant();
  resExpo.setConstant();
  nbbH.setConstant();
  nOthers.setConstant();
  expRateBjet.setConstant(); p0expBjet.setConstant();
  p1Bjet.setConstant(); p2Bjet.setConstant(); p3Bjet.setConstant(); p4Bjet.setConstant();
  expRatePho.setConstant(); p0expPho.setConstant();
  p1Pho.setConstant(); p2Pho.setConstant(); p3Pho.setConstant(); p4Pho.setConstant();
  
  //Add all 2D stuff to workspace
  model2Dpdf.graphVizTree("fullModel.dot");
  ws->import(sigBjetData, Rename("sigBjetData"));
  ws->import(resBjetDataExt, Rename("resBjetDataExt"));
  wsGenerate->import(nonresBjetData, Rename("nonresBjetData"));
  ws->import(sigPhoData, Rename("sigPhoData"));
  ws->import(resPhoData, Rename("resPhoData"));
  wsGenerate->import(nonresPhoData, Rename("nonresPhoData"));
  ws->import(sigPDFBjetCut);
  ws->import(resPDFBjetExt);
  ws->import(sigPDFPhoCut);
  ws->import(resPDFPhoCut);
  
  ws->import(model2Dpdf);
  wsGenerate->import(model2Dgenerate);
  
  //-------------------------------------------------------------
  //Plot the fitting functions for signal+background
  //=============================================================
  if (plotOption == 1) {
    RooPlot *framex = 0;
    
    //Plot the Crystal Ball Fit for Signal Bjet
    cv = new TCanvas("cv","cv",800,600);
    framex = massBjetCut.frame(Bins(50));
    sigBjetData.plotOn(framex, DrawOption("B"), DataError(RooAbsData::None), FillColor(kRed));
    sigPDFBjetCut.plotOn(framex);
    framex->SetTitle(""); framex->SetXTitle("M_{bb} [GeV/c^{2}]");  framex->SetYTitle("Number of Events");
    framex->Draw();
    tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.048);
    tex->SetTextFont(42);
    tex->DrawLatex(0.64, 0.84, Form("#bar{M}_{bb} = %.2f", sigMeanBjet.getVal()));
    tex->DrawLatex(0.64, 0.79, Form("#sigma_{bb} = %.2f", sigSigmaBjet.getVal()));
    //tex->DrawLatex(0.64, 0.74, Form("#alpha_{bb} = %.2f", sigAlpha.getVal()));
    //tex->DrawLatex(0.64, 0.69, Form("n_{bb} = %.2f", sigPower.getVal()));
    tex->Draw();
    cv->Update();
    cv->SaveAs("Plots/AllSignalBkgd/Fits/InitialFits/sigBjetCutCB.gif");
    
    //Plot the Crystal Ball Fit for Bjet Resonant Background
    cv = new TCanvas("cv","cv",800,600);
    framex = massBjetExt.frame(Bins(50));
    resBjetDataExt.plotOn(framex, DrawOption("B"), DataError(RooAbsData::None), FillColor(kOrange));
    resPDFBjetExt.plotOn(framex);
    framex->SetTitle(""); framex->SetXTitle("M_{bb} [GeV/c^{2}]");  framex->SetYTitle("Number of Events");
    framex->Draw();
    tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.048);
    tex->SetTextFont(42);
    tex->DrawLatex(0.64, 0.84, Form("#bar{M}_{bb} = %.2f", resMeanBjet.getVal()));
    tex->DrawLatex(0.64, 0.79, Form("#sigma_{bb} = %.2f", resSigmaBjet.getVal()));
    //tex->DrawLatex(0.64, 0.74, Form("#alpha_{bb} = %.2f", resAlpha.getVal()));
    //tex->DrawLatex(0.64, 0.69, Form("n_{bb} = %.2f", resPower.getVal()));
    tex->Draw();
    cv->Update();
    cv->SaveAs("Plots/AllSignalBkgd/Fits/InitialFits/resBjetExtCB.gif");
    
    //Plot Exponential Fit to Non-Resonant Background Bjet 
    cv = new TCanvas("cv","cv",800,600);
    framex = massBjet.frame(Bins(50));
    nonresBjetData.plotOn(framex, DrawOption("B"), DataError(RooAbsData::None), FillColor(kGreen));
    nonresPDFBjet.plotOn(framex);
    framex->SetTitle(""); framex->SetXTitle("M_{bb} [GeV/c^{2}]");  framex->SetYTitle("Number of Events");
    framex->Draw();
    tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.048);
    tex->SetTextFont(42);
    tex->DrawLatex(0.64, 0.84, Form("#lambda_{bb} = %.3f", expRateBjet.getVal()));
    tex->Draw();
    cv->Update();
    cv->SaveAs("Plots/AllSignalBkgd/Fits/InitialFits/nonresBjetExpo.gif");
    
    //Plot the Expo*Poly4 Fit of the Non-Resonant Background Bjet
    cv = new TCanvas("cv","cv",800,600);
    framex = massBjet.frame(Bins(50));
    nonresBjetData.plotOn(framex, DrawOption("B"), DataError(RooAbsData::None), FillColor(kGreen));
    nonresPDFprodBjet.plotOn(framex);
    framex->SetTitle(""); framex->SetXTitle("M_{bb} [GeV/c^{2}]");  framex->SetYTitle("Number of Events");
    framex->Draw();
    tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.048);
    tex->SetTextFont(42);
    tex->DrawLatex(0.44, 0.83, Form("e^{#lambda_{bb}x} (1 + %.1ex + %.1ex^{2}", p1Bjet.getVal(), p2Bjet.getVal()));
    tex->DrawLatex(0.52, 0.77, Form(" + %.1ex^{3} + %.1ex^{4})", p3Bjet.getVal(), p4Bjet.getVal()));
    tex->DrawLatex(0.64, 0.71, Form("#lambda_{bb} = %.3f", p0expBjet.getVal()));
    tex->Draw();
    cv->Update();
    cv->SaveAs("Plots/AllSignalBkgd/Fits/InitialFits/nonresBjetExpoPoly.gif");
    
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
    
    //Plot the Crystal Ball Fit for Resonant Photon
    cv = new TCanvas("cv","cv",800,600);
    framex = massPhoCut.frame(Bins(50));
    sigPhoData.plotOn(framex, DrawOption("B"), DataError(RooAbsData::None), FillColor(kOrange));
    sigPDFPhoCut.plotOn(framex);
    framex->SetTitle(""); framex->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");  framex->SetYTitle("Number of Events");
    framex->Draw();
    tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.048);
    tex->SetTextFont(42);
    tex->DrawLatex(0.64, 0.84, Form("#bar{M}_{#gamma#gamma} = %.2f", resMeanPho.getVal()));
    tex->DrawLatex(0.64, 0.79, Form("#sigma_{#gamma#gamma} = %.2f", resSigmaPho.getVal()));
    //tex->DrawLatex(0.64, 0.74, Form("#alpha_{#gamma#gamma} = %.2f", resAlpha.getVal()));
    //tex->DrawLatex(0.64, 0.69, Form("n_{#gamma#gamma} = %.2f", resPower.getVal()));
    tex->Draw();
    cv->Update();
    cv->SaveAs("Plots/AllSignalBkgd/Fits/InitialFits/resPhoCutCB.gif");
    
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
    
    //Plot the Expo*Poly4 Fit of the Non-Resonant Background Photon
    cv = new TCanvas("cv","cv",800,600);
    framex = massPho.frame(Bins(50));
    nonresPhoData.plotOn(framex, DrawOption("B"), DataError(RooAbsData::None), FillColor(kGreen));
    nonresPDFprodPho.plotOn(framex);
    framex->SetTitle(""); framex->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");  framex->SetYTitle("Number of Events");
    framex->Draw();
    tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.048);
    tex->SetTextFont(42);
    tex->DrawLatex(0.44, 0.83, Form("e^{#lambda_{#gamma#gamma}x} (1 + %.1ex + %.1ex^{2}", p1Pho.getVal(), p2Pho.getVal()));
    tex->DrawLatex(0.52, 0.77, Form(" + %.1ex^{3} + %.1ex^{4})", p3Pho.getVal(), p4Pho.getVal()));
    tex->DrawLatex(0.64, 0.71, Form("#lambda_{#gamma#gamma} = %.3f", p0expPho.getVal()));
    tex->Draw();
    cv->Update();
    cv->SaveAs("Plots/AllSignalBkgd/Fits/InitialFits/nonresPhoExpoPoly.gif"); 
    
  }//end of plotting
}


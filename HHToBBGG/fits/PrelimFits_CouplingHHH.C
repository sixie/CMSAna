//
// root -l CMSAna/HHToBBGG/fits/FitMassTwoD.C+'("/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_diphotonMass.root", "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_dibjetMass.root", "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_twoDMass.root", 1, 1)'
//
// for plotOption argument: 0 = don't plot, 1 = plot and save
//
//
// when changing from _tight to regular, remember to change expected yields
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

void AddModels(RooWorkspace *ws, const string inputfilePho, const string inputfileBjet, const string inputfileTwoD, Int_t plotOption, Int_t constBkg);
void plotInitialFit(TH1F *hist, TF1 *func, const string outputName);
void MakeCBPlots(RooWorkspace *ws);

//-------------------------------------------------------------
//Main macro for generating data and fitting
//=============================================================  
void PrelimFits(const string inputfilePho = "/afs/cern.ch/user/s/sixie/work/public/Phase2Upgrade/HHToBBGG/MassFit/InputMassHistograms/noPU/HHToBBGG_SignalBkgd_AfterCuts_diphotonMass.root", const string inputfileBjet = "/afs/cern.ch/user/s/sixie/work/public/Phase2Upgrade/HHToBBGG/MassFit/InputMassHistograms/noPU/HHToBBGG_SignalBkgd_AfterCuts_dibjetMass.root", const string inputfileTwoD = "/afs/cern.ch/user/s/sixie/work/public/Phase2Upgrade/HHToBBGG/MassFit/InputMassHistograms/noPU/HHToBBGG_SignalBkgd_AfterCuts_twoDMass.root", Int_t plotOption = 0, Int_t constBkg = 0) {

  TRandom3 *randomnumber = new TRandom3(1200);
  RooWorkspace *ws = new RooWorkspace("MassFitWorkspace");
  AddModels(ws, inputfilePho, inputfileBjet, inputfileTwoD, plotOption, constBkg);
  if (plotOption == 1) MakeCBPlots(ws);
  ws->writeToFile("FitWorkspace.root",kTRUE);
}

void AddModels(RooWorkspace *ws, const string inputfilePho, const string inputfileBjet, const string inputfileTwoD, Int_t plotOption, Int_t constbkg) {

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
  
  //Preliminary Fits of Signal + Background and Store Parameters
  //fitting diphoton signal and background
  TF1 *sigPhoFunc = new TF1("Signal diPhoton Function","gaus",100,150);
  TF1 *resPhoFunc = new TF1("Resonant Background diPhoton Function","gaus",100,150);
  TF1 *nonresPhoExpo = new TF1("Non-Resonant Background diPhoton Function","expo",100,150);
  sigPhoMass->Fit(sigPhoFunc,"R");
  resPhoMass->Fit(resPhoFunc,"R");
  nonresPhoMass->Fit(nonresPhoExpo,"R");
  
  //fitting dibjet signal and background
  TF1 *sigBjetFunc = new TF1("Signal diBjet Function","gaus",70,200);
  TF1 *resBjetFuncExt = new TF1("Resonant Background diBjet Function","gaus",40,200);
  TF1 *nonresBjetExpo = new TF1("Non-Resonant Background diBjet Function","expo",70,200);
  sigBjetMass->Fit(sigBjetFunc,"R");
  resBjetMass->Fit(resBjetFuncExt,"R");
  nonresBjetMass->Fit(nonresBjetExpo,"R");
  
  //storing parameters
  Double_t parSigPho[3]; Double_t parBkgPho[5];
  sigPhoFunc->GetParameters(parSigPho);
  resPhoFunc->GetParameters(&parBkgPho[0]);
  nonresPhoExpo->GetParameters(&parBkgPho[3]);
  
  Double_t parSigBjet[3]; Double_t parBkgBjet[5];
  sigBjetFunc->GetParameters(&parSigBjet[0]);
  resBjetFuncExt->GetParameters(&parBkgBjet[0]);
  nonresBjetExpo->GetParameters(&parBkgBjet[3]);
  
  //plot preliminary fits
  if (plotOption == 1) {
    plotInitialFit(sigPhoMass,sigPhoFunc,"sigPhoHistFit");
  	plotInitialFit(resPhoMass,resPhoFunc,"resPhoHistFit");
  	plotInitialFit(nonresPhoMass,nonresPhoExpo,"nonresPhoHistFitExpo");
  	
  	plotInitialFit(sigBjetMass,sigBjetFunc,"sigBjetHistFit");
  	plotInitialFit(resBjetMassExt,resBjetFuncExt,"resBjetHistExtFit");
  	plotInitialFit(nonresBjetMass,nonresBjetExpo,"nonresBjetHistFitExpo");
  }
  
  //-------------------------------------------------------------
  //Make Variables for 2D Fits
  //=============================================================
    
  //Variables for diBjet
  RooRealVar massBjet("massBjet","M_{bb}",70.,200.,"GeV/c^{2}");
  RooRealVar massBjetExt("massBjetExt","M_{bb}",40.,200.,"GeV/c^{2}");
  RooRealVar massBjetCut("massBjetCut","M_{bb}",70.,170.,"GeV/c^{2}");
  //Signal diBjet
  RooRealVar sigMeanBjet("sigMeanBjet", "#bar{M}_{bb}", parSigBjet[1], 119., 131., "GeV/c^{2}");
  RooRealVar sigSigmaBjet("sigSigmaBjet", "#sigma_{bb}", parSigBjet[2], 12., 20., "GeV/c^{2}");
  RooRealVar sigAlpha("sigAlpha","#alpha_{bb} (cut)", 1., 0.1, 2.);
  RooRealVar sigPower("sigPower", "n_{bb} (power)", 3., 1., 5.);
  //Resonant Background diBjet
  RooRealVar resMeanBjet("resMeanBjet", "#bar{M}_{bb} (ResBkg)", parBkgBjet[1], 80., 100., "GeV/c^{2}");
  RooRealVar resSigmaBjet("resSigmaBjet", "#sigma_{bb} (ResBkg)", parBkgBjet[2], 5., 15., "GeV/c^{2}");
  RooRealVar resAlpha("resAlpha","#alpha_{bb} (ResBkg) (cut)", -1.1, -2., -0.1);
  RooRealVar resPower("resPower", "n_{bb} (ResBkg) (power)", 3., 1., 5.);
  RooRealVar resExpo("resExpo", "#lambda_{bb} (ResBkg)", -0.01, -0.1, 0.0);
  //Non-Resonant Background diBjet
  RooRealVar expRateBjet("expRateBjet", "#lambda_{bb}", parBkgBjet[4], -1., 0.);
  if (constbkg == 1) expRateBjet.setConstant();
  //Weights for Resonant Background
  RooRealVar nbbH("nbbH", "# bbH events", 3.2, 0,5);
  RooRealVar nOthers("nOthers", "# resonant background events - bbH", 24.2, 0,50);
  
  //Variables for diPhoton (everything fixed except for expo rate)
  RooRealVar massPho("massPho","M_{#gamma#gamma}",100.,150.,"GeV/c^{2}");
  RooRealVar meanPho("meanPho", "#bar{M}_{#gamma#gamma}", parSigPho[1], "GeV/c^{2}");
  RooRealVar sigmaPho("sigmaPho", "#sigma_{#gamma#gamma}", parSigPho[2], "GeV/c^{2}");
  RooRealVar expRatePho("expRatePho", "#lambda_{#gamma#gamma}", parBkgPho[4], -1., 0.);
  if (constbkg == 1) expRatePho.setConstant();
  
  //Weights for the 2D fit
  RooRealVar HHHCoupling("HHHCoupling", "#lambda_{HHH}", 1, -5,5);
  RooFormulaVar nsig("nsig","16.3 * (68.5935 - 44.6022*HHHCoupling+9.13082*HHHCoupling*HHHCoupling) / 33.12212 ", RooArgList(HHHCoupling));

  //RooRealVar nsig("N (Sig)", "# signal events", 16.3, -50,100);

  RooRealVar nres("N (ResBkg)", "# resonant background events", 27.4, -50,100);
  RooRealVar nnonres("N (NonResBkg)", "# non-resonant background events", 281.3, 0,1000);
  
  //-------------------------------------------------------------
  //Make PDFs for Signal + Background
  //=============================================================
  //PDFs for diPhoton
  RooGaussian sigPDFPho("sigPDFPho", "Signal PDF #gamma#gamma", massPho, meanPho, sigmaPho);
  RooGaussian resPDFPho("resPDFPho", "Resonant Background PDF #gamma#gamma", massPho, meanPho, sigmaPho);
  RooExponential nonresPDFPho("nonresPDFPho", "Non-Resonant Background PDF #gamma#gamma", massPho, expRatePho);
  
  //PDFs for diBjet
  //Signal
  RooCBShape sigPDFBjetCut("sigPDFBjetCut", "Signal PDF bb", massBjetCut, sigMeanBjet, sigSigmaBjet, sigAlpha, sigPower);
  //Preliminary Fit for Cut Signal
  RooDataHist sigBjetData("sigBjetData", "sigBjetData", massBjetCut, sigBjetMass);
  RooFitResult *sigCBResult = sigPDFBjetCut.fitTo(sigBjetData, RooFit::Extended(kTRUE), RooFit::Strategy(2));
  RooCBShape sigPDFBjet("sigPDFBjet", "Signal PDF bb", massBjet, sigMeanBjet, sigSigmaBjet, sigAlpha, sigPower);
  //Extended Resonant Background
  RooCBShape resPDFBjetCBExt("resPDFBjetCBExt", "Resonant CB Background PDF bb", massBjetExt, resMeanBjet, resSigmaBjet, resAlpha, resPower);
  RooExponential resPDFExpoExt("resPDFExpoExt", "Resonant Exponential Background PDF bb", massBjetExt, resExpo);
  RooAddPdf resPDFBjetExt("resPDFBjetExt", "Extended Resonant Background PDF bb", RooArgList(resPDFBjetCBExt, resPDFExpoExt), RooArgList(nOthers, nbbH));
  //Preliminary Fit for Extended Resonant Background
  RooDataHist resBjetDataExt("resBjetDataExt", "resBjetDataExt", massBjetExt, resBjetMassExt);
  RooFitResult *resCBResult = resPDFBjetExt.fitTo(resBjetDataExt, RooFit::Extended(kTRUE), RooFit::Strategy(2));
  //Resonant Background
  RooCBShape resPDFBjetCB("resPDFBjetCB", "Resonant CB Background PDF bb", massBjet, resMeanBjet, resSigmaBjet, resAlpha, resPower);
  RooExponential resPDFExpo("resPDFExp", "Resonant Exponential Background PDF bb", massBjet, resExpo);
  RooAddPdf resPDFBjet("resPDFBjet", "Resonant Background PDF bb", RooArgList(resPDFBjetCB, resPDFExpo), RooArgList(nOthers, nbbH));
  //Non-Resonant background
  RooExponential nonresPDFBjet("nonresPDFBjet", "Non-Resonant Background PDF bb", massBjet, expRateBjet);
  
  //PDFs for diPhoton*diBjet (2D)
  RooProdPdf sig2Dpdf("sig2Dpdf", "2D Signal PDF", RooArgList(sigPDFPho, sigPDFBjet));
  RooProdPdf res2Dpdf("res2Dpdf", "2D Resonant Background PDF", RooArgList(resPDFPho, resPDFBjet));
  RooProdPdf nonres2Dpdf("nonres2Dpdf", "2D Non-Resonant Background PDF", RooArgList(nonresPDFPho, nonresPDFBjet));
  RooAddPdf model2Dpdf("model2Dpdf", "2D Signal+Background PDF", RooArgList(sig2Dpdf, res2Dpdf, nonres2Dpdf), RooArgList(nsig, nres, nnonres));
  
  //variables to set constant
  sigMeanBjet.setConstant();
  sigSigmaBjet.setConstant();
  sigAlpha.setConstant();
  sigPower.setConstant();
  resMeanBjet.setConstant();
  resSigmaBjet.setConstant();
  resAlpha.setConstant();
  resPower.setConstant();
  resExpo.setConstant();
  nbbH.setConstant();
  nOthers.setConstant();
  
  //Add all 2D stuff to workspace
  model2Dpdf.graphVizTree("fullModel.dot");
  ws->import(sigBjetData, Rename("sigBjetData"));
  ws->import(sigPDFBjetCut);
  ws->import(resBjetDataExt, Rename("resBjetDataExt"));
  ws->import(resPDFBjetExt);
  ws->import(model2Dpdf);
  
  cout << "#bar{M}_{bb} (signal) = " << sigMeanBjet.getVal() << endl;
  cout << "#sigma_{bb} (signal) = " << sigSigmaBjet.getVal() << endl;
  cout << "#bar{M}_{bb} (res) = " << resMeanBjet.getVal() << endl;
  cout << "#sigma_{bb} (res) = " << resSigmaBjet.getVal() << endl;
  cout << "#lambda_{bb} = " << expRateBjet.getVal() << endl;
  cout << "#bar{M}_{#gamma#gamma} =" << meanPho.getVal() << endl;
  cout << "#sigma_{#gamma#gamma} = " << sigmaPho.getVal() << endl;
  cout << "#lambda_{#gamma#gamma} = " << expRatePho.getVal() << endl;
}

void plotInitialFit(TH1F *hist, TF1 *func, const string outputName) {
  cv = new TCanvas("cv","cv", 800,600);
  cv->SetFillStyle(4000);
  hist->Draw("HIST"); 
  hist->SetStats(kFALSE);
  std::cout << hist->Integral() << std::endl;
  func->SetLineColor(kBlue);
  func->SetLineWidth(3);
  func->Draw("SAME");
  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.048);
  tex->SetTextFont(42);
  if (outputName == "sigPhoHistFit") {
  	tex->DrawLatex(0.64, 0.84, Form("#bar{M}_{#gamma#gamma} = %.2f", func->GetParameter(1)));
  	tex->DrawLatex(0.64, 0.79, Form("#sigma_{#gamma#gamma} = %.2f", func->GetParameter(2)));
  }
  else if (outputName == "resPhoHistFit") {
  	tex->DrawLatex(0.64, 0.84, Form("#bar{M}_{#gamma#gamma} = %.2f", func->GetParameter(1)));
  	tex->DrawLatex(0.64, 0.79, Form("#sigma_{#gamma#gamma} = %.2f", func->GetParameter(2)));
  }
  else if (outputName == "nonresPhoHistFitExpo") {
  	tex->DrawLatex(0.64, 0.84, Form("#lambda_{#gamma#gamma} = %.3f", func->GetParameter(1)));
  }
  else if (outputName == "sigBjetHistFit") {
  	tex->DrawLatex(0.64, 0.84, Form("#bar{M}_{bb} = %.2f", func->GetParameter(1)));
  	tex->DrawLatex(0.64, 0.79, Form("#sigma_{bb} = %.2f", func->GetParameter(2)));
  }
  else if (outputName == "resBjetHistExtFit") {
  	tex->DrawLatex(0.64, 0.84, Form("#bar{M}_{bb} = %.2f", func->GetParameter(1)));
  	tex->DrawLatex(0.64, 0.79, Form("#sigma_{bb} = %.2f", func->GetParameter(2)));
  }
  else if (outputName == "nonresBjetHistFitExpo") {
  	tex->DrawLatex(0.64, 0.84, Form("#lambda_{bb} = %.3f", func->GetParameter(1)));
  }
  tex->Draw();
  cv->Update();
  cv->SaveAs(("Plots/AllSignalBkgd/Fits/InitialFits/"+outputName+".gif").c_str());
}

//-------------------------------------------------------------
//Plot the model fitting function for signal+background
//=============================================================
void MakeCBPlots(RooWorkspace *ws) {
  
  RooPlot* framex = 0;
  RooPlot* framey = 0;
  
  //Import the PDF's
  RooAbsPdf *sigPDFBjetCut = ws->pdf("sigPDFBjetCut");
  RooAbsPdf *resPDFBjetExt = ws->pdf("resPDFBjetExt");
  //x-axis variables
  RooRealVar *massBjetExt = ws->var("massBjetExt");
  RooRealVar *massBjetCut = ws->var("massBjetCut");
  //sig variables
  RooRealVar *sigMeanBjet = ws->var("sigMeanBjet");
  RooRealVar *sigSigmaBjet = ws->var("sigSigmaBjet");
  RooRealVar *sigAlpha = ws->var("sigAlpha");
  RooRealVar *sigPower = ws->var("sigPower");
  //res bkg variables
  RooRealVar *resMeanBjet = ws->var("resMeanBjet");
  RooRealVar *resSigmaBjet = ws->var("resSigmaBjet");
  RooRealVar *resAlpha = ws->var("resAlpha");
  RooRealVar *resPower = ws->var("resPower");
  RooRealVar *resExpo = ws->var("resExpo");
  RooRealVar *nbbH = ws->var("nbbH");
  RooRealVar *nOthers = ws->var("nOthers");
  //simulated data
  RooDataHist *sigBjetData = (RooDataHist *)ws->data("sigBjetData");
  RooDataHist *resBjetDataExt = (RooDataHist *)ws->data("resBjetDataExt");
  
  //Plot the Crystal Ball Fit for Signal Bjet
  cv = new TCanvas("cv","cv",800,600);
  framex = massBjetCut->frame(Bins(50));
  sigBjetData->plotOn(framex, DrawOption("B"), DataError(RooAbsData::None), FillColor(kRed));
  sigPDFBjetCut->plotOn(framex);
  framex->SetTitle(""); framex->SetXTitle("M_{bb} [GeV/c^{2}]");  framex->SetYTitle("Number of Events");
  framex->Draw();
  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.048);
  tex->SetTextFont(42);
  tex->DrawLatex(0.64, 0.84, Form("#bar{M}_{bb} = %.2f", sigMeanBjet->getVal()));
  tex->DrawLatex(0.64, 0.79, Form("#sigma_{bb} = %.2f", sigSigmaBjet->getVal()));
  //tex->DrawLatex(0.64, 0.74, Form("#alpha_{bb} = %.2f", sigAlpha->getVal()));
  //tex->DrawLatex(0.64, 0.69, Form("n_{bb} = %.2f", sigPower->getVal()));
  tex->Draw();
  cv->Update();
  cv->SaveAs("Plots/AllSignalBkgd/Fits/InitialFits/sigBjetHistCutFitCB.gif");
  
  //Plot the Crystal Ball Fit for Bjet Resonant Background
  cv = new TCanvas("cv","cv",800,600);
  framex = massBjetExt->frame(Bins(50));
  resBjetDataExt->plotOn(framex, DrawOption("B"), DataError(RooAbsData::None), FillColor(kOrange));
  resPDFBjetExt->plotOn(framex);
  framex->SetTitle(""); framex->SetXTitle("M_{bb} [GeV/c^{2}]");  framex->SetYTitle("Number of Events");
  framex->Draw();
  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.048);
  tex->SetTextFont(42);
  tex->DrawLatex(0.64, 0.84, Form("#bar{M}_{bb} = %.2f", resMeanBjet->getVal()));
  tex->DrawLatex(0.64, 0.79, Form("#sigma_{bb} = %.2f", resSigmaBjet->getVal()));
  //tex->DrawLatex(0.64, 0.74, Form("#alpha_{bb} = %.2f", resAlpha->getVal()));
  //tex->DrawLatex(0.64, 0.69, Form("n_{bb} = %.2f", resPower->getVal()));
  tex->Draw();
  cv->Update();
  cv->SaveAs("Plots/AllSignalBkgd/Fits/InitialFits/resBjetHistExtFitCB.gif");
}



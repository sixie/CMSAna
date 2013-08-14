//
// root -l CMSAna/HHToBBGG/fits/FitMassOneD.C+'("/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_dibjetMassPhoWin.root", 0, 0, 1)'
//
//
// for plotOption argument: 0 = don't plot, 1 = plot and save
// for storeOption argument: 0 = don't store, 1 = store and save
//
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
int alreadyPlotted = 0;

void plotInitialFit(TH1F *hist, TF1 *func, const string outputName);
void AddModels(RooWorkspace *ws, const string inputfileBjet, Int_t plotOption);
void MakePlots(RooWorkspace *ws, RooFitResult *fitResultBjet);

void FitMassOneD(const string inputfileBjet = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_dibjetMassPhoWin.root", Int_t plotOption = 1, Int_t storeOption = 1, Int_t NToys = 10000) {

  TRandom3 *randomnumber = new TRandom3(1200);
  RooWorkspace* ws = new RooWorkspace("CMSAna/HHToBBGG/fits/Workspace");
  AddModels(ws, inputfileBjet, plotOption);
  
  //Import variables from workspace
  RooAbsPdf *modelBjet = ws->pdf("modelBjet");
  RooRealVar *massBjet = ws->var("massBjet");
  RooRealVar *nsigBjet = ws->var("N_{bb} (Sig)"); RooRealVar constNsig(*nsigBjet);
  RooRealVar *nresBjet = ws->var("N_{bb} (ResBkg)"); RooRealVar constNres(*nresBjet);
  RooRealVar *nnonresBjet = ws->var("N_{bb} (NonResBkg)"); RooRealVar constNnonres(*nnonresBjet);
  RooRealVar *expRateBjet = ws->var("expRateBjet"); RooRealVar constexpBjet(*expRateBjet);
  
  //Variables to set constant
  RooRealVar *sigMeanBjet = ws->var("sigMeanBjet"); sigMeanBjet->setConstant();
  RooRealVar *sigSigmaBjet = ws->var("sigSigmaBjet"); sigSigmaBjet->setConstant();
  RooRealVar *sigAlpha = ws->var("sigAlpha"); sigAlpha->setConstant();
  RooRealVar *sigPower = ws->var("sigPower"); sigPower->setConstant();
  RooRealVar *resMeanBjet = ws->var("resMeanBjet"); resMeanBjet->setConstant();
  RooRealVar *resSigmaBjet = ws->var("resSigmaBjet"); resSigmaBjet->setConstant();
  RooRealVar *resAlpha = ws->var("resAlpha"); resAlpha->setConstant();
  RooRealVar *resPower = ws->var("resPower"); resPower->setConstant();
  
  //Create TTree to store the resulting yield data
  TFile *f = new TFile("CMSAna/HHToBBGG/data/MassFitResults/MassFitOneD.root","RECREATE");
  TTree *resultTree = new TTree("resultTree", "Parameter results from fitting");
  Float_t nsigOut, nresOut, nnonresOut;
  Float_t nsigStd, nresStd, nnonresStd;
  
  resultTree->Branch("nsigOut",&nsigOut, "nsigOut/F");
  resultTree->Branch("nresOut",&nresOut, "nresOut/F");
  resultTree->Branch("nnonresOut",&nnonresOut, "nnonresOut/F");
  resultTree->Branch("nsigStd",&nsigStd, "nsigStd/F");
  resultTree->Branch("nresStd",&nresStd, "nresStd/F");
  resultTree->Branch("nnonresStd",&nnonresStd, "nnonresStd/F");
  
  //-------------------------------------------------------------
  //Generate Toy MC experiment data and fit
  //=============================================================  		
  for(UInt_t t=0; t < NToys; ++t) {
    
    nsigBjet->setVal(constNsig.getVal()); nresBjet->setVal(constNres.getVal()); nnonresBjet->setVal(constNnonres.getVal());
    expRateBjet->setVal(constexpBjet.getVal());
    RooDataSet *pseudoDataBjet = modelBjet->generate(RooArgList(*massBjet), randomnumber->Poisson(95.76));//45.66, 95.76
    RooFitResult *fitResultBjet = modelBjet->fitTo(*pseudoDataBjet, RooFit::Save(), RooFit::Extended(kTRUE), RooFit::Strategy(2));
    ws->import(*pseudoDataBjet, kTRUE);
    ws->import(*pseudoDataBjet, Rename("pseudoDataBjet"));
    
    if (plotOption == 1) MakePlots(ws, fitResultBjet);
    //-------------------------------------------------------------
    //Store fit parameters into ROOT file
    //=============================================================
    if (storeOption == 1) {
      //Save variables into separate branches of root tree
  		nsigOut = nsigBjet->getVal();
  		nresOut = nresBjet->getVal();
  		nnonresOut = nnonresBjet->getVal();
  		nsigStd = nsigBjet->getPropagatedError(*fitResultBjet);
  		nresStd = nresBjet->getPropagatedError(*fitResultBjet);
  		nnonresStd = nnonresBjet->getPropagatedError(*fitResultBjet);
  		
  		resultTree->Fill();
    }
  }
  
  //Write to the TTree and close it
  resultTree->Write();
  f->Close();
  delete ws;
  
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
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.032);
  tex->SetTextFont(2);
  if (outputName == "sigBjetHistFitPhoWin") {
  	tex->DrawLatex(0.72, 0.86, Form("#bar{M}_{bb} = %.2f", func->GetParameter(1)));
  	tex->DrawLatex(0.72, 0.82, Form("#sigma_{bb} = %.2f", func->GetParameter(2)));
  }
  else if (outputName == "resBjetHistExtFitPhoWin") {
  	tex->DrawLatex(0.72, 0.86, Form("#bar{M}_{bb} = %.2f", func->GetParameter(1)));
  	tex->DrawLatex(0.72, 0.82, Form("#sigma_{bb} = %.2f", func->GetParameter(2)));
  }
  else if (outputName == "nonresBjetHistFitExpoPhoWin") {
  	tex->DrawLatex(0.72, 0.86, Form("#lambda_{bb} = %.3f", func->GetParameter(1)));
  }
  tex->Draw();
  cv->Update();
  cv->SaveAs(("Plots/AllSignalBkgd/Fits/"+outputName+".gif").c_str());
}

void AddModels(RooWorkspace *ws, const string inputfileBjet, Int_t plotOption) {

  gStyle->SetOptFit(0111);
  
  //-------------------------------------------------------------
  //Make Histograms
  //=============================================================

  TFile *fileBjet = new TFile(inputfileBjet.c_str(), "READ");
  TH1F *sigBjetMass = (TH1F*)fileBjet->Get("sigBjetMass");
  TH1F *resBjetMass = (TH1F*)fileBjet->Get("resBjetMass");
  TH1F *resBjetMassExt = (TH1F*)fileBjet->Get("resBjetMassExtPhoWin");
  TH1F *nonresBjetMass = (TH1F*)fileBjet->Get("nonresBjetMass");
  
  //-------------------------------------------------------------
  //Fit Signal + Background and store parameters
  //=============================================================
  
  //fitting dibjet signal and background with diphoton mass window of 120-130 GeV
  TF1 *sigBjetFunc = new TF1("Signal diBjet Function","gaus",70,200);
  TF1 *resBjetFuncExt = new TF1("Resonant Background diBjet Function","gaus",40,200);
  TF1 *nonresBjetExpo = new TF1("Non-Resonant Background diBjet Function","expo",70,200);
  sigBjetMass->Fit(sigBjetFunc,"R");
  resBjetMass->Fit(resBjetFuncExt,"R");
  nonresBjetMass->Fit(nonresBjetExpo,"R");
  
  //storing parameters
  Double_t parSigBjet[3]; Double_t parBkgBjet[5];
  sigBjetFunc->GetParameters(&parSigBjet[0]);
  resBjetFuncExt->GetParameters(&parBkgBjet[0]);
  nonresBjetExpo->GetParameters(&parBkgBjet[3]);
  
  if (plotOption == 1) {
    plotInitialFit(sigBjetMass,sigBjetFunc,"sigBjetHistFitPhoWin");
  	plotInitialFit(resBjetMassExt,resBjetFuncExt,"resBjetHistExtFitPhoWin");
  	plotInitialFit(nonresBjetMass,nonresBjetExpo,"nonresBjetHistFitExpoPhoWin");
  
  }

  //-------------------------------------------------------------
  //Make Variables for Fits
  //=============================================================
  
  //Variables for 1D diBjet with diPhoton Mass Window of 120-130 GeV
  RooRealVar massBjet("massBjet","M_{bb}",70.,200.,"GeV/c^{2}");
  RooRealVar massBjetExt("massBjetExt","M_{bb} (Ext)",40.,200.,"GeV/c^{2}");
  //Signal diBjet
  RooRealVar sigMeanBjet("sigMeanBjet", "#bar{M}_{bb}", parSigBjet[1], 119., 131., "GeV/c^{2}");
  RooRealVar sigSigmaBjet("sigSigmaBjet", "#sigma_{bb}", parSigBjet[2], 12., 20., "GeV/c^{2}");
  RooRealVar sigAlpha("sigAlpha","#alpha_{bb} (cut)", 1., 0.5, 1.5);
  RooRealVar sigPower("sigPower", "n_{bb} (power)", 3., 1., 5.);
  //Resonant Background diBjet
  RooRealVar resMeanBjet("resMeanBjet", "#bar{M}_{bb} (ResBkg)", parBkgBjet[1], 85., 90., "GeV/c^{2}");
  RooRealVar resSigmaBjet("resSigmaBjet", "#sigma_{bb} (ResBkg)", parBkgBjet[2], 3.5, 6.5, "GeV/c^{2}");
  RooRealVar resAlpha("resAlpha","#alpha_{bb} (ResBkg) (cut)", -1.1, -2., -0.1);
  RooRealVar resPower("resPower", "n_{bb} (ResBkg) (power)", 3., 1., 5.);
  //Non-Resonant Background diBjet
  RooRealVar expRateBjet("expRateBjet", "#lambda_{bb}", parBkgBjet[4], -0.1, 0.0);
  //Weights for 1D fit
  RooRealVar nsigBjet("N_{bb} (Sig)", "# signal events bb", 15.85, 0, 1000); //11.86, 15.85
  RooRealVar nresBjet("N_{bb} (ResBkg)", "# resonant background events bb", 26.2, 0, 1000); //16.7, 26.2
  RooRealVar nnonresBjet("N_{bb} (NonResBkg)", "# non-resonant background events bb", 53.7, 0, 1000); //17.1, 53.7

  //-------------------------------------------------------------
  //Make PDFs for Signal + Background
  //=============================================================
  
  //PDFs for 1D diBjet with diPhoton Mass Window of 120-130 GeV
  RooCBShape sigPDFBjet("sigPDFBjet", "Signal PDF bb", massBjet, sigMeanBjet, sigSigmaBjet, sigAlpha, sigPower);
  RooDataHist sigBjetData("sigBjetData", "sigBjetData", massBjet, sigBjetMass);
  RooFitResult *CBResult = sigPDFBjet.fitTo(sigBjetData);
  RooCBShape resPDFBjetExt("resPDFBjetExt", "Extended Resonant Background PDF bb", massBjetExt, resMeanBjet, resSigmaBjet, resAlpha, resPower);
  RooDataHist resBjetDataExt("resBjetDataExt", "resBjetDataExt", massBjetExt, resBjetMassExt);
  RooFitResult *resCBResult = resPDFBjetExt.fitTo(resBjetDataExt);
  RooCBShape resPDFBjet("resPDFBjet", "Resonant Background PDF bb", massBjet, resMeanBjet, resSigmaBjet, resAlpha, resPower);
  RooExponential nonresPDFBjet("nonresPDFBjet", "Non-Resonant Background PDF bb", massBjet, expRateBjet);
  
  RooAddPdf modelBjet("modelBjet", "Signal+Background PDF bb", RooArgList(sigPDFBjet, resPDFBjet, nonresPDFBjet), RooArgList(nsigBjet,nresBjet,nnonresBjet));
    
  //Add all 1D stuff to workspace
  modelBjet.graphVizTree("fullModel.dot");
  ws->import(sigBjetData, Rename("sigBjetData"));
  ws->import(resBjetDataExt, Rename("resBjetDataExt"));
  ws->import(resPDFBjetExt);
  ws->import(modelBjet);
}


//-------------------------------------------------------------
//Plot the model fitting function for signal+background
//=============================================================
void MakePlots(RooWorkspace *ws, RooFitResult *fitResultBjet) {
  
  //For this to actually work, need to initialize a bunch of variables.... have fun
  RooPlot *plot = 0;
  RooPlot* framex = 0;
  RooPlot* framey = 0;
  
  //Import yield variables
  RooRealVar *nsigBjet = ws->var("N_{bb} (Sig)");
  RooRealVar *nresBjet = ws->var("N_{bb} (ResBkg)");
  RooRealVar *nnonresBjet = ws->var("N_{bb} (NonResBkg)");
  //Select and plot only one experiment
  if (! (nsigBjet->getVal() >= 15.5 && nsigBjet->getVal() <= 16.2
  				&& nresBjet->getVal() >= 25.8 && nresBjet->getVal() <= 26.5))
  				//&& nnonresBjet->getVal() >= 16.9 && nresBjet->getVal() <= 17.3))
  	return;
  if (alreadyPlotted) return;
  alreadyPlotted = 1;
  
  RooAbsPdf *modelBjet = ws->pdf("modelBjet");
  RooAbsPdf *sigPDFBjet = ws->pdf("sigPDFBjet");
  RooAbsPdf *resPDFBjet = ws->pdf("resPDFBjet");
  RooAbsPdf *resPDFBjetExt = ws->pdf("resPDFBjetExt");
  RooAbsPdf *nonresPDFBjet = ws->pdf("nonresPDFBjet");
  //x-axis variables
  RooRealVar *massBjet = ws->var("massBjet");
  RooRealVar *massBjetExt = ws->var("massBjetExt");
  //signal variables
  RooRealVar *sigMeanBjet = ws->var("sigMeanBjet");
  RooRealVar *sigSigmaBjet = ws->var("sigSigmaBjet");
  RooRealVar *sigAlpha = ws->var("sigAlpha");
  RooRealVar *sigPower = ws->var("sigPower");
  //res bkg variables
  RooRealVar *resMeanBjet = ws->var("resMeanBjet");
  RooRealVar *resSigmaBjet = ws->var("resSigmaBjet");
  RooRealVar *resAlpha = ws->var("resAlpha");
  RooRealVar *resPower = ws->var("resPower");
  //simulated data
  RooDataHist *sigBjetData = (RooDataHist *)ws->data("sigBjetData");
  RooDataHist *resBjetDataExt = (RooDataHist *)ws->data("resBjetDataExt");
  RooDataSet *pseudoDataBjet = (RooDataSet *)ws->data("pseudoDataBjet");
  
  //Plot the Crystal Ball Fit for Signal Bjet
  cv = new TCanvas("cv","cv",800,600);
  framex = massBjet->frame(Bins(50));
  sigBjetData->plotOn(framex, DrawOption("B"), DataError(RooAbsData::None), FillColor(kRed));
  sigPDFBjet->plotOn(framex);
  framex->Draw();
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.032);
  tex->SetTextFont(2);
  tex->DrawLatex(0.72, 0.86, Form("#bar{M}_{bb} = %.2f", sigMeanBjet->getVal()));
  tex->DrawLatex(0.72, 0.82, Form("#sigma_{bb} = %.2f", sigSigmaBjet->getVal()));
  tex->DrawLatex(0.72, 0.78, Form("#alpha_{bb} = %.2f", sigAlpha->getVal()));
  tex->DrawLatex(0.72, 0.74, Form("n_{bb} = %.2f", sigPower->getVal()));
  tex->Draw();
  cv->Update();
  cv->SaveAs("Plots/AllSignalBkgd/Fits/sigBjetHistFitCBPhoWin.gif");
  
  //Plot the Crystal Ball Fit for Bjet Resonant Background
  cv = new TCanvas("cv","cv",800,600);
  framex = massBjetExt->frame(Bins(50));
  resBjetDataExt->plotOn(framex, DrawOption("B"), DataError(RooAbsData::None), FillColor(kOrange));
  resPDFBjetExt->plotOn(framex);
  framex->Draw();
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.032);
  tex->SetTextFont(2);
  tex->DrawLatex(0.72, 0.86, Form("#bar{M}_{bb} = %.2f", resMeanBjet->getVal()));
  tex->DrawLatex(0.72, 0.82, Form("#sigma_{bb} = %.2f", resSigmaBjet->getVal()));
  tex->DrawLatex(0.72, 0.78, Form("#alpha_{bb} = %.2f", resAlpha->getVal()));
  tex->DrawLatex(0.72, 0.74, Form("n_{bb} = %.2f", resPower->getVal()));
  tex->Draw();
  cv->Update();
  cv->SaveAs("Plots/AllSignalBkgd/Fits/resBjetExtHistFitCBPhoWin.gif");
  
  //Plot of 1D fits for massBjet
  framex = massBjet->frame(Bins(35));
  pseudoDataBjet->plotOn(framex);
  modelBjet->plotOn(framex, VisualizeError(*fitResultBjet), FillStyle(3001));
  modelBjet->plotOn(framex);
  modelBjet->plotOn(framex,Components("sigPDFBjet"), LineStyle(kDashed), LineColor(kRed));
  modelBjet->plotOn(framex,Components("resPDFBjet"), LineStyle(kDashed), LineColor(kOrange));
  modelBjet->plotOn(framex,Components("nonresPDFBjet"), LineStyle(kDashed), LineColor(kGreen));
    
  cv = new TCanvas("cv","cv",800,600);
  framex->Draw();
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.024);
  tex->SetTextFont(2);
  tex->DrawLatex(0.62, 0.86, Form("N_{bb} (Sig) = %.2f +/- %.2f", nsigBjet->getVal(), nsigBjet->getPropagatedError(*fitResultBjet)));
  tex->DrawLatex(0.62, 0.82, Form("N_{bb} (ResBkg) = %.2f +/- %.2f", nresBjet->getVal(), nresBjet->getPropagatedError(*fitResultBjet)));
  tex->DrawLatex(0.62, 0.78, Form("N_{bb} (NonResBkg) = %.2f +/- %.2f", nnonresBjet->getVal(), nnonresBjet->getPropagatedError(*fitResultBjet)));
  tex->Draw();
  cv->Update();
  cv->SaveAs(Form("Plots/AllSignalBkgd/Fits/oneDimFitPhoWin_%d.gif",0));
  
}



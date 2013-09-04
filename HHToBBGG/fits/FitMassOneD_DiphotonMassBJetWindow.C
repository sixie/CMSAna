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

//  void MakePlots(RooWorkspace *ws, RooFitResult *fitResultBjet);

void FitMassOneD_DiphotonMassBJetWindow(const string workspaceFile = "CMSAna/HHToBBGG/data/FitWorkspace_DiphotonMassInBJetWindow.root", Int_t plotOption = 1, Int_t storeOption = 1, Int_t NToys = 5000) {

  TRandom3 *randomnumber = new TRandom3(1200);
  TFile *wsFile = new TFile (workspaceFile.c_str(), "READ");
  RooWorkspace *ws = (RooWorkspace*)wsFile->Get("MassFitWorkspace");
  
  //Import variables from workspace
  RooAbsPdf *modelpdf = ws->pdf("modelpdf");
  RooRealVar *massPho = ws->var("massPho");
  RooRealVar *nsig = ws->var("N (Sig)"); RooRealVar constNsig(*nsig);
  RooRealVar *nnonres = ws->var("N (NonResBkg)"); RooRealVar constNnonres(*nnonres);
  RooRealVar *expRatePho = ws->var("expRatePho"); RooRealVar constexpPho(*expRatePho);

  //for floating bkg exponential parameter
  //expRatePho->setConstant(false);
  
  //Create TTree to store the resulting yield data
  TFile *f = new TFile("CMSAna/HHToBBGG/data/MassFitResults/MassFitOneD_DiphotonMassBJetWindow.root","RECREATE");
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
    
    nsig->setVal(constNsig.getVal()); nnonres->setVal(constNnonres.getVal());
    RooDataSet *pseudoData = modelpdf->generate(RooArgList(*massPho), randomnumber->Poisson(118.2));
    RooFitResult *fitResult = modelpdf->fitTo(*pseudoData, RooFit::Save(), RooFit::Extended(kTRUE), RooFit::Strategy(2));
    ws->import(*pseudoData, kTRUE);
    ws->import(*pseudoData, Rename("pseudoData"));
    
//     if (plotOption == 1) MakePlots(ws, fitResult);
    //-------------------------------------------------------------
    //Store fit parameters into ROOT file
    //=============================================================
    if (storeOption == 1) {
      //Save variables into separate branches of root tree
      nsigOut = nsig->getVal();
      nresOut = 0;
      nnonresOut = nnonres->getVal();
      nsigStd = nsig->getPropagatedError(*fitResult);
      nresStd = 0;
      nnonresStd = nnonres->getPropagatedError(*fitResult);
      
      resultTree->Fill();
    }
  }
  
  //Write to the TTree and close it
  resultTree->Write();
  f->Close();
  delete ws;
  
}




// //-------------------------------------------------------------
// //Plot the model fitting function for signal+background
// //=============================================================
// void MakePlots(RooWorkspace *ws, RooFitResult *fitResult) {
  
//   //For this to actually work, need to initialize a bunch of variables.... have fun
//   RooPlot *plot = 0;
//   RooPlot* framex = 0;
//   RooPlot* framey = 0;
  
//   //Import yield variables
//   RooRealVar *nsig = ws->var("NSig");
//   RooRealVar *nnonres = ws->var("NNonResBkg");
//   //Select and plot only one experiment

// //   if (! (nsigBjet->getVal() >= 15.5 && nsigBjet->getVal() <= 16.2
// //   				&& nresBjet->getVal() >= 25.8 && nresBjet->getVal() <= 26.5))
// //   				//&& nnonresBjet->getVal() >= 16.9 && nresBjet->getVal() <= 17.3))
// //   	return;

//   if (alreadyPlotted) return;
//   alreadyPlotted = 1;
  
//   RooAbsPdf *modelBjet = ws->pdf("modelBjet");
//   RooAbsPdf *sigPDFBjet = ws->pdf("sigPDFBjet");
//   RooAbsPdf *resPDFBjet = ws->pdf("resPDFBjet");
//   RooAbsPdf *resPDFBjetExt = ws->pdf("resPDFBjetExt");
//   RooAbsPdf *nonresPDFBjet = ws->pdf("nonresPDFBjet");
//   //x-axis variables
//   RooRealVar *massBjet = ws->var("massBjet");
//   RooRealVar *massBjetExt = ws->var("massBjetExt");
//   //signal variables
//   RooRealVar *sigMeanBjet = ws->var("sigMeanBjet");
//   RooRealVar *sigSigmaBjet = ws->var("sigSigmaBjet");
//   RooRealVar *sigAlpha = ws->var("sigAlpha");
//   RooRealVar *sigPower = ws->var("sigPower");
//   //res bkg variables
//   RooRealVar *resMeanBjet = ws->var("resMeanBjet");
//   RooRealVar *resSigmaBjet = ws->var("resSigmaBjet");
//   RooRealVar *resAlpha = ws->var("resAlpha");
//   RooRealVar *resPower = ws->var("resPower");
//   //simulated data
//   RooDataHist *sigBjetData = (RooDataHist *)ws->data("sigBjetData");
//   RooDataHist *resBjetDataExt = (RooDataHist *)ws->data("resBjetDataExt");
//   RooDataSet *pseudoDataBjet = (RooDataSet *)ws->data("pseudoDataBjet");
  
//   //Plot the Crystal Ball Fit for Signal Bjet
//   cv = new TCanvas("cv","cv",800,600);
//   framex = massBjet->frame(Bins(50));
//   sigBjetData->plotOn(framex, DrawOption("B"), DataError(RooAbsData::None), FillColor(kRed));
//   sigPDFBjet->plotOn(framex);
//   framex->Draw();
//   TLatex *tex = new TLatex();
//   tex->SetNDC();
//   tex->SetTextSize(0.032);
//   tex->SetTextFont(2);
//   tex->DrawLatex(0.72, 0.86, Form("#bar{M}_{bb} = %.2f", sigMeanBjet->getVal()));
//   tex->DrawLatex(0.72, 0.82, Form("#sigma_{bb} = %.2f", sigSigmaBjet->getVal()));
//   tex->DrawLatex(0.72, 0.78, Form("#alpha_{bb} = %.2f", sigAlpha->getVal()));
//   tex->DrawLatex(0.72, 0.74, Form("n_{bb} = %.2f", sigPower->getVal()));
//   tex->Draw();
//   cv->Update();
//   cv->SaveAs("Plots/AllSignalBkgd/Fits/sigBjetHistFitCBPhoWin.gif");
  
//   //Plot the Crystal Ball Fit for Bjet Resonant Background
//   cv = new TCanvas("cv","cv",800,600);
//   framex = massBjetExt->frame(Bins(50));
//   resBjetDataExt->plotOn(framex, DrawOption("B"), DataError(RooAbsData::None), FillColor(kOrange));
//   resPDFBjetExt->plotOn(framex);
//   framex->Draw();
//   TLatex *tex = new TLatex();
//   tex->SetNDC();
//   tex->SetTextSize(0.032);
//   tex->SetTextFont(2);
//   tex->DrawLatex(0.72, 0.86, Form("#bar{M}_{bb} = %.2f", resMeanBjet->getVal()));
//   tex->DrawLatex(0.72, 0.82, Form("#sigma_{bb} = %.2f", resSigmaBjet->getVal()));
//   tex->DrawLatex(0.72, 0.78, Form("#alpha_{bb} = %.2f", resAlpha->getVal()));
//   tex->DrawLatex(0.72, 0.74, Form("n_{bb} = %.2f", resPower->getVal()));
//   tex->Draw();
//   cv->Update();
//   cv->SaveAs("Plots/AllSignalBkgd/Fits/resBjetExtHistFitCBPhoWin.gif");
  
//   //Plot of 1D fits for massBjet
//   framex = massBjet->frame(Bins(35));
//   pseudoDataBjet->plotOn(framex);
//   modelBjet->plotOn(framex, VisualizeError(*fitResultBjet), FillStyle(3001));
//   modelBjet->plotOn(framex);
//   modelBjet->plotOn(framex,Components("sigPDFBjet"), LineStyle(kDashed), LineColor(kRed));
//   modelBjet->plotOn(framex,Components("resPDFBjet"), LineStyle(kDashed), LineColor(kOrange));
//   modelBjet->plotOn(framex,Components("nonresPDFBjet"), LineStyle(kDashed), LineColor(kGreen));
    
//   cv = new TCanvas("cv","cv",800,600);
//   framex->Draw();
//   TLatex *tex = new TLatex();
//   tex->SetNDC();
//   tex->SetTextSize(0.024);
//   tex->SetTextFont(2);
//   tex->DrawLatex(0.62, 0.86, Form("N_{bb} (Sig) = %.2f +/- %.2f", nsigBjet->getVal(), nsigBjet->getPropagatedError(*fitResultBjet)));
//   tex->DrawLatex(0.62, 0.82, Form("N_{bb} (ResBkg) = %.2f +/- %.2f", nresBjet->getVal(), nresBjet->getPropagatedError(*fitResultBjet)));
//   tex->DrawLatex(0.62, 0.78, Form("N_{bb} (NonResBkg) = %.2f +/- %.2f", nnonresBjet->getVal(), nnonresBjet->getPropagatedError(*fitResultBjet)));
//   tex->Draw();
//   cv->Update();
//   cv->SaveAs(Form("Plots/AllSignalBkgd/Fits/oneDimFitPhoWin_%d.gif",0));
  
// }



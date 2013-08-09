//
// root -l CMSAna/HHToBBGG/fits/FitMassTwoD.C+'("/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_diphotonMass.root", "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_dibjetMass.root", "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_twoDMass.root", 0, 0, 1)'
//
// for plotOption argument: 0 = don't plot, 1 = plot and save
// for storeOption argument: 0 = don't store, 1 = store and save
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

void plotInitialFit(TH1F *hist, TF1 *func, const string outputName);
void AddModels(RooWorkspace *ws, const string inputfilePho, const string inputfileBjet, const string inputfileTwoD, Int_t plotOption);
void MakePlots(RooWorkspace *ws, RooFitResult *fitResult2D);

void FitMassTwoD(const string inputfilePho = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_diphotonMass_tight.root", const string inputfileBjet = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_dibjetMass_tight.root", const string inputfileTwoD = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_twoDMass_tight.root", Int_t plotOption = 0, Int_t storeOption = 1, Int_t NToys = 10000) {

  TRandom3 *randomnumber = new TRandom3(1200);
  RooWorkspace* ws = new RooWorkspace("CMSAna/HHToBBGG/fits/Workspace");
  AddModels(ws, inputfilePho, inputfileBjet, inputfileTwoD, plotOption);
  
  //Import variables from workspace
  RooAbsPdf *model2Dpdf = ws->pdf("model2Dpdf");
  RooRealVar *massBjet = ws->var("massBjet");
  RooRealVar *massPho = ws->var("massPho");
  RooRealVar *nsig = ws->var("N (Sig)"); RooRealVar constNsig(*nsig);
  RooRealVar *nres = ws->var("N (ResBkg)"); RooRealVar constNres(*nres);
  RooRealVar *nnonres = ws->var("N (NonResBkg)"); RooRealVar constNnonres(*nnonres);
  RooRealVar *expRateBjet = ws->var("expRateBjet"); RooRealVar constexpBjet(*expRateBjet);
  RooRealVar *expRatePho = ws->var("expRatePho"); RooRealVar constexpPho(*expRatePho);
  
    
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
  TFile *f = new TFile("CMSAna/HHToBBGG/data/MassFitResults/MassFitTwoD_tight.root","RECREATE");
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
    
    nsig->setVal(constNsig.getVal()); nres->setVal(constNres.getVal()); nnonres->setVal(constNnonres.getVal());
    expRateBjet->setVal(constexpBjet.getVal()); expRatePho->setVal(constexpPho.getVal());
    RooDataSet *pseudoData2D = model2Dpdf->generate(RooArgList(*massBjet,*massPho), randomnumber->Poisson(125.5));
    RooFitResult *fitResult2D = model2Dpdf->fitTo(*pseudoData2D, RooFit::Save(), RooFit::Extended(kTRUE), RooFit::Strategy(2));
    ws->import(*pseudoData2D, Rename("pseudoData2D"));
    
    if (plotOption == 1) MakePlots(ws, fitResult2D);
    //-------------------------------------------------------------
    //Store fit parameters into ROOT file
    //=============================================================
    if (storeOption == 1) {
      //Save variables into separate branches of root tree
  		nsigOut = nsig->getVal();
  		nresOut = nres->getVal();
  		nnonresOut = nnonres->getVal();
  		nsigStd = nsig->getPropagatedError(*fitResult2D);
  		nresStd = nres->getPropagatedError(*fitResult2D);
  		nnonresStd = nnonres->getPropagatedError(*fitResult2D);
  		
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
  cv->SaveAs(("Plots/AllSignalBkgd/Fits/"+outputName+"_tight.gif").c_str());
}

void AddModels(RooWorkspace *ws, const string inputfilePho, const string inputfileBjet, const string inputfileTwoD, Int_t plotOption) {

  gStyle->SetOptFit(0111);
  
  //-------------------------------------------------------------
  //Make Histograms
  //=============================================================
  
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
  
  std::cout << "Correlation Factors {sig, res, nonres}: {" << sigMassTwoD->GetCorrelationFactor() << ", " << resMassTwoD->GetCorrelationFactor() << ", " << nonresMassTwoD->GetCorrelationFactor() << "}"<< std::endl;
  
  //-------------------------------------------------------------
  //Fit Signal + Background and store parameters
  //=============================================================
  
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
  nonresBjetExpo->GetParameters(&parBkgPho[3]);
  
  Double_t parSigBjet[3]; Double_t parBkgBjet[5];
  sigBjetFunc->GetParameters(&parSigBjet[0]);
  resBjetFuncExt->GetParameters(&parBkgBjet[0]);
  nonresBjetExpo->GetParameters(&parBkgBjet[3]);
  
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
  RooRealVar massBjetExt("massBjetExt","M_{bb} (Ext)",40.,200.,"GeV/c^{2}");
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
  //Non-Resonant Background diBjet
  RooRealVar expRateBjet("expRateBjet", "#lambda_{bb}", parBkgBjet[4], -0.1, 0.0);
  
  //Variables for diPhoton (everything fixed except for expo rate)
  RooRealVar massPho("massPho","M_{#gamma#gamma}",100.,150.,"GeV/c^{2}");
  RooRealVar meanPho("meanPho", "#bar{M}_{#gamma#gamma}", parSigPho[1], "GeV/c^{2}");
  RooRealVar sigmaPho("sigmaPho", "#sigma_{#gamma#gamma}", parSigPho[2], "GeV/c^{2}");
  RooRealVar expRatePho("expRatePho", "#lambda_{#gamma#gamma}", parBkgPho[4], -0.1, 0.0);
  
  //Weights for the 2D fit
  RooRealVar nsig("N (Sig)", "# signal events", 12.3, 0,100);
  RooRealVar nres("N (ResBkg)", "# resonant background events", 15.6, 0,100);
  RooRealVar nnonres("N (NonResBkg)", "# non-resonant background events", 95.6, 0,1000);

  //-------------------------------------------------------------
  //Make PDFs for Signal + Background
  //=============================================================
    
  //PDFs for diPhoton
  RooGaussian sigPDFPho("sigPDFPho", "Signal PDF #gamma#gamma", massPho, meanPho, sigmaPho);
  RooGaussian resPDFPho("resPDFPho", "Resonant Background PDF #gamma#gamma", massPho, meanPho, sigmaPho);
  RooExponential nonresPDFPho("nonresPDFPho", "Non-Resonant Background PDF #gamma#gamma", massPho, expRatePho);
    
  //PDFs for diBjet
  RooCBShape sigPDFBjet("sigPDFBjet", "Signal PDF bb", massBjet, sigMeanBjet, sigSigmaBjet, sigAlpha, sigPower);
  RooDataHist sigBjetData("sigBjetData", "sigBjetData", massBjet, sigBjetMass);
  RooFitResult *sigCBResult = sigPDFBjet.fitTo(sigBjetData);
  RooCBShape resPDFBjetExt("resPDFBjetExt", "Extended Resonant Background PDF bb", massBjetExt, resMeanBjet, resSigmaBjet, resAlpha, resPower);
  RooDataHist resBjetDataExt("resBjetDataExt", "resBjetDataExt", massBjetExt, resBjetMassExt);
  RooFitResult *resCBResult = resPDFBjetExt.fitTo(resBjetDataExt);
  RooCBShape resPDFBjet("resPDFBjet", "Resonant Background PDF bb", massBjet, resMeanBjet, resSigmaBjet, resAlpha, resPower);
  RooExponential nonresPDFBjet("nonresPDFBjet", "Non-Resonant Background PDF bb", massBjet, expRateBjet);
    
  //PDFs for diPhoton*diBjet (2D)
  RooProdPdf sig2Dpdf("sig2Dpdf", "2D Signal PDF", RooArgList(sigPDFPho, sigPDFBjet));
  RooProdPdf res2Dpdf("res2Dpdf", "2D Resonant Background PDF", RooArgList(resPDFPho, resPDFBjet));
  RooProdPdf nonres2Dpdf("nonres2Dpdf", "2D Non-Resonant Background PDF", RooArgList(nonresPDFPho, nonresPDFBjet));
  RooAddPdf model2Dpdf("model2Dpdf", "2D Signal+Background PDF", RooArgList(sig2Dpdf, res2Dpdf, nonres2Dpdf), RooArgList(nsig, nres, nnonres));
    
  //Add all 2D stuff to workspace
  model2Dpdf.graphVizTree("fullModel.dot");
  ws->import(sigBjetData, Rename("sigBjetData"));
  ws->import(resBjetDataExt, Rename("resBjetDataExt"));
  ws->import(resPDFBjetExt);
  ws->import(model2Dpdf);
}

//-------------------------------------------------------------
//Plot the model fitting function for signal+background
//=============================================================
void MakePlots(RooWorkspace *ws, RooFitResult *fitResult2D) {
  
  RooPlot *plot = 0;
  RooPlot* framex = 0;
  RooPlot* framey = 0;
  
  //Import the PDF's
  RooAbsPdf *model2Dpdf = ws->pdf("model2Dpdf");
  RooAbsPdf *sigPDFPho = ws->pdf("sigPDFPho");
  RooAbsPdf *resPDFPho = ws->pdf("resPDFPho");
  RooAbsPdf *nonresPDFPho = ws->pdf("nonresPDFPho");
  RooAbsPdf *sigPDFBjet = ws->pdf("sigPDFBjet");
  RooAbsPdf *resPDFBjet = ws->pdf("resPDFBjet");
  RooAbsPdf *resPDFBjetExt = ws->pdf("resPDFBjetExt");
  RooAbsPdf *nonresPDFBjet = ws->pdf("nonresPDFBjet");
  //x-axis variables
  RooRealVar *massPho = ws->var("massPho");
  RooRealVar *massBjet = ws->var("massBjet");
  RooRealVar *massBjetExt = ws->var("massBjetExt");
  //yields
  RooRealVar *nsig = ws->var("N (Sig)");
  RooRealVar *nres = ws->var("N (ResBkg)");
  RooRealVar *nnonres = ws->var("N (NonResBkg)");
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
  //simulated data
  RooDataHist *sigBjetData = (RooDataHist *)ws->data("sigBjetData");
  RooDataHist *resBjetDataExt = (RooDataHist *)ws->data("resBjetDataExt");
  RooDataSet *pseudoData2D = (RooDataSet *)ws->data("pseudoData2D");
  
  
  //Plot the Crystal Ball Fit for Signal Bjet
  cv = new TCanvas("cv","cv",800,600);
  framex = massBjet->frame(Bins(50));
  sigBjetData->plotOn(framex, DrawOption("B"), DataError(RooAbsData::None), FillColor(kRed));
  sigPDFBjet->paramOn(framex, Parameters(RooArgList(*sigMeanBjet, *sigSigmaBjet, *sigAlpha, *sigPower)), Format("NE"), Layout(0.59,0.89,0.89));
  sigPDFBjet->plotOn(framex);
  framex->Draw();
  cv->SaveAs("Plots/AllSignalBkgd/Fits/sigBjetHistFitCB_tight.gif");
  
  //Plot the Crystal Ball Fit for Bjet Resonant Background
  cv = new TCanvas("cv","cv",800,600);
  framex = massBjetExt->frame(Bins(50));
  resBjetDataExt->plotOn(framex, DrawOption("B"), DataError(RooAbsData::None), FillColor(kOrange));
  resPDFBjetExt->paramOn(framex, RooFit::Parameters(RooArgList(*resMeanBjet, *resSigmaBjet, *resAlpha, *resPower)), Format("NE"), Layout(0.59,0.89,0.89));
  resPDFBjetExt->plotOn(framex);
  framex->Draw();
  cv->SaveAs("Plots/AllSignalBkgd/Fits/resBjetHistExtFitCB_tight.gif");
  
  
  //Plot of 2D generated data and fits
  TH1 *data2d = pseudoData2D->createHistogram("2D Data", *massBjet,Binning(25), YVar(*massPho,Binning(25)));
  data2d->SetStats(0);
  TH1 *fit2d = model2Dpdf->createHistogram("2D Fit", *massBjet, YVar(*massPho));
  fit2d->SetStats(0);
  
  cv = new TCanvas("cv","cv",1600,600);
  cv->Divide(2);
  cv->cd(1); gPad->SetLeftMargin(0.15); data2d->Draw("LEGO2");
  data2d->GetXaxis()->SetTitleOffset(2);
  data2d->GetYaxis()->SetTitleOffset(2);
  data2d->GetZaxis()->SetTitleOffset(1.75);
  
  cv->cd(2); gPad->SetLeftMargin(0.15); fit2d->Draw("SURF1");
  fit2d->GetXaxis()->SetTitleOffset(2);
  fit2d->GetYaxis()->SetTitleOffset(2);
  fit2d->GetZaxis()->SetTitleOffset(1.75);
  cv->SaveAs(Form("Plots/AllSignalBkgd/Fits/twoDimensionalFits_%d_tight.gif",0));
  
  //Plot of massBjet and massPho projections from 2D fit
  framex = massBjet->frame(Bins(50)); 
  pseudoData2D->plotOn(framex);
  model2Dpdf->plotOn(framex, VisualizeError(*fitResult2D), FillStyle(3001));
  model2Dpdf->plotOn(framex);
  model2Dpdf->plotOn(framex,Components("sigPDFBjet"), LineStyle(kDashed), LineColor(kRed));
  model2Dpdf->plotOn(framex,Components("resPDFBjet"), LineStyle(kDashed), LineColor(kOrange));
  model2Dpdf->plotOn(framex,Components("nonresPDFBjet"), LineStyle(kDashed), LineColor(kGreen));
  
  framey = massPho->frame(Bins(50)); 
  pseudoData2D->plotOn(framey);
  model2Dpdf->plotOn(framey, VisualizeError(*fitResult2D), FillStyle(3001));
  model2Dpdf->plotOn(framey);
  model2Dpdf->plotOn(framey,Components("sigPDFPho"), LineStyle(kDashed), LineColor(kRed));
  model2Dpdf->plotOn(framey,Components("resPDFPho"), LineStyle(kDashed), LineColor(kOrange));
  model2Dpdf->plotOn(framey,Components("nonresPDFPho"), LineStyle(kDashed), LineColor(kGreen));
  
  
  cv = new TCanvas("cv","cv",1600,600);
  cv->Divide(2);
  cv->cd(1); framex->Draw();
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.024);
  tex->SetTextFont(2);
  tex->DrawLatex(0.64, 0.86, Form("N (Sig) = %.2f +/- %.2f", nsig->getVal(), nsig->getPropagatedError(*fitResult2D)));
  tex->DrawLatex(0.64, 0.82, Form("N (ResBkg) = %.2f +/- %.2f", nres->getVal(), nres->getPropagatedError(*fitResult2D)));
  tex->DrawLatex(0.64, 0.78, Form("N (NonResBkg) = %.2f +/- %.2f", nnonres->getVal(), nnonres->getPropagatedError(*fitResult2D)));
  tex->Draw();
  cv->Update();
  cv->cd(2); framey->Draw();
  tex->DrawLatex(0.64, 0.86, Form("N (Sig) = %.2f +/- %.2f", nsig->getVal(), nsig->getPropagatedError(*fitResult2D)));
  tex->DrawLatex(0.64, 0.82, Form("N (ResBkg) = %.2f +/- %.2f", nres->getVal(), nres->getPropagatedError(*fitResult2D)));
  tex->DrawLatex(0.64, 0.78, Form("N (NonResBkg) = %.2f +/- %.2f", nnonres->getVal(), nnonres->getPropagatedError(*fitResult2D)));
  tex->Draw();
  cv->Update();
  cv->SaveAs(Form("Plots/AllSignalBkgd/Fits/projectionFits_%d_tight.gif",0));
}



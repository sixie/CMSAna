//
// root -l CMSAna/HHToBBGG/fits/fitMassGG.C+'("/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_diphotonMass.root","/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_dibjetMass.root")'
//
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

void plotInitialFit(TH1F *hist, TF1 *func, const string outputName) {
  cv = new TCanvas("cv","cv", 800,600);
  cv->SetFillStyle(4000);
  hist->Draw("HIST");
  std::cout << hist->Integral() << std::endl;
  func->SetLineColor(kBlack);
  func->SetLineWidth(3);
  func->Draw("SAME");
  cv->SaveAs(("Plots/AllSignalBkgd/Fits/"+outputName+".gif").c_str());
}

void fitMassGG(const string inputfilePho = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_diphotonMass.root", const string inputfileBjet = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_dibjetMass.root", const string inputfileBjetPhoWin = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_dibjetMassPhoWin.root", const string inputfileTwoD = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_twoDMass.root", Int_t NToys = 1) {

  TRandom3 *randomnumber = new TRandom3(1200);
  TRandom3 *randomnumber1 = new TRandom3(1200);
  TRandom3 *randomnumber2 = new TRandom3(1200);
  gStyle->SetOptFit(0111);
  TFile *filePho = new TFile(inputfilePho.c_str(), "READ");
  TH1F *sigPhoMass = (TH1F*)filePho->Get("sigPhoMass");
  TH1F *bkgPhoMass = (TH1F*)filePho->Get("bkgPhoMass");
  TH1F *resPhoMass = (TH1F*)filePho->Get("resPhoMass");
  TH1F *nonresPhoMass = (TH1F*)filePho->Get("nonresPhoMass");
  TH1F *allPhoMass = (TH1F*)filePho->Get("allPhoMass");

  TFile *fileBjet = new TFile(inputfileBjet.c_str(), "READ");
  TH1F *sigBjetMass = (TH1F*)fileBjet->Get("sigBjetMass");
  TH1F *bkgBjetMass = (TH1F*)fileBjet->Get("bkgBjetMass");
  TH1F *resBjetMass = (TH1F*)fileBjet->Get("resBjetMass");
  TH1F *nonresBjetMass = (TH1F*)fileBjet->Get("nonresBjetMass");
  TH1F *allBjetMass = (TH1F*)fileBjet->Get("allBjetMass");

  TFile *fileBjetPhoWin = new TFile(inputfileBjetPhoWin.c_str(), "READ");
  TH1F *sigBjetMassPhoWin = (TH1F*)fileBjetPhoWin->Get("sigBjetMass");
  TH1F *resBjetMassPhoWin = (TH1F*)fileBjetPhoWin->Get("resBjetMass");
  TH1F *nonresBjetMassPhoWin = (TH1F*)fileBjetPhoWin->Get("nonresBjetMass");
  TH1F *allBjetMassPhoWin = (TH1F*)fileBjetPhoWin->Get("allBjetMass");
  
  TFile *fileTwoD = new TFile(inputfileTwoD.c_str(), "READ");
  TH2F *sigMassTwoD = (TH2F*)fileTwoD->Get("sigMassPhoBjet");
  TH2F *resMassTwoD = (TH2F*)fileTwoD->Get("resMassPhoBjet");
  TH2F *nonresMassTwoD = (TH2F*)fileTwoD->Get("nonresMassPhoBjet");
  
  std::cout << "Correlation Factors {sig, res, nonres}: {" << sigMassTwoD->GetCorrelationFactor() << ", " << resMassTwoD->GetCorrelationFactor() << ", " << nonresMassTwoD->GetCorrelationFactor() << "}"<< std::endl;
  //-------------------------------------------------------------
  //Fit Signal + Background and store parameters
  //=============================================================
  
  //fitting diphoton signal and background
  TF1 *sigPhoFunc = new TF1("Signal diPhoton Function","gaus",105,145);
  TF1 *resPhoFunc = new TF1("Resonant Background diPhoton Function","gaus",105,145);
  TF1 *nonresPhoExpo = new TF1("Non-Resonant Background diPhoton Function","expo",105,145);
  sigPhoMass->Fit(sigPhoFunc,"R");
  resPhoMass->Fit(resPhoFunc,"R");
  nonresPhoMass->Fit(nonresPhoExpo,"R");
  
  //fitting dibjet signal and background
  TF1 *sigBjetFunc = new TF1("Signal diBjet Function","gaus",70,190);
  TF1 *resBjetFunc = new TF1("Resonant Background diBjet Function","landau",70,190);
  TF1 *nonresBjetExpo = new TF1("Non-Resonant Background diBjet Function","expo",70,190);
  sigBjetMass->Fit(sigBjetFunc,"R");
  resBjetMass->Fit(resBjetFunc,"R");
  nonresBjetMass->Fit(nonresBjetExpo,"R");
  
  //fitting dibjet signal and background with diphoton mass window of 122-128 GeV
  TF1 *sigBjetFuncPhoWin = new TF1("Signal diBjet Function with diPhoton Mass Window","gaus",70,190);
  TF1 *resBjetFuncPhoWin = new TF1("Resonant Background diBjet Function with diPhoton Mass Window","landau",70,190);
  TF1 *nonresBjetExpoPhoWin = new TF1("Non-Resonant Background diBjet Function with diPhoton Mass Window","expo",70,190);
  sigBjetMassPhoWin->Fit(sigBjetFuncPhoWin,"R");
  resBjetMassPhoWin->Fit(resBjetFuncPhoWin,"R");
  nonresBjetMassPhoWin->Fit(nonresBjetExpoPhoWin,"R");
  
  //storing parameters
  Double_t parSigPho[3]; Double_t parBkgPho[9];
  sigPhoFunc->GetParameters(parSigPho);
  resPhoFunc->GetParameters(&parBkgPho[0]);
  nonresBjetExpo->GetParameters(&parBkgPho[3]);
  
  Double_t parSigBjet[6]; Double_t parBkgBjet[10];
  sigBjetFunc->GetParameters(&parSigBjet[0]);
  sigBjetFuncPhoWin->GetParameters(&parSigBjet[3]);
  resBjetFunc->GetParameters(&parBkgBjet[0]);
  nonresBjetExpo->GetParameters(&parBkgBjet[3]);
  resBjetFuncPhoWin->GetParameters(&parBkgBjet[5]);
  nonresBjetExpoPhoWin->GetParameters(&parBkgBjet[8]);
  
  /*
  plotInitialFit(sigPhoMass,sigPhoFunc,"sigPhoHistFit");
  plotInitialFit(resPhoMass,resPhoFunc,"resPhoHistFit");
  plotInitialFit(nonresPhoMass,nonresPhoExpo,"nonresPhoHistFitExpo");
  
  plotInitialFit(sigBjetMass,sigBjetFunc,"sigBjetHistFit");
  plotInitialFit(resBjetMass,resBjetFunc,"resBjetHistFit");
  plotInitialFit(nonresBjetMass,nonresBjetExpo,"nonresBjetHistFitExpo");
  
  plotInitialFit(sigBjetMassPhoWin,sigBjetFuncPhoWin,"sigBjetHistFitPhoWin");
  plotInitialFit(resBjetMassPhoWin,resBjetFuncPhoWin,"resBjetHistFitPhoWin");
  plotInitialFit(nonresBjetMassPhoWin,nonresBjetExpoPhoWin,"nonresBjetHistFitExpoPhoWin");
  */

  //-------------------------------------------------------------
  //Generate Toy MC experiment data and fit
  //=============================================================
  for(UInt_t t=0; t < NToys; ++t) {

    //-------------------------------------------------------------
    //Make PDFs for Signal + Background
    //=============================================================
    //RooWorkspace *w = new RooWorkspace("MassGGFitWorkspace");
  
    RooRealVar numEvents("numEvents","Number of Events",0.,30.);
  
    //Variables for diPhoton
    RooRealVar massPho("M_{#gamma#gamma}","M_{#gamma#gamma}",105.,145.,"GeV/c^{2}");
    RooRealVar meanPho("#bar{M}_{#gamma#gamma}", "Signal Mean #gamma#gamma", parSigPho[1], 123., 125.5, "GeV/c^{2}");
    RooRealVar sigSigmaPho("#sigma_{#gamma#gamma}", "Signal Width #gamma#gamma", parSigPho[2], 1.0, 2.0, "GeV/c^{2}");
    RooRealVar resSigmaPho("#sigma_{#gamma#gamma} (ResBkg)", "Resonant Width #gamma#gamma", parBkgPho[2], 1.1, 1.65, "GeV/c^{2}");
    RooRealVar expRatePho("#lambda_{#gamma#gamma}", "Exponential Rate #gamma#gamma", parBkgPho[4], -0.1, 0.0);
    RooRealVar nsigPho("N_{#gamma#gamma} (Sig)", "# signal events #gamma#gamma", 16.5, 12.5, 20.5);
    RooRealVar nresPho("N_{#gamma#gamma} (ResBkg)", "# resonant background events #gamma#gamma", 22., 17.3, 26.7);
    RooRealVar nnonresPho("N_{#gamma#gamma} (NonResBkg)", "# non-resonant background events #gamma#gamma", 146.5, 134.5, 158.5);
  
    //PDFs for diPhoton
    RooGaussian sigPDFPho("sigPDFPho", "Signal PDF #gamma#gamma", massPho, meanPho, sigSigmaPho);
    RooGaussian resPDFPho("resPDFPho", "Resonant Background PDF #gamma#gamma", massPho, meanPho, resSigmaPho);
    RooExponential nonresPDFPho("nonresPDFPho", "Non-Resonant Background PDF #gamma#gamma", massPho, expRatePho);
    RooAddPdf modelPho("modelPho", "Signal+Background PDF #gamma#gamma", RooArgList(sigPDFPho, resPDFPho, nonresPDFPho), RooArgList(nsigPho, nresPho, nnonresPho));
    
    //Variables for diBjet
    RooRealVar massBjet("M_{bb}","M_{bb}",70.,190.,"GeV/c^{2}");
    RooRealVar sigMeanBjet("#bar{M}_{bb}", "Signal Mean bb", parSigBjet[1], 119., 131., "GeV/c^{2}");
    RooRealVar sigSigmaBjet("#sigma_{bb}", "Signal Width bb", parSigBjet[2], 12., 20., "GeV/c^{2}");
    RooRealVar sigAlpha("#alpha_{bb}","Crystal Ball Cut", 1., 0.5, 1.5);
    RooRealVar sigPower("n_{bb}", "Crystal Ball Power", 3., 1., 5.);
    RooRealVar resMeanBjet("#bar{M}_{bb} (ResBkg)", "Resonant Mean bb", parBkgBjet[1], 85., 90., "GeV/c^{2}");
    RooRealVar resSigmaBjet("#sigma_{bb} (ResBkg)", "Resonant Width bb", parBkgBjet[2], 3.5, 6.5, "GeV/c^{2}");
    RooRealVar expRateBjet("#lambda_{bb}", "Exponential Rate bb", parBkgBjet[4], -0.1, 0.0);
    RooRealVar nsigBjet("N_{bb} (Sig)", "# signal events bb", 16.5, 12.5, 20.5);
    RooRealVar nresBjet("N_{bb} (ResBkg)", "# resonant background events bb", 22., 17.3, 26.7);
    RooRealVar nnonresBjet("N_{bb} (NonResBkg)", "# non-resonant background events bb", 146.5, 134.5, 158.5);
    
    //PDFs for diBjet
    RooCBShape sigPDFBjet("signalPDFBjet", "Signal PDF bb", massBjet, sigMeanBjet, sigSigmaBjet, sigAlpha, sigPower);
    RooDataHist sigBjetData("sigBjetData", "sigBjetData", massBjet, sigBjetMass);
    RooFitResult *CBResult = sigPDFBjet.fitTo(sigBjetData);
    RooLandau resPDFBjet("resPDFBjet", "Resonant Background PDF bb", massBjet, resMeanBjet, resSigmaBjet);
    RooExponential nonresPDFBjet("nonresPDFBjet", "Non-Resonant Background PDF bb", massBjet, expRateBjet);
    RooAddPdf modelBjet("modelBjet", "Signal+Background PDF bb", RooArgList(sigPDFBjet, resPDFBjet, nonresPDFBjet), RooArgList(nsigBjet,nresBjet,nnonresBjet));
    
    //Variables for diBjet with diPhoton Mass Window of 122-128 GeV
    RooRealVar massBjetPhoWin("M_{bb}","M_{bb}",70.,190.,"GeV/c^{2}");
    RooRealVar sigMeanBjetPhoWin("#bar{M}_{bb}", "Signal Mean bb", parSigBjet[4], 119., 131., "GeV/c^{2}");
    RooRealVar sigSigmaBjetPhoWin("#sigma_{bb}", "Signal Width bb", parSigBjet[5], 12., 20., "GeV/c^{2}");
    RooRealVar sigAlphaPhoWin("#alpha_{bb}","Crystal Ball Cut", 1., 0.5, 1.5);
    RooRealVar sigPowerPhoWin("n_{bb}", "Crystal Ball Power", 3., 1., 5.);
    RooRealVar resMeanBjetPhoWin("#bar{M}_{bb} (ResBkg)", "Resonant Mean bb", parBkgBjet[6], 85., 90., "GeV/c^{2}");
    RooRealVar resSigmaBjetPhoWin("#sigma_{bb} (ResBkg)", "Resonant Width bb", parBkgBjet[7], 3.5, 6.5, "GeV/c^{2}");
    RooRealVar expRateBjetPhoWin("#lambda_{bb}", "Exponential Rate bb", parBkgBjet[9], -0.1, 0.0);
    RooRealVar nsigBjetPhoWin("N_{bb} (Sig)", "# signal events bb", 14.6, 10.8, 18.4);
    RooRealVar nresBjetPhoWin("N_{bb} (ResBkg)", "# resonant background events bb", 18.2, 14.0, 22.4);
    RooRealVar nnonresBjetPhoWin("N_{bb} (NonResBkg)", "# non-resonant background events bb", 19.2, 14.8, 23.6);
    
    //PDFs for diBjet with diPhoton Mass Window of 122-128 GeV
    RooCBShape sigPDFBjetPhoWin("signalPDFBjet", "Signal PDF bb", massBjetPhoWin, sigMeanBjetPhoWin, sigSigmaBjetPhoWin, sigAlphaPhoWin, sigPowerPhoWin);
    RooDataHist sigBjetDataPhoWin("sigBjetData", "sigBjetData", massBjetPhoWin, sigBjetMassPhoWin);
    RooFitResult *CBResultPhoWin = sigPDFBjetPhoWin.fitTo(sigBjetDataPhoWin);
    RooLandau resPDFBjetPhoWin("resPDFBjet", "Resonant Background PDF bb", massBjetPhoWin, resMeanBjetPhoWin, resSigmaBjetPhoWin);
    RooExponential nonresPDFBjetPhoWin("nonresPDFBjet", "Non-Resonant Background PDF bb", massBjetPhoWin, expRateBjetPhoWin);
    RooAddPdf modelBjetPhoWin("modelBjet", "Signal+Background PDF bb", RooArgList(sigPDFBjetPhoWin, resPDFBjetPhoWin, nonresPDFBjetPhoWin), RooArgList(nsigBjetPhoWin,nresBjetPhoWin,nnonresBjetPhoWin));
    
    //PDFs for diPhoton*diBjet (2D)
    RooProdPdf model2Dpdf("model2Dpdf", "2D Signal+Background PDF", RooArgList(modelPho, modelBjet));
    
    
    //-------------------------------------------------------------
    //Generate data
    //=============================================================
  
    RooDataSet *pseudoDataPho = modelPho.generate(massPho, randomnumber->Poisson(163.));
    RooDataSet *pseudoDataBjetPhoWin = modelBjetPhoWin.generate(massBjetPhoWin, randomnumber1->Poisson(52.));
    RooDataSet *pseudoData2D = model2Dpdf.generate(RooArgList(massBjet,massPho), randomnumber2->Poisson(163.));
    //pseudoData->write(Form("Plots/toy%d.txt",t));

    RooFitResult *fitResultPho=0;
    RooFitResult *fitResultBjetPhoWin=0;
    RooFitResult *fitResult2D = 0;
    fitResultPho = modelPho.fitTo(*pseudoDataPho, RooFit::Save(), RooFit::Extended(kTRUE), RooFit::Strategy(2));
    fitResultBjetPhoWin = modelBjetPhoWin.fitTo(*pseudoDataBjetPhoWin, RooFit::Save(), RooFit::Extended(kTRUE), RooFit::Strategy(2));
    fitResult2D = model2Dpdf.fitTo(*pseudoData2D, RooFit::Save(), RooFit::Extended(kTRUE), RooFit::Strategy(2));
    
    //-------------------------------------------------------------
    //Store fit parameters into ROOT file
    //=============================================================
    /*
    TFile *outputfile = new TFile ("CMSAna/HHToBBGG/data/MassFitResults/Mass_Fit.root", "UPDATE");
    
    //declare variables of interest
    float massPhoArr[2] = {massPho.getVal(), massPho.getError()};
    float meanPhoArr[2] = {meanPho.getVal(), meanPho.getError()};
    float sigSigmaPhoArr[2] = {sigSigmaPho.getVal(), sigSigmaPho.getError()};
    float sigMeanBjetArr[2] = {sigMeanBjet.getVal(), sigMeanBjet.getError()};
    float sigSigmaBjetArr[2] = {sigSigmaBjet.getVal(), sigSigmaBjet.getError()};
    
    //save variables into separate branches of root tree
    TTree *outTree = new TTree(Form("fit_result_%d",t),Form("fit_result_%d",t));
    outTree->Branch("massPho",&massPhoArr, "Par/F:Err/F");
    outTree->Branch("meanPho",&meanPhoArr, "Par/F:Err/F");
    outTree->Branch("sigSigmaPho",&sigSigmaPhoArr, "Par/F:Err/F");
    outTree->Branch("sigMeanBjet",&sigMeanBjetArr, "Par/F:Err/F");
    outTree->Branch("sigSigmaBjet",&sigSigmaBjetArr, "Par/F:Err/F");
    
    outTree->Fill();
    outputfile->WriteTObject(outTree, outTree->GetName(), "WriteDelete");
    outputfile->Close();
    */
    //-------------------------------------------------------------
    //Plot the model fitting function for signal+background
    //=============================================================
    
    RooPlot *plot = 0;
    RooPlot* framex = 0;
    RooPlot* framey = 0;
    /*
    //Plot the Crystal Ball Fit for Signal Bjet
    cv = new TCanvas("cv","cv",800,600);
    framex = massBjet.frame(Bins(50));
    sigBjetData.plotOn(framex, DrawOption("B"), DataError(RooAbsData::None), FillColor(kRed));
    sigPDFBjet.paramOn(framex, Parameters(RooArgList(sigMeanBjet, sigSigmaBjet, sigAlpha, sigPower)), Format("NE"), Layout(0.59,0.89,0.89));
    sigPDFBjet.plotOn(framex);
    framex->Draw();
    cv->SaveAs("Plots/AllSignalBkgd/Fits/sigBjetHistFitCBPhoWin.gif");
    */
    
    //Plot of 1D fits for massBjet and massPho
    framex = massBjetPhoWin.frame(Bins(25));
    pseudoDataBjetPhoWin->plotOn(framex);
    modelBjetPhoWin.paramOn(framex, Parameters(RooArgList(nsigBjetPhoWin, nresBjetPhoWin, nnonresBjetPhoWin)), Format("NE"), Layout(0.59,0.89,0.89));
    framex->getAttText()->SetTextSize(.030);
    modelBjetPhoWin.plotOn(framex, VisualizeError(*fitResultBjetPhoWin), FillStyle(3001));
    modelBjetPhoWin.plotOn(framex);
    modelBjetPhoWin.plotOn(framex,Components(sigPDFBjetPhoWin), LineStyle(kDashed), LineColor(kRed));
    modelBjetPhoWin.plotOn(framex,Components(resPDFBjetPhoWin), LineStyle(kDashed), LineColor(kOrange));
    modelBjetPhoWin.plotOn(framex,Components(nonresPDFBjetPhoWin), LineStyle(kDashed), LineColor(kGreen));
    
    framey = massPho.frame(Bins(50));
    pseudoDataPho->plotOn(framey);
    modelPho.paramOn(framey, Parameters(RooArgList(nnonresPho, nsigPho, nresPho)), Format("NE"), Layout(0.59,0.89,0.89));
    framey->getAttText()->SetTextSize(.030);
    modelPho.plotOn(framey, VisualizeError(*fitResultPho), FillStyle(3001));
    modelPho.plotOn(framey);
    modelPho.plotOn(framey,Components(sigPDFPho), LineStyle(kDashed), LineColor(kRed));
    modelPho.plotOn(framey,Components(resPDFPho), LineStyle(kDashed), LineColor(kOrange));
    modelPho.plotOn(framey,Components(nonresPDFPho), LineStyle(kDashed), LineColor(kGreen));
    
    cv = new TCanvas("cv","cv",1600,600);
    cv->Divide(2); 
		cv->cd(1); framex->Draw();
		cv->cd(2); framey->Draw();
    cv->SaveAs(Form("Plots/AllSignalBkgd/Fits/test_oneDimensionalFits_%d.gif",t));
    
    //Plot of 2D generated data and fits
    TH1 *data2d = pseudoData2D->createHistogram("2D Data", massBjet,Binning(22), YVar(massPho,Binning(22)));
    data2d->SetStats(0);
    TH1 *fit2d = model2Dpdf.createHistogram("2D Fit", massBjet, YVar(massPho));
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
    cv->SaveAs(Form("Plots/AllSignalBkgd/Fits/test_twoDimensionalFits_%d.gif",t));
    
    //Plot of massBjet and massPho projections from 2D fit
    framex = massBjet.frame(Bins(50)); 
    pseudoData2D->plotOn(framex);
    model2Dpdf.paramOn(framex, Parameters(RooArgList(nsigBjet, nresBjet, nnonresBjet)), Format("NE"), Layout(0.59,0.89,0.89));
    framex->getAttText()->SetTextSize(.030);
    model2Dpdf.plotOn(framex, VisualizeError(*fitResult2D), FillStyle(3001));
		model2Dpdf.plotOn(framex);
    model2Dpdf.plotOn(framex,Components(sigPDFBjet), LineStyle(kDashed), LineColor(kRed));
    model2Dpdf.plotOn(framex,Components(resPDFBjet), LineStyle(kDashed), LineColor(kOrange));
    model2Dpdf.plotOn(framex,Components(nonresPDFBjet), LineStyle(kDashed), LineColor(kGreen));
		framex->Draw();
		
    framey = massPho.frame(Bins(50)); 
    pseudoData2D->plotOn(framey);
    model2Dpdf.paramOn(framey, Parameters(RooArgList(nnonresPho, nsigPho, nresPho)), Format("NE"), Layout(0.59,0.89,0.89));
    framey->getAttText()->SetTextSize(.030);
    model2Dpdf.plotOn(framey, VisualizeError(*fitResult2D), FillStyle(3001));
		model2Dpdf.plotOn(framey);
    model2Dpdf.plotOn(framey,Components(sigPDFPho), LineStyle(kDashed), LineColor(kRed));
    model2Dpdf.plotOn(framey,Components(resPDFPho), LineStyle(kDashed), LineColor(kOrange));
    model2Dpdf.plotOn(framey,Components(nonresPDFPho), LineStyle(kDashed), LineColor(kGreen));
		framex->Draw();
		
    cv = new TCanvas("cv","cv",1600,600);
		cv->Divide(2); 
		cv->cd(1); framex->Draw(); 
		cv->cd(2); framey->Draw(); 
    cv->SaveAs(Form("Plots/AllSignalBkgd/Fits/test_projectionFits_%d.gif",t));
    /*
    cv = new TCanvas("cv","cv",800,600);
    cv->SetLeftMargin(0.15);
    TH2* hcorr = fitResult2D->correlationHist();
    hcorr->Draw("COLZ");
    hcorr->SetStats(0);
    cv->SaveAs(Form("Plots/AllSignalBkgd/Fits/correlation_%d.gif",t));
    
    /*
    cv = new TCanvas("cv","cv",800,600);
    cv->SetLeftMargin(0.15);
    RooPlot* frame = new RooPlot(massBjet, massPho);
    frame->SetTitle("Covariance between massBjet and massPho") ;
    fitResult2D->plotOn(frame, massBjet, massPho, "ME12VHB") ;
    frame->Draw();
    */
  }
}

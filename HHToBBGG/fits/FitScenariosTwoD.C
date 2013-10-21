//
// root -l CMSAna/HHToBBGG/fits/FitScenariosTwoD.C+'("/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_diphotonMass.root", "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_dibjetMass.root", "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_twoDMass.root", "pho", 5000)'
//
// for scanOption: pho, jet, lum
//
//


#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TStyle.h>                 // class to handle ROOT plotting style
#include <TLatex.h>									// class for latex
#include <TFile.h>                  // file handle class
#include <TGraph.h>									// plotting class
#include <TGraphErrors.h>						// plotting error class
#include <TLine.h>									// plotting lines
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
TGraphErrors *gr = 0;
TLine *line = 0;
TLatex *tex = 0;
int alreadyPlotted = 0;

Float_t phoResolution[13] = {0.9/1.53, 1.15/1.53, 1.4/1.53, 1.65/1.53, 1.9/1.53, 2.15/1.53, 2.4/1.53, 2.65/1.53, 2.9/1.53, 3.15/1.53, 3.4/1.53, 3.65/1.53, 3.9/1.53};
Float_t bjetResolution[13] = {10./14.31, 12.5/14.31, 15.0/14.31, 17.5/14.31, 20./14.31, 22.5/14.31, 25./14.31, 27.5/14.31, 30./14.31, 32.5/14.31, 35.0/14.31, 37.5/14.31, 40.0/14.31};
Float_t luminosity[13] = {1.5/3., 2./3., 2.5/3., 3./3., 4./3., 5./3., 6./3., 7./3., 8./3., 9./3., 10./3., 11./3., 12./3.}; 
Float_t fakerate[7] = {0, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03}; 

void AddModels(RooWorkspace *ws, const string inputfilePho, const string inputfileBjet, const string inputfileTwoD, const string scanOption, Int_t s);

void FitScenariosTwoD(const string inputfilePho = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_diphotonMass.root", const string inputfileBjet = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_dibjetMass.root", const string inputfileTwoD = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_twoDMass.root", const string scanOption = "lum", Int_t NToys = 500, Int_t s = 5) {
	
  TRandom3 *randomNumber = new TRandom3(1200);
  RooWorkspace* ws = new RooWorkspace("CMSAna/HHToBBGG/fits/Workspace");
  AddModels(ws, inputfilePho, inputfileBjet, inputfileTwoD, scanOption, s);
  
  //Import variables from workspace
  RooAbsPdf *model2Dpdf = ws->pdf("model2Dpdf");
  RooRealVar *massBjet = ws->var("massBjet");
  RooRealVar *massPho = ws->var("massPho");
  RooRealVar *nsig = ws->var("N (Sig)"); RooRealVar constNsig(*nsig);
  RooRealVar *nres = ws->var("N (ResBkg)"); RooRealVar constNres(*nres);
  RooRealVar *nnonres = ws->var("N (NonResBkg)"); RooRealVar constNnonres(*nnonres);
  
  //Create TTree to store the resulting yield data
  TFile *f = new TFile(Form(("CMSAna/HHToBBGG/data/MassFitResults/ResolutionAnalysis/WithEndcap_LumiTimes2/"+scanOption+"FitScenario%d.root").c_str(), s), "RECREATE");
//   TFile *f = new TFile(Form(("CMSAna/HHToBBGG/data/MassFitResults/ResolutionAnalysis/WithEndcap/"+scanOption+"FitScenario%d.root").c_str(), s), "RECREATE");
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
    cout << t << "|" << s << endl;
    nsig->setVal(constNsig.getVal());
    nres->setVal(constNres.getVal());
    nnonres->setVal(constNnonres.getVal());
    Float_t ran;
    //OptionA
//     if (scanOption == "lum") ran = randomNumber->Poisson(194.82*luminosity[s]);
//     else ran = randomNumber->Poisson(194.82);
//     //NoEndcaps
//     if (scanOption == "lum") ran = randomNumber->Poisson(120.91*luminosity[s]);
//     else ran = randomNumber->Poisson(120.91);
    //NoEndcaps (lumi*2)
    if (scanOption == "lum") ran = randomNumber->Poisson(241.82*luminosity[s]);
    else if (scanOption == "photonFakerate") ran = randomNumber->Poisson(2*(112.18 + 254.66*( fakerate[s] / 0.003 )));
    else ran = randomNumber->Poisson(241.82);   

    RooDataSet *pseudoData2D = model2Dpdf->generate(RooArgList(*massBjet,*massPho), ran);
    RooFitResult *fitResult2D = model2Dpdf->fitTo(*pseudoData2D, RooFit::Save(), RooFit::Extended(kTRUE), RooFit::Strategy(2));
    ws->import(*pseudoData2D, kTRUE);
    ws->import(*pseudoData2D, Rename("pseudoData2D"));
    
    //Save variables into separate branches of root tree
    nsigOut = nsig->getVal();
    nresOut = nres->getVal();
    nnonresOut = nnonres->getVal();
    nsigStd = nsig->getPropagatedError(*fitResult2D);
    nresStd = nres->getPropagatedError(*fitResult2D);
    nnonresStd = nnonres->getPropagatedError(*fitResult2D);
    
    resultTree->Fill();      
  }
  
  //Write to the TTree and close it
  resultTree->Write();
}

void AddModels(RooWorkspace *ws, const string inputfilePho, const string inputfileBjet, const string inputfileTwoD, const string scanOption, Int_t s) {

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
  
  //-------------------------------------------------------------
  //Make Variables for 2D Fits
  //=============================================================
    
  //Variables for diBjet
  RooRealVar massBjet("massBjet","M_{bb}",70.,200.,"GeV/c^{2}");
  RooRealVar massBjetExt("massBjetExt","M_{bb}",40.,200.,"GeV/c^{2}");
  RooRealVar massBjetCut("massBjetCut","M_{bb}",70.,170.,"GeV/c^{2}");
  //Signal diBjet
  RooRealVar sigMeanBjet("sigMeanBjet", "#bar{M}_{bb}", 125., 119., 131., "GeV/c^{2}");
  RooRealVar sigSigmaBjet("sigSigmaBjet", "#sigma_{bb}", 14., 5., 50., "GeV/c^{2}");
  RooRealVar sigAlphaBjet("sigAlphaBjet","#alpha_{bb} (cut)", 1., 0.1, 2.);
  RooRealVar sigPowerBjet("sigPowerBjet", "n_{bb} (power)", 3., 1., 5.);
  //Resonant Background diBjet
  RooRealVar resMeanBjet("resMeanBjet", "#bar{M}_{bb} (ResBkg)", 90., 80., 100., "GeV/c^{2}");
  RooRealVar resSigmaBjet("resSigmaBjet", "#sigma_{bb} (ResBkg)", 10., 2., 50., "GeV/c^{2}");
  RooRealVar resAlphaBjet("resAlphaBjet","#alpha_{bb} (ResBkg) (cut)", -1.1, -2., -0.1);
  RooRealVar resPowerBjet("resPowerBjet", "n_{bb} (ResBkg) (power)", 3., 1., 5.);
  RooRealVar resExpo("resExpo", "#lambda_{bb} (ResBkg)", -.01, -.1, 0.);
  //Non-Resonant Background diBjet
  RooRealVar expRateBjet("expRateBjet", "#lambda_{bb}", -.018, -.1, 0.);
  //Weights for Resonant Background
  RooRealVar nbbH("nbbH", "# bbH events", 3.2, 0,5);
  RooRealVar nOthers("nOthers", "# resonant background events - bbH", 24.2, 0,50);
  
  //Variables for diPhoton
  RooRealVar massPho("massPho","M_{#gamma#gamma}",100.,150.,"GeV/c^{2}");
  RooRealVar massPhoCut("massPhoCut","M_{#gamma#gamma} (Cut)",110.,135.,"GeV/c^{2}");
  //Signal diPhoton
  RooRealVar sigMeanPho("sigMeanPho", "#bar{M}_{#gamma#gamma}", 125., 120., 130., "GeV/c^{2}");
  RooRealVar sigSigmaPho("sigSigmaPho", "#sigma_{#gamma#gamma}", 1.5, 0.5, 4.5, "GeV/c^{2}");
  RooRealVar sigAlphaPho("sigAlphaPho","#alpha_{#gamma#gamma} (cut)", 1., 0.1, 2.);
  RooRealVar sigPowerPho("sigPowerPho", "n_{#gamma#gamma} (power)", 3., 1., 5.);
  //Resonant Background diPhoton
  RooRealVar resMeanPho("resMeanPho", "#bar{M}_{#gamma#gamma} (ResBkg)", 125., 120., 130., "GeV/c^{2}");
  RooRealVar resSigmaPho("resSigmaPho", "#sigma_{#gamma#gamma} (ResBkg)", 1.5, 0.5, 4.5, "GeV/c^{2}");
  RooRealVar resAlphaPho("resAlphaPho","#alpha_{#gamma#gamma} (ResBkg) (cut)", 1., 0.1, 2.);
  RooRealVar resPowerPho("resPowerPho", "n_{#gamma#gamma} (ResBkg) (power)", 3., 1., 5.);
  //Non-Resonant Background diPhoton
  RooRealVar expRatePho("expRatePho", "#lambda_{#gamma#gamma}", -.027, -.1, 0.);
  
  //Weights for the 2D fit
//OptionA
//   RooRealVar nsig("N (Sig)", "# signal events", 5.36, -100,1000);
//   RooRealVar nres("N (ResBkg)", "# resonant background events", 8.86, -50,1000);
//   RooRealVar nnonres("N (NonResBkg)", "# non-resonant background events", 180.6, 0,5000);

// //NoEndcap
//   RooRealVar nsig("N (Sig)", "# signal events", 4.93, -100,1000);
//   RooRealVar nres("N (ResBkg)", "# resonant background events", 7.88, -50,1000);
//   RooRealVar nnonres("N (NonResBkg)", "# non-resonant background events", 108.1, 0,5000);

//NoEndcap (Lumi * 2)
   RooRealVar nsig("N (Sig)", "# signal events", 9.86, -100,1000);
   RooRealVar nres("N (ResBkg)", "# resonant background events", 15.76, -50,1000);
   RooRealVar nnonres("N (NonResBkg)", "# non-resonant background events", 216.2, 0,5000);


  if (scanOption == "pho") {
//     nsig.setVal(nsig.getVal()*phoResolution[s]);
//     nres.setVal(nres.getVal()*phoResolution[s]);
//     nnonres.setVal(nnonres.getVal()*phoResolution[s]);
  }
  else if (scanOption == "jet") {
//     nsig.setVal(nsig.getVal()*bjetResolution[s]);
//     nres.setVal(nres.getVal()*bjetResolution[s]);
//     nnonres.setVal(nnonres.getVal()*bjetResolution[s]);
  } else if (scanOption == "photonFakerate") {
    nsig.setVal( 2*6.11 );
    nres.setVal( 2*10.29 );
    nnonres.setVal( 2*(95.78 + 254.66*( fakerate[s] / 0.003 ))  );
    cout << "Photon Fake rate scan: " << fakerate[s] << " | nnonres = " << nnonres.getVal() << "\n";
  }
  else {
    nsig.setVal(nsig.getVal()*luminosity[s]);
    nres.setVal(nres.getVal()*luminosity[s]);
    nnonres.setVal(nnonres.getVal()*luminosity[s]);
  }
  
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
  //Non-Resonant Background
  RooExponential nonresPDFPho("nonresPDFPho", "Non-Resonant Background PDF #gamma#gamma", massPho, expRatePho);
  RooDataHist nonresPhoData("nonresPhoData", "nonresPhoData", massPho, nonresPhoMass);
  RooFitResult *nonresPhoResult = nonresPDFPho.fitTo(nonresPhoData, RooFit::Strategy(2));
    
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
  //Non-Resonant Background
  RooExponential nonresPDFBjet("nonresPDFBjet", "Non-Resonant Background PDF bb", massBjet, expRateBjet);
  RooDataHist nonresBjetData("nonresBjetData", "nonresBjetData", massBjet, nonresBjetMass);
  RooFitResult *nonresBjetResult = nonresPDFBjet.fitTo(nonresBjetData, RooFit::Strategy(2));
    
  //PDFs for diPhoton*diBjet (2D)
  RooProdPdf sig2Dpdf("sig2Dpdf", "2D Signal PDF", RooArgList(sigPDFPho, sigPDFBjet));
  RooProdPdf res2Dpdf("res2Dpdf", "2D Resonant Background PDF", RooArgList(resPDFPho, resPDFBjet));
  RooProdPdf nonres2Dpdf("nonres2Dpdf", "2D Non-Resonant Background PDF", RooArgList(nonresPDFPho, nonresPDFBjet));
  RooAddPdf model2Dpdf("model2Dpdf", "2D Signal+Background PDF", RooArgList(sig2Dpdf, res2Dpdf, nonres2Dpdf), RooArgList(nsig, nres, nnonres));
  
  //Variables to set constant
  sigMeanBjet.setConstant();  sigMeanPho.setConstant();
  sigSigmaBjet.setConstant(); sigSigmaPho.setConstant();
  sigAlphaBjet.setConstant(); sigAlphaPho.setConstant();
  sigPowerBjet.setConstant(); sigPowerPho.setConstant();
  resMeanBjet.setConstant();  resMeanPho.setConstant();
  resSigmaBjet.setConstant(); resSigmaPho.setConstant();
  resAlphaBjet.setConstant(); resAlphaPho.setConstant();
  resPowerBjet.setConstant(); resPowerPho.setConstant();
  resExpo.setConstant(); nbbH.setConstant(); nOthers.setConstant();
  expRateBjet.setConstant();  expRatePho.setConstant();
  

  //by default use high pileup degraded jet energy resolution
  sigSigmaBjet.setVal(sigSigmaBjet.getVal()*(20./14.31));
  resSigmaBjet.setVal(resSigmaBjet.getVal()*(20./14.31));

  if (scanOption == "jet") {
    cout << "Jet Energy Resolution Scan : " << s << " : " << sigSigmaBjet.getVal() << " * " << bjetResolution[s] << " = " << sigSigmaBjet.getVal()*bjetResolution[s] << "\n";
    sigSigmaBjet.setVal(sigSigmaBjet.getVal()*bjetResolution[s]);
    resSigmaBjet.setVal(resSigmaBjet.getVal()*bjetResolution[s]);
  }
  if (scanOption == "pho") {
    sigSigmaPho.setVal(sigSigmaPho.getVal()*phoResolution[s]);
    resSigmaPho.setVal(resSigmaPho.getVal()*phoResolution[s]);
  }
  
  //Add all 2D stuff to workspace
  ws->import(model2Dpdf);
}



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
Float_t endcapPhotonFakerate[6] = {0/0.003, 0.0001/0.003, 0.0003/0.003, 0.001/0.003, 0.003/0.003, 0.01/0.003}; 

Float_t phoEffRatio[6] = {0.9, 1.0, 1.1, 1.2, 1.3, 1.4};
Float_t btagEffRatio[6] = {0.9, 1.0, 1.1, 1.2, 1.3, 1.4};



//double phoEffRatio = 1.25;
//double btagEffRatio = 1.25;


void AddModels(RooWorkspace *ws, const string inputfilePho, const string inputfileBjet, const string inputfileTwoD, const string scanOption, Int_t s);

void FitScenariosTwoD_RequiredPerformance_card(const string inputfilePho = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_diphotonMass.root", const string inputfileBjet = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_dibjetMass.root", const string inputfileTwoD = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_twoDMass.root", const string scanOption = "lum", Int_t NToys = 500, Int_t s = 5) {
	
  

  TRandom3 *randomNumber = new TRandom3(1200);
  TFile *file = new TFile(Form("HHToBBGGWorkspace_%s_%d.root",scanOption.c_str(),s),"RECREATE");
  RooWorkspace* ws = new RooWorkspace("CMSAna/HHToBBGG/fits/Workspace");
  AddModels(ws, inputfilePho, inputfileBjet, inputfileTwoD, scanOption, s);

  //Import variables from workspace
  RooAbsPdf *model2Dpdf = ws->pdf("model2Dpdf");
  RooAbsPdf *model2Dpdf_cat0 = ws->pdf("model2Dpdf_cat0");
  RooAbsPdf *model2Dpdf_cat1 = ws->pdf("model2Dpdf_cat1");
  RooRealVar *massBjet = ws->var("massBjet");
  RooRealVar *massPho = ws->var("massPho");
  RooCategory *sampleCategory = ws->cat("sampleCategory");

  RooRealVar *nsig = ws->var("nsig"); RooRealVar constNsig(*nsig);
  RooRealVar *nres_cat0 = ws->var("nres_cat0"); RooRealVar constNres_cat0(*nres_cat0);
  RooRealVar *nnonres_cat0 = ws->var("nnonres_cat0"); RooRealVar constNnonres_cat0(*nnonres_cat0);
  RooRealVar *nres_cat1 = ws->var("nres_cat1"); RooRealVar constNres_cat1(*nres_cat1);
  RooRealVar *nnonres_cat1 = ws->var("nnonres_cat1"); RooRealVar constNnonres_cat1(*nnonres_cat1);
  
  //Create TTree to store the resulting yield data
  TFile *f = new TFile(Form(("CMSAna/HHToBBGG/data/MassFitResults/ResolutionAnalysis/OptionB_Categories/tmp/"+scanOption+"FitScenario%d.root").c_str(), s), "RECREATE");
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
  //Yield Information
  //=============================================================     
  double myPhoEffRatio = 1.00;
  double myBtagEffRatio = 1.00;
  double myEndcapPhotonFakerateRatio = 1.0;
  if (scanOption == "phoEff") {
    myPhoEffRatio = phoEffRatio[s];
    myBtagEffRatio = 1.25;
  }
  if (scanOption == "btagEff") {
    myBtagEffRatio = btagEffRatio[s];
    myPhoEffRatio = 1.25;
  }
  if (scanOption == "endcapPhotonFakerate") myEndcapPhotonFakerateRatio = endcapPhotonFakerate[s];
  
  double totalYield = 
    6.11*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio +
    7.88*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio +      
    31.1*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio
    + 26.4*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio
    + 48.7*myBtagEffRatio*myBtagEffRatio*myPhoEffRatio
    + 1.8*myBtagEffRatio*myBtagEffRatio
    +
    2.41*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio +
    17.4*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio
    + 17.2*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio
    + 206.0*myBtagEffRatio*myBtagEffRatio*myPhoEffRatio*myEndcapPhotonFakerateRatio
    + 1.7*myBtagEffRatio*myBtagEffRatio
    ;
  double cat0Yield = 6.11*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio +
    7.88*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio +      
    31.1*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio
    + 26.4*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio
    + 48.7*myBtagEffRatio*myBtagEffRatio*myPhoEffRatio
    + 1.8*myBtagEffRatio*myBtagEffRatio;
  double cat1Yield = 2.41*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio +
    17.4*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio
    + 17.2*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio
    + 206.0*myBtagEffRatio*myBtagEffRatio*myPhoEffRatio*myEndcapPhotonFakerateRatio
    + 1.7*myBtagEffRatio*myBtagEffRatio;
  double cat0Fraction = cat0Yield / totalYield;
  double cat1Fraction = cat1Yield / totalYield;

  //-------------------------------------------------------------
  //Generate Toy MC experiment data and fit
  //=============================================================     
  for(UInt_t t=0; t < NToys; ++t) {
    cout << t << "|" << s << endl;
    nsig->setVal(constNsig.getVal());
    nres_cat0->setVal(constNres_cat0.getVal());
    nnonres_cat0->setVal(constNnonres_cat0.getVal());
    nres_cat1->setVal(constNres_cat1.getVal());
    nnonres_cat1->setVal(constNnonres_cat1.getVal());
    Float_t ran;

    //if (scanOption == "lum") ran = randomNumber->Poisson(totalYield*luminosity[s]);
    //else ran = randomNumber->Poisson(totalYield);

    if (scanOption == "lum") ran = totalYield*luminosity[s];
    else ran = totalYield;


    RooDataSet *pseudoData2D = model2Dpdf->generate(RooArgList(*massBjet,*massPho, *sampleCategory), ran);
    RooFitResult *fitResult2D = model2Dpdf->fitTo(*pseudoData2D, RooFit::Save(), RooFit::Extended(kTRUE), RooFit::Strategy(2));
    ws->import(*pseudoData2D, kTRUE);
    ws->import(*pseudoData2D, Rename("pseudoData2D"));

    RooDataSet *pseudoData2D_cat0 = model2Dpdf_cat0->generate(RooArgList(*massBjet,*massPho), ran*cat0Fraction);
    RooDataSet *pseudoData2D_cat1 = model2Dpdf_cat1->generate(RooArgList(*massBjet,*massPho), ran*cat1Fraction);
    ws->import(*pseudoData2D_cat0, kTRUE);
    ws->import(*pseudoData2D_cat0, Rename("pseudoData2D_cat0"));
    ws->import(*pseudoData2D_cat1, kTRUE);
    ws->import(*pseudoData2D_cat1, Rename("pseudoData2D_cat1"));
    
    //Save variables into separate branches of root tree
    nsigOut = nsig->getVal();
    nresOut = nres_cat0->getVal();
    nnonresOut = nnonres_cat0->getVal();
    nsigStd = nsig->getPropagatedError(*fitResult2D);
    nresStd = nres_cat0->getPropagatedError(*fitResult2D);
    nnonresStd = nnonres_cat0->getPropagatedError(*fitResult2D);
    
    resultTree->Fill();      
  }
  
  //Write to the TTree and close it
  resultTree->Write();

  file->WriteTObject(ws, Form("HHToBBGGWorkspace",scanOption.c_str(),s), "WriteDelete");
  file->Close();
  delete file;




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
  //category 0: barrel-barrel photons
  //category 1: barrel-endcap & endcap-endcap photons

  //********************** Category0 *******************************
  //Variables for diBjet
  RooRealVar massBjet("massBjet","M_{bb}",70.,200.,"GeV/c^{2}");
  RooRealVar massBjetExt_cat0("massBjetExt_cat0","M_{bb}",40.,200.,"GeV/c^{2}");
  RooRealVar massBjetCut_cat0("massBjetCut_cat0","M_{bb}",70.,170.,"GeV/c^{2}");
  //Signal diBjet
  RooRealVar sigMeanBjet_cat0("sigMeanBjet_cat0", "#bar{M}_{bb}", 125., 119., 131., "GeV/c^{2}");
  RooRealVar sigSigmaBjet_cat0("sigSigmaBjet_cat0", "#sigma_{bb}", 14., 5., 50., "GeV/c^{2}");
  RooRealVar sigAlphaBjet_cat0("sigAlphaBjet_cat0","#alpha_{bb} (cut)", 1., 0.1, 2.);
  RooRealVar sigPowerBjet_cat0("sigPowerBjet_cat0", "n_{bb} (power)", 3., 1., 5.);
  //Resonant Background diBjet
  RooRealVar resMeanBjet_cat0("resMeanBjet_cat0", "#bar{M}_{bb} (ResBkg)", 90., 80., 100., "GeV/c^{2}");
  RooRealVar resSigmaBjet_cat0("resSigmaBjet_cat0", "#sigma_{bb} (ResBkg)", 10., 2., 50., "GeV/c^{2}");
  RooRealVar resAlphaBjet_cat0("resAlphaBjet_cat0","#alpha_{bb} (ResBkg) (cut)", -1.1, -2., -0.1);
  RooRealVar resPowerBjet_cat0("resPowerBjet_cat0", "n_{bb} (ResBkg) (power)", 3., 1., 5.);
  RooRealVar resExpo_cat0("resExpo_cat0", "#lambda_{bb} (ResBkg)", -.01, -.1, 0.);
  //Non-Resonant Background diBjet
  RooRealVar expRateBjet_cat0("expRateBjet_cat0", "#lambda_{bb}", -.018, -.1, 0.);
  //Weights for Resonant Background
  RooRealVar nbbH_cat0("nbbH_cat0", "# bbH events", 3.2, 0,5);
  RooRealVar nOthers_cat0("nOthers_cat0", "# resonant background events - bbH", 24.2, 0,50);
  
  //Variables for diPhoton
  RooRealVar massPho("massPho","M_{#gamma#gamma}",100.,150.,"GeV/c^{2}");
  RooRealVar massPhoCut_cat0("massPhoCut_cat0","M_{#gamma#gamma} (Cut)",110.,135.,"GeV/c^{2}");
  //Signal diPhoton
  RooRealVar sigMeanPho_cat0("sigMeanPho_cat0", "#bar{M}_{#gamma#gamma}", 125., 120., 130., "GeV/c^{2}");
  RooRealVar sigSigmaPho_cat0("sigSigmaPho_cat0", "#sigma_{#gamma#gamma}", 1.5, 0.5, 4.5, "GeV/c^{2}");
  RooRealVar sigAlphaPho_cat0("sigAlphaPho_cat0","#alpha_{#gamma#gamma} (cut)", 1., 0.1, 2.);
  RooRealVar sigPowerPho_cat0("sigPowerPho_cat0", "n_{#gamma#gamma} (power)", 3., 1., 5.);
  //Resonant Background diPhoton
  RooRealVar resMeanPho_cat0("resMeanPho_cat0", "#bar{M}_{#gamma#gamma} (ResBkg)", 125., 120., 130., "GeV/c^{2}");
  RooRealVar resSigmaPho_cat0("resSigmaPho_cat0", "#sigma_{#gamma#gamma} (ResBkg)", 1.5, 0.5, 4.5, "GeV/c^{2}");
  RooRealVar resAlphaPho_cat0("resAlphaPho_cat0","#alpha_{#gamma#gamma} (ResBkg) (cut)", 1., 0.1, 2.);
  RooRealVar resPowerPho_cat0("resPowerPho_cat0", "n_{#gamma#gamma} (ResBkg) (power)", 3., 1., 5.);
  //Non-Resonant Background diPhoton
  RooRealVar expRatePho_cat0("expRatePho_cat0", "#lambda_{#gamma#gamma}", -.027, -.1, 0.);
  

  //********************** Category1 *******************************
  //Variables for diBjet
  RooRealVar massBjetExt_cat1("massBjetExt_cat1","M_{bb}",40.,200.,"GeV/c^{2}");
  RooRealVar massBjetCut_cat1("massBjetCut_cat1","M_{bb}",70.,170.,"GeV/c^{2}");
  //Signal diBjet
  RooRealVar sigMeanBjet_cat1("sigMeanBjet_cat1", "#bar{M}_{bb}", 125., 119., 131., "GeV/c^{2}");
  RooRealVar sigSigmaBjet_cat1("sigSigmaBjet_cat1", "#sigma_{bb}", 14., 5., 50., "GeV/c^{2}");
  RooRealVar sigAlphaBjet_cat1("sigAlphaBjet_cat1","#alpha_{bb} (cut)", 1., 0.1, 2.);
  RooRealVar sigPowerBjet_cat1("sigPowerBjet_cat1", "n_{bb} (power)", 3., 1., 5.);
  //Resonant Background diBjet
  RooRealVar resMeanBjet_cat1("resMeanBjet_cat1", "#bar{M}_{bb} (ResBkg)", 90., 80., 100., "GeV/c^{2}");
  RooRealVar resSigmaBjet_cat1("resSigmaBjet_cat1", "#sigma_{bb} (ResBkg)", 10., 2., 50., "GeV/c^{2}");
  RooRealVar resAlphaBjet_cat1("resAlphaBjet_cat1","#alpha_{bb} (ResBkg) (cut)", -1.1, -2., -0.1);
  RooRealVar resPowerBjet_cat1("resPowerBjet_cat1", "n_{bb} (ResBkg) (power)", 3., 1., 5.);
  RooRealVar resExpo_cat1("resExpo_cat1", "#lambda_{bb} (ResBkg)", -.01, -.1, 0.);
  //Non-Resonant Background diBjet
  RooRealVar expRateBjet_cat1("expRateBjet_cat1", "#lambda_{bb}", -.018, -.1, 0.);
  //Weights for Resonant Background
  RooRealVar nbbH_cat1("nbbH_cat1", "# bbH events", 3.2, 0,5);
  RooRealVar nOthers_cat1("nOthers_cat1", "# resonant background events - bbH", 24.2, 0,50);
  
  //Variables for diPhoton
  RooRealVar massPhoCut_cat1("massPhoCut_cat1","M_{#gamma#gamma} (Cut)",110.,135.,"GeV/c^{2}");
  //Signal diPhoton
  RooRealVar sigMeanPho_cat1("sigMeanPho_cat1", "#bar{M}_{#gamma#gamma}", 125., 120., 130., "GeV/c^{2}");
  RooRealVar sigSigmaPho_cat1("sigSigmaPho_cat1", "#sigma_{#gamma#gamma}", 1.5, 0.5, 4.5, "GeV/c^{2}");
  RooRealVar sigAlphaPho_cat1("sigAlphaPho_cat1","#alpha_{#gamma#gamma} (cut)", 1., 0.1, 2.);
  RooRealVar sigPowerPho_cat1("sigPowerPho_cat1", "n_{#gamma#gamma} (power)", 3., 1., 5.);
  //Resonant Background diPhoton
  RooRealVar resMeanPho_cat1("resMeanPho_cat1", "#bar{M}_{#gamma#gamma} (ResBkg)", 125., 120., 130., "GeV/c^{2}");
  RooRealVar resSigmaPho_cat1("resSigmaPho_cat1", "#sigma_{#gamma#gamma} (ResBkg)", 1.5, 0.5, 4.5, "GeV/c^{2}");
  RooRealVar resAlphaPho_cat1("resAlphaPho_cat1","#alpha_{#gamma#gamma} (ResBkg) (cut)", 1., 0.1, 2.);
  RooRealVar resPowerPho_cat1("resPowerPho_cat1", "n_{#gamma#gamma} (ResBkg) (power)", 3., 1., 5.);
  //Non-Resonant Background diPhoton
  RooRealVar expRatePho_cat1("expRatePho_cat1", "#lambda_{#gamma#gamma}", -.027, -.1, 0.);


  //****************************
  //yields
  //****************************
  double myPhoEffRatio = 1.00;
  double myBtagEffRatio = 1.00;
  double myEndcapPhotonFakerateRatio = 1.0;
  if (scanOption == "phoEff") {
    myPhoEffRatio = phoEffRatio[s];
    myBtagEffRatio = 1.25;
  }
  if (scanOption == "btagEff") {
    myBtagEffRatio = btagEffRatio[s];
    myPhoEffRatio = 1.25;
  }
  if (scanOption == "endcapPhotonFakerate") myEndcapPhotonFakerateRatio = endcapPhotonFakerate[s];
  

  //With Endcap photons, but assuming factor of 10 reduction in endcap fake photons
  RooRealVar nsig("nsig", "# signal events", 6.11*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio, -100,1000); //0.807*6.11 in cat0, 0.193*6.11 in cat1
  RooFormulaVar nsig_cat0("nsig_cat0","nsig*0.807",RooArgList(nsig));
  RooFormulaVar nsig_cat1("nsig_cat1","nsig*0.193",RooArgList(nsig));
  RooRealVar nres_cat0("nres_cat0", "# resonant background events cat0", 7.88*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio, -50,1000);
  RooRealVar nnonres_cat0("nnonres_cat0", "# non-resonant background events cat0", 
                          31.1*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio
                          + 26.4*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio
                          + 48.7*myBtagEffRatio*myBtagEffRatio*myPhoEffRatio
                          + 1.8*myBtagEffRatio*myBtagEffRatio
                          , 0,5000);  
  RooRealVar nres_cat1("nres_cat1", "# resonant background events cat1", 2.41*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio, -50,1000);
  RooRealVar nnonres_cat1("nnonres_cat1", "# non-resonant background events cat1", 
                          17.4*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio
                          + 17.2*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio
                          + 206.0*myBtagEffRatio*myBtagEffRatio*myPhoEffRatio*myEndcapPhotonFakerateRatio
                          + 1.7*myBtagEffRatio*myBtagEffRatio
                          , 0,5000);
  
  double lumiRatio = 1.0;

  if (scanOption == "pho") {
  }
  else if (scanOption == "jet") {
  }
  else if (scanOption == "lum") {
    lumiRatio = luminosity[s];
    nsig.setVal(nsig.getVal()*luminosity[s]);
    nres_cat0.setVal(nres_cat0.getVal()*luminosity[s]);
    nnonres_cat0.setVal(nnonres_cat0.getVal()*luminosity[s]);
    nres_cat1.setVal(nres_cat1.getVal()*luminosity[s]);
    nnonres_cat1.setVal(nnonres_cat1.getVal()*luminosity[s]);
  } else {
    nsig.setVal(nsig.getVal());
    nres_cat0.setVal(nres_cat0.getVal());
    nnonres_cat0.setVal(nnonres_cat0.getVal());
    nres_cat1.setVal(nres_cat1.getVal());
    nnonres_cat1.setVal(nnonres_cat1.getVal());
  }
  
  //-------------------------------------------------------------
  //Make PDFs for Signal + Background
  //=============================================================
  
    
  //PDFs for diPhoton
  //Preliminary Fit for Cut Signal
  RooCBShape sigPDFPhoCut_cat0("sigPDFPhoCut_cat0", "Signal PDF #gamma#gamma", massPhoCut_cat0, sigMeanPho_cat0, sigSigmaPho_cat0, sigAlphaPho_cat0, sigPowerPho_cat0);
  RooDataHist sigPhoData_cat0("sigPhoData_cat0", "sigPhoData cat0", massPhoCut_cat0, sigPhoMass);
  RooFitResult *sigPhoCBResult_cat0 = sigPDFPhoCut_cat0.fitTo(sigPhoData_cat0, RooFit::Strategy(2));
  //Signal
  RooCBShape sigPDFPho_cat0("sigPDFPho_cat0", "Signal PDF #gamma#gamma cat0", massPho, sigMeanPho_cat0, sigSigmaPho_cat0, sigAlphaPho_cat0, sigPowerPho_cat0);
  //Preliminary Fit for Cut Resonant Background 
  RooCBShape resPDFPhoCut_cat0("resPDFPhoCut", "Resonant Background Cut PDF #gamma#gamma cat0", massPhoCut_cat0, resMeanPho_cat0, resSigmaPho_cat0, resAlphaPho_cat0, resPowerPho_cat0);
  RooDataHist resPhoData_cat0("resPhoData_cat0", "resPhoData_cat0", massPhoCut_cat0, resPhoMass);
  RooFitResult *resPhoCBResult_cat0 = resPDFPhoCut_cat0.fitTo(resPhoData_cat0, RooFit::Strategy(2));
  //Resonant Background
  RooCBShape resPDFPho_cat0("resPDFPho_cat0", "Resonant Background PDF #gamma#gamma", massPho, resMeanPho_cat0, resSigmaPho_cat0, resAlphaPho_cat0, resPowerPho_cat0);
  //Non-Resonant Background
  RooExponential nonresPDFPho_cat0("nonresPDFPho_cat0", "Non-Resonant Background PDF #gamma#gamma cat0", massPho, expRatePho_cat0);
  RooDataHist nonresPhoData_cat0("nonresPhoData_cat0", "nonresPhoData_cat0", massPho, nonresPhoMass);
  RooFitResult *nonresPhoResult_cat0 = nonresPDFPho_cat0.fitTo(nonresPhoData_cat0, RooFit::Strategy(2));
    
  //PDFs for diBjet
  //Preliminary Fit for Cut Signal
  RooCBShape sigPDFBjetCut_cat0("sigPDFBjetCut_cat0", "Signal PDF bb cat0", massBjetCut_cat0, sigMeanBjet_cat0, sigSigmaBjet_cat0, sigAlphaBjet_cat0, sigPowerBjet_cat0);
  RooDataHist sigBjetData_cat0("sigBjetData_cat0", "sigBjetData_cat0", massBjetCut_cat0, sigBjetMass);
  RooFitResult *sigBjetCBResult_cat0 = sigPDFBjetCut_cat0.fitTo(sigBjetData_cat0, RooFit::Strategy(2));
  //Signal
  RooCBShape sigPDFBjet_cat0("sigPDFBjet_cat0", "Signal PDF bb cat0", massBjet, sigMeanBjet_cat0, sigSigmaBjet_cat0, sigAlphaBjet_cat0, sigPowerBjet_cat0);
  //Preliminary Fit for Extended Resonant Background
  RooCBShape resPDFBjetCBExt_cat0("resPDFBjetCBExt_cat0", "Resonant CB Background PDF bb cat0", massBjetExt_cat0, resMeanBjet_cat0, resSigmaBjet_cat0, resAlphaBjet_cat0, resPowerBjet_cat0);
  RooExponential resPDFExpoExt_cat0("resPDFExpoExt_cat0", "Resonant Exponential Background PDF bb cat0", massBjetExt_cat0, resExpo_cat0);
  RooAddPdf resPDFBjetExt_cat0("resPDFBjetExt_cat0", "Extended Resonant Background PDF bb cat0", RooArgList(resPDFBjetCBExt_cat0, resPDFExpoExt_cat0), RooArgList(nOthers_cat0, nbbH_cat0));
  RooDataHist resBjetDataExt_cat0("resBjetDataExt_cat0", "resBjetDataExt_cat0", massBjetExt_cat0, resBjetMassExt);
  RooFitResult *resBjetCBResult_cat0 = resPDFBjetExt_cat0.fitTo(resBjetDataExt_cat0, RooFit::Strategy(2));
  //Resonant Background
  RooCBShape resPDFBjetCB_cat0("resPDFBjetCB_cat0", "Resonant CB Background PDF bb cat0", massBjet, resMeanBjet_cat0, resSigmaBjet_cat0, resAlphaBjet_cat0, resPowerBjet_cat0);
  RooExponential resPDFExpo_cat0("resPDFExp_cat0", "Resonant Exponential Background PDF bb cat0", massBjet, resExpo_cat0);
  RooAddPdf resPDFBjet_cat0("resPDFBjet_cat0", "Resonant Background PDF bb cat0", RooArgList(resPDFBjetCB_cat0, resPDFExpo_cat0), RooArgList(nOthers_cat0, nbbH_cat0));
  //Non-Resonant Background
  RooExponential nonresPDFBjet_cat0("nonresPDFBjet_cat0", "Non-Resonant Background PDF bb cat0", massBjet, expRateBjet_cat0);
  RooDataHist nonresBjetData_cat0("nonresBjetData_cat0", "nonresBjetData_cat0", massBjet, nonresBjetMass);
  RooFitResult *nonresBjetResult_cat0 = nonresPDFBjet_cat0.fitTo(nonresBjetData_cat0, RooFit::Strategy(2));
    

  //Variables to set constant
  sigMeanBjet_cat0.setConstant();  sigMeanPho_cat0.setConstant();
  sigSigmaBjet_cat0.setConstant(); sigSigmaPho_cat0.setConstant();
  sigAlphaBjet_cat0.setConstant(); sigAlphaPho_cat0.setConstant();
  sigPowerBjet_cat0.setConstant(); sigPowerPho_cat0.setConstant();
  resMeanBjet_cat0.setConstant();  resMeanPho_cat0.setConstant();
  resSigmaBjet_cat0.setConstant(); resSigmaPho_cat0.setConstant();
  resAlphaBjet_cat0.setConstant(); resAlphaPho_cat0.setConstant();
  resPowerBjet_cat0.setConstant(); resPowerPho_cat0.setConstant();
  resExpo_cat0.setConstant(); nbbH_cat0.setConstant(); nOthers_cat0.setConstant();
  expRateBjet_cat0.setConstant();  expRatePho_cat0.setConstant();
  


    
  //PDFs for diPhoton
  //Preliminary Fit for Cut Signal
  RooCBShape sigPDFPhoCut_cat1("sigPDFPhoCut_cat1", "Signal PDF #gamma#gamma", massPhoCut_cat1, sigMeanPho_cat1, sigSigmaPho_cat1, sigAlphaPho_cat1, sigPowerPho_cat1);
  RooDataHist sigPhoData_cat1("sigPhoData_cat1", "sigPhoData cat1", massPhoCut_cat1, sigPhoMass);
  RooFitResult *sigPhoCBResult_cat1 = sigPDFPhoCut_cat1.fitTo(sigPhoData_cat1, RooFit::Strategy(2));
  //Signal
  RooCBShape sigPDFPho_cat1("sigPDFPho_cat1", "Signal PDF #gamma#gamma cat1", massPho, sigMeanPho_cat1, sigSigmaPho_cat1, sigAlphaPho_cat1, sigPowerPho_cat1);
  //Preliminary Fit for Cut Resonant Background 
  RooCBShape resPDFPhoCut_cat1("resPDFPhoCut", "Resonant Background Cut PDF #gamma#gamma cat1", massPhoCut_cat1, resMeanPho_cat1, resSigmaPho_cat1, resAlphaPho_cat1, resPowerPho_cat1);
  RooDataHist resPhoData_cat1("resPhoData_cat1", "resPhoData_cat1", massPhoCut_cat1, resPhoMass);
  RooFitResult *resPhoCBResult_cat1 = resPDFPhoCut_cat1.fitTo(resPhoData_cat1, RooFit::Strategy(2));
  //Resonant Background
  RooCBShape resPDFPho_cat1("resPDFPho_cat1", "Resonant Background PDF #gamma#gamma", massPho, resMeanPho_cat1, resSigmaPho_cat1, resAlphaPho_cat1, resPowerPho_cat1);
  //Non-Resonant Background
  RooExponential nonresPDFPho_cat1("nonresPDFPho_cat1", "Non-Resonant Background PDF #gamma#gamma cat1", massPho, expRatePho_cat1);
  RooDataHist nonresPhoData_cat1("nonresPhoData_cat1", "nonresPhoData_cat1", massPho, nonresPhoMass);
  RooFitResult *nonresPhoResult_cat1 = nonresPDFPho_cat1.fitTo(nonresPhoData_cat1, RooFit::Strategy(2));
    
  //PDFs for diBjet
  //Preliminary Fit for Cut Signal
  RooCBShape sigPDFBjetCut_cat1("sigPDFBjetCut_cat1", "Signal PDF bb cat1", massBjetCut_cat1, sigMeanBjet_cat1, sigSigmaBjet_cat1, sigAlphaBjet_cat1, sigPowerBjet_cat1);
  RooDataHist sigBjetData_cat1("sigBjetData_cat1", "sigBjetData_cat1", massBjetCut_cat1, sigBjetMass);
  RooFitResult *sigBjetCBResult_cat1 = sigPDFBjetCut_cat1.fitTo(sigBjetData_cat1, RooFit::Strategy(2));
  //Signal
  RooCBShape sigPDFBjet_cat1("sigPDFBjet_cat1", "Signal PDF bb cat1", massBjet, sigMeanBjet_cat1, sigSigmaBjet_cat1, sigAlphaBjet_cat1, sigPowerBjet_cat1);
  //Preliminary Fit for Extended Resonant Background
  RooCBShape resPDFBjetCBExt_cat1("resPDFBjetCBExt_cat1", "Resonant CB Background PDF bb cat1", massBjetExt_cat1, resMeanBjet_cat1, resSigmaBjet_cat1, resAlphaBjet_cat1, resPowerBjet_cat1);
  RooExponential resPDFExpoExt_cat1("resPDFExpoExt_cat1", "Resonant Exponential Background PDF bb cat1", massBjetExt_cat1, resExpo_cat1);
  RooAddPdf resPDFBjetExt_cat1("resPDFBjetExt_cat1", "Extended Resonant Background PDF bb cat1", RooArgList(resPDFBjetCBExt_cat1, resPDFExpoExt_cat1), RooArgList(nOthers_cat1, nbbH_cat1));
  RooDataHist resBjetDataExt_cat1("resBjetDataExt_cat1", "resBjetDataExt_cat1", massBjetExt_cat1, resBjetMassExt);
  RooFitResult *resBjetCBResult_cat1 = resPDFBjetExt_cat1.fitTo(resBjetDataExt_cat1, RooFit::Strategy(2));
  //Resonant Background
  RooCBShape resPDFBjetCB_cat1("resPDFBjetCB_cat1", "Resonant CB Background PDF bb cat1", massBjet, resMeanBjet_cat1, resSigmaBjet_cat1, resAlphaBjet_cat1, resPowerBjet_cat1);
  RooExponential resPDFExpo_cat1("resPDFExp_cat1", "Resonant Exponential Background PDF bb cat1", massBjet, resExpo_cat1);
  RooAddPdf resPDFBjet_cat1("resPDFBjet_cat1", "Resonant Background PDF bb cat1", RooArgList(resPDFBjetCB_cat1, resPDFExpo_cat1), RooArgList(nOthers_cat1, nbbH_cat1));
  //Non-Resonant Background
  RooExponential nonresPDFBjet_cat1("nonresPDFBjet_cat1", "Non-Resonant Background PDF bb cat1", massBjet, expRateBjet_cat1);
  RooDataHist nonresBjetData_cat1("nonresBjetData_cat1", "nonresBjetData_cat1", massBjet, nonresBjetMass);
  RooFitResult *nonresBjetResult_cat1 = nonresPDFBjet_cat1.fitTo(nonresBjetData_cat1, RooFit::Strategy(2));
    

  //Variables to set constant
  sigMeanBjet_cat1.setConstant();  sigMeanPho_cat1.setConstant();
  sigSigmaBjet_cat1.setConstant(); sigSigmaPho_cat1.setConstant();
  sigAlphaBjet_cat1.setConstant(); sigAlphaPho_cat1.setConstant();
  sigPowerBjet_cat1.setConstant(); sigPowerPho_cat1.setConstant();
  resMeanBjet_cat1.setConstant();  resMeanPho_cat1.setConstant();
  resSigmaBjet_cat1.setConstant(); resSigmaPho_cat1.setConstant();
  resAlphaBjet_cat1.setConstant(); resAlphaPho_cat1.setConstant();
  resPowerBjet_cat1.setConstant(); resPowerPho_cat1.setConstant();
  resExpo_cat1.setConstant(); nbbH_cat1.setConstant(); nOthers_cat1.setConstant();
  expRateBjet_cat1.setConstant();  expRatePho_cat1.setConstant();
  

  //by default use degraded jet energy resolution
  if (scanOption != "phoEff" && scanOption != "btagEff" && scanOption != "jet") {
     sigSigmaBjet_cat0.setVal(sigSigmaBjet_cat0.getVal()*(20./14.31));
     resSigmaBjet_cat0.setVal(resSigmaBjet_cat0.getVal()*(20./14.31));
     sigSigmaBjet_cat1.setVal(sigSigmaBjet_cat1.getVal()*(20./14.31));
     resSigmaBjet_cat1.setVal(resSigmaBjet_cat1.getVal()*(20./14.31));
  }
  
  if (scanOption == "lum") {
  }

  if (scanOption == "jet") {
    cout << "Jet Energy Resolution Scan cat0: " << s << " : " << sigSigmaBjet_cat0.getVal() << " * " << bjetResolution[s] << " = " << sigSigmaBjet_cat0.getVal()*bjetResolution[s] << "\n";
    sigSigmaBjet_cat0.setVal(sigSigmaBjet_cat0.getVal()*bjetResolution[s]);
    resSigmaBjet_cat0.setVal(resSigmaBjet_cat0.getVal()*bjetResolution[s]);
    sigSigmaBjet_cat1.setVal(sigSigmaBjet_cat1.getVal()*bjetResolution[s]);
    resSigmaBjet_cat1.setVal(resSigmaBjet_cat1.getVal()*bjetResolution[s]);
  }
  if (scanOption == "pho") {
    cout << "Pho Energy Resolution Scan cat0: " << s << " : " << sigSigmaPho_cat0.getVal() << " * " << phoResolution[s] << " = " << sigSigmaPho_cat0.getVal()*phoResolution[s] << "\n";
    sigSigmaPho_cat0.setVal(sigSigmaPho_cat0.getVal()*phoResolution[s]);
    resSigmaPho_cat0.setVal(resSigmaPho_cat0.getVal()*phoResolution[s]);
    sigSigmaPho_cat1.setVal(sigSigmaPho_cat1.getVal()*phoResolution[s]);
    resSigmaPho_cat1.setVal(resSigmaPho_cat1.getVal()*phoResolution[s]);
  }

  RooCategory sampleCategory("sampleCategory","sampleCategory") ;
  sampleCategory.defineType("cat0") ;
  sampleCategory.defineType("cat1") ;

  //PDFs for diPhoton*diBjet (2D)
  RooProdPdf sig2Dpdf_cat0("sig2Dpdf_cat0", "2D Signal PDF cat0", RooArgList(sigPDFPho_cat0, sigPDFBjet_cat0));
  RooProdPdf res2Dpdf_cat0("res2Dpdf_cat0", "2D Resonant Background PDF cat0", RooArgList(resPDFPho_cat0, resPDFBjet_cat0));
  RooProdPdf nonres2Dpdf_cat0("nonres2Dpdf_cat0", "2D Non-Resonant Background PDF cat0", RooArgList(nonresPDFPho_cat0, nonresPDFBjet_cat0));
  RooAddPdf model2Dpdf_cat0("model2Dpdf_cat0", "2D Signal+Background PDF cat0", RooArgList(sig2Dpdf_cat0, res2Dpdf_cat0, nonres2Dpdf_cat0), RooArgList(nsig_cat0, nres_cat0, nnonres_cat0));

  RooProdPdf sig2Dpdf_cat1("sig2Dpdf_cat1", "2D Signal PDF cat1", RooArgList(sigPDFPho_cat1, sigPDFBjet_cat1));
  RooProdPdf res2Dpdf_cat1("res2Dpdf_cat1", "2D Resonant Background PDF cat1", RooArgList(resPDFPho_cat1, resPDFBjet_cat1));
  RooProdPdf nonres2Dpdf_cat1("nonres2Dpdf_cat1", "2D Non-Resonant Background PDF cat1", RooArgList(nonresPDFPho_cat1, nonresPDFBjet_cat1));
  RooAddPdf model2Dpdf_cat1("model2Dpdf_cat1", "2D Signal+Background PDF cat1", RooArgList(sig2Dpdf_cat1, res2Dpdf_cat1, nonres2Dpdf_cat1), RooArgList(nsig_cat1, nres_cat1, nnonres_cat1));
  
  RooSimultaneous model2Dpdf("model2Dpdf","model2Dpdf",sampleCategory) ;   
  model2Dpdf.addPdf(model2Dpdf_cat0,"cat0") ;
  model2Dpdf.addPdf(model2Dpdf_cat1,"cat1") ;
  
  //Add all 2D stuff to workspace
  ws->import(model2Dpdf);


  //************************************************
  //Make Datacard
  //************************************************
  ofstream datacard;
  datacard.open(Form("HHToBBGG_card_%s_%d.txt",scanOption.c_str(),s));
  datacard << "imax 2 number of bins" << endl;
  datacard << "jmax 2 number of processes minus 1" << endl;
  datacard << "kmax 6 number of nuisance parameters" << endl;
  datacard << "----------------------------------------------------------------------------------------------------------------------------------" << endl;
  datacard << Form("shapes bkg_nonres ch0       HHToBBGGWorkspace_%s_%d.root HHToBBGGWorkspace:nonres2Dpdf_cat0",scanOption.c_str(),s) << endl;
  datacard << Form("shapes bkg_res    ch0       HHToBBGGWorkspace_%s_%d.root HHToBBGGWorkspace:res2Dpdf_cat0",scanOption.c_str(),s) << endl;
  datacard << Form("shapes hhbbgg     ch0       HHToBBGGWorkspace_%s_%d.root HHToBBGGWorkspace:sig2Dpdf_cat0",scanOption.c_str(),s) << endl;
  datacard << Form("shapes data_obs   ch0       HHToBBGGWorkspace_%s_%d.root HHToBBGGWorkspace:pseudoData2D_cat0",scanOption.c_str(),s) << endl;
  datacard << Form("shapes bkg_nonres ch1       HHToBBGGWorkspace_%s_%d.root HHToBBGGWorkspace:nonres2Dpdf_cat1",scanOption.c_str(),s) << endl;
  datacard << Form("shapes bkg_res    ch1       HHToBBGGWorkspace_%s_%d.root HHToBBGGWorkspace:res2Dpdf_cat1",scanOption.c_str(),s) << endl;
  datacard << Form("shapes hhbbgg     ch1       HHToBBGGWorkspace_%s_%d.root HHToBBGGWorkspace:sig2Dpdf_cat1",scanOption.c_str(),s) << endl;
  datacard << Form("shapes data_obs   ch1       HHToBBGGWorkspace_%s_%d.root HHToBBGGWorkspace:pseudoData2D_cat1",scanOption.c_str(),s) << endl;
  datacard << "" << endl;
  datacard << "----------------------------------------------------------------------------------------------------------------------------------" << endl;
  datacard << "bin          ch0          ch1" << endl;
  datacard << "observation  -1.0         -1.0" << endl;
  datacard << "----------------------------------------------------------------------------------------------------------------------------------" << endl;
  datacard <<      "bin                                       ch0          ch0          ch0          ch1          ch1          ch1        " << endl;
  datacard <<      "process                                   hhbbgg       bkg_res      bkg_nonres   hhbbgg       bkg_res      bkg_nonres " << endl;
  datacard <<      "process                                   0            1            2            0            1            2          " << endl;
  datacard << Form("rate                                      %.2f         %.2f         %.2f         %.2f         %.2f         %.2f       ",
                   0.807*nsig.getVal()*lumiRatio,
                   nres_cat0.getVal()*lumiRatio,
                   nnonres_cat0.getVal()*lumiRatio,
                   0.193*nsig.getVal()*lumiRatio,
                   nres_cat1.getVal()*lumiRatio,
                   nnonres_cat1.getVal()*lumiRatio)
           << endl;
  datacard <<      "----------------------------------------------------------------------------------------------------------------------------------" << endl;
  datacard <<      "CMS_eff_g_7TeV          lnN               1.04         1.04         -            1.04         1.04         -          " << endl;
  datacard <<      "CMS_eff_b_7TeV          lnN               1.06         1.06         -            1.06         1.06         -          " << endl;
  datacard <<      "CMS_resbkgcat0          lnN               -            5.0          -            -            -            -      " << endl;    
  datacard <<      "CMS_resbkgcat1          lnN               -            -            -            -            5.0          -      " << endl;    
  datacard <<      "CMS_nonresbkgcat0       lnN               -            -            5.0          -            -            -     " << endl;     
  datacard <<      "CMS_nonresbkgcat1       lnN               -            -            -            -            -            5.0   " << endl;       
  datacard <<      "" << endl;


}



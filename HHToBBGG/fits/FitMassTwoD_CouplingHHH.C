//
// root -l -b -q CMSAna/HHToBBGG/fits/FitMassTwoD.C+'("/afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/FitWorkspace.root", "MassFitTwoD.root", 1, 1, 0, 5000)'
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
TLatex *tex = 0;
int alreadyPlotted = 0;

void MakePlots(RooWorkspace *ws, RooFitResult *fitResult2D);

//-------------------------------------------------------------
//Main macro for generating data and fitting
//=============================================================  
void FitMassTwoD(const string workspaceFile = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/FitWorkspace_asdf.root", const string outputTree = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/MassFitResults/MassFitTwoD_asdf.root", Int_t plotOption = 0, Int_t storeOption = 1, Int_t constBkg = 1, Int_t NToys = 5000) {

  TRandom3 *randomnumber = new TRandom3(1200);
  TFile *wsFile = new TFile (workspaceFile.c_str(), "READ");
  RooWorkspace *ws = (RooWorkspace*)wsFile->Get("MassFitWorkspace");
  
  //Import variables from workspace
  RooAbsPdf *model2Dpdf = ws->pdf("model2Dpdf");
  RooRealVar *massBjet = ws->var("massBjet");
  RooRealVar *massPho = ws->var("massPho");
  RooRealVar *HHHCoupling = ws->var("HHHCoupling"); RooRealVar constHHHCoupling(*HHHCoupling);
  RooRealVar *nres = ws->var("N (ResBkg)"); RooRealVar constNres(*nres);
  RooRealVar *nnonres = ws->var("N (NonResBkg)"); RooRealVar constNnonres(*nnonres);
  RooRealVar *expRateBjet = ws->var("expRateBjet"); RooRealVar constexpBjet(*expRateBjet);
  RooRealVar *expRatePho = ws->var("expRatePho"); RooRealVar constexpPho(*expRatePho);
  RooFormulaVar *nsig = (RooFormulaVar*)ws->function("nsig"); 

  if (constBkg == 0) { 
    expRateBjet->setConstant(false);
    expRatePho->setConstant(false);
  } else {
    expRateBjet->setConstant(true);
    expRatePho->setConstant(true);
  }

  //Variables to set constant
  RooRealVar *sigMeanBjet = ws->var("sigMeanBjet"); sigMeanBjet->setConstant();
  RooRealVar *sigSigmaBjet = ws->var("sigSigmaBjet"); sigSigmaBjet->setConstant();
  RooRealVar *sigAlpha = ws->var("sigAlpha"); sigAlpha->setConstant();
  RooRealVar *sigPower = ws->var("sigPower"); sigPower->setConstant();
  RooRealVar *resMeanBjet = ws->var("resMeanBjet"); resMeanBjet->setConstant();
  RooRealVar *resSigmaBjet = ws->var("resSigmaBjet"); resSigmaBjet->setConstant();
  RooRealVar *resAlpha = ws->var("resAlpha"); resAlpha->setConstant();
  RooRealVar *resPower = ws->var("resPower"); resPower->setConstant();
  RooRealVar *resExpo = ws->var("resExpo"); resExpo->setConstant();
  RooRealVar *nbbH = ws->var("nbbH"); nbbH->setConstant();
  RooRealVar *nOthers = ws->var("nOthers"); nOthers->setConstant();
  
  //Create TTree to store the resulting yield data
  TFile *f = new TFile(outputTree.c_str(), "RECREATE");
  TTree *resultTree = new TTree("resultTree", "Parameter results from fitting");
  Float_t nsigOut, nresOut, nnonresOut;
  Float_t nsigStd, nresStd, nnonresStd;
  Float_t HHHCouplingOut, HHHCouplingStd;
  
  resultTree->Branch("nsigOut",&nsigOut, "nsigOut/F");
  resultTree->Branch("nresOut",&nresOut, "nresOut/F");
  resultTree->Branch("nnonresOut",&nnonresOut, "nnonresOut/F");
  resultTree->Branch("nsigStd",&nsigStd, "nsigStd/F");
  resultTree->Branch("nresStd",&nresStd, "nresStd/F");
  resultTree->Branch("nnonresStd",&nnonresStd, "nnonresStd/F");
  resultTree->Branch("HHHCouplingOut",&HHHCouplingOut, "HHHCouplingOut/F");
  resultTree->Branch("HHHCouplingStd",&HHHCouplingStd, "HHHCouplingStd/F");
  
  //Generate Toy MC experiment data and fits
  for(UInt_t t=0; t < NToys; ++t) {
    cout << "Toy #" << t << endl;
    HHHCoupling->setVal(constHHHCoupling.getVal()); nres->setVal(constNres.getVal()); nnonres->setVal(constNnonres.getVal());
    expRateBjet->setVal(constexpBjet.getVal()); expRatePho->setVal(constexpPho.getVal());

    Float_t ran = randomnumber->Poisson(325.);
    RooDataSet *pseudoData2D = model2Dpdf->generate(RooArgList(*massBjet,*massPho), ran);
    RooFitResult *fitResult2D = model2Dpdf->fitTo(*pseudoData2D, RooFit::Save(), RooFit::Extended(kTRUE), RooFit::Strategy(2));

//     if (t == 1763) {
//     	ws->import(*pseudoData2D, kTRUE);
//     	ws->import(*pseudoData2D, Rename("pseudoData2D"));
//     }
//     if (plotOption == 1) MakePlots(ws, fitResult2D);
    

//     cout << "DEBUG: " << constexpBjet.getVal() << " , " << constexpPho.getVal() << " | " << expRateBjet->getVal() << " " << expRatePho->getVal() << "\n";
    cout << "DEBUG: " << HHHCoupling->getVal() << " " << HHHCoupling->getPropagatedError(*fitResult2D) << "\n";


    //Store fit parameters into ROOT file
    if (storeOption == 1) {
      //Save variables into separate branches of root tree
      
      nsigOut = nsig->getVal();
      nresOut = nres->getVal();
      nnonresOut = nnonres->getVal();
      nsigStd = nsig->getPropagatedError(*fitResult2D);
      nresStd = nres->getPropagatedError(*fitResult2D);
      nnonresStd = nnonres->getPropagatedError(*fitResult2D);
      HHHCouplingOut = HHHCoupling->getVal();
      HHHCouplingStd = HHHCoupling->getPropagatedError(*fitResult2D);

      //cout << "\n\n\n\n\n\n" << nsigOut << " | " << nresOut << " | " << nnonresOut << " | " << ran  << "\n\n\n\n\n" << endl;
      resultTree->Fill();
    }
  }
  
  //Write to the TTree and close it
  resultTree->Write();
  f->Close();
}

//-------------------------------------------------------------
//Plot the model fitting function for signal+background
//=============================================================
void MakePlots(RooWorkspace *ws, RooFitResult *fitResult2D) {
  
  RooPlot* framex = 0;
  RooPlot* framey = 0;
  
  //Import yield variables
  RooRealVar *nsig = ws->var("N (Sig)");
  RooRealVar *nres = ws->var("N (ResBkg)");
  RooRealVar *nnonres = ws->var("N (NonResBkg)");
  //Select and plot only one experiment
  if (! (nsig->getVal() >= 15.6 && nsig->getVal() <= 17.0
  				&& nres->getVal() >= 27.9 && nres->getVal() <= 28.0
  				&& nnonres->getVal() <= 287.))
  	return;
  if (alreadyPlotted) return;  
  alreadyPlotted = 1;
  
  //Import the PDF's
  RooAbsPdf *model2Dpdf = ws->pdf("model2Dpdf");
  RooAbsPdf *sigPDFPho = ws->pdf("sigPDFPho");
  RooAbsPdf *resPDFPho = ws->pdf("resPDFPho");
  RooAbsPdf *nonresPDFPho = ws->pdf("nonresPDFPho");
  RooAbsPdf *sigPDFBjet = ws->pdf("sigPDFBjet");
  RooAbsPdf *sigPDFBjetCut = ws->pdf("sigPDFBjetCut");
  RooAbsPdf *resPDFBjet = ws->pdf("resPDFBjet");
  RooAbsPdf *resPDFBjetExt = ws->pdf("resPDFBjetExt");
  RooAbsPdf *nonresPDFBjet = ws->pdf("nonresPDFBjet");
  //x-axis variables
  RooRealVar *massPho = ws->var("massPho");
  RooRealVar *massBjet = ws->var("massBjet");
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
  RooDataSet *pseudoData2D = (RooDataSet *)ws->data("pseudoData2D");
  
  //Plot of 2D generated data and fits
  TH1 *data2d = pseudoData2D->createHistogram("2D Data", *massBjet,Binning(25), YVar(*massPho,Binning(25)));
  data2d->SetStats(0);
  TH1 *fit2d = model2Dpdf->createHistogram("2D Fit", *massBjet, YVar(*massPho));
  fit2d->SetStats(0);
  
  cv = new TCanvas("cv","cv",1600,600);
  cv->Divide(2);
  cv->cd(1); gPad->SetLeftMargin(0.15); data2d->Draw("LEGO2");
  data2d->SetTitle("");
  data2d->GetXaxis()->SetTitleOffset(2); data2d->GetXaxis()->SetTitle("M_{bb} [GeV/c^{2}]");
  data2d->GetYaxis()->SetTitleOffset(2.2); data2d->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
  data2d->GetZaxis()->SetTitleOffset(1.75);
  
  cv->cd(2); gPad->SetLeftMargin(0.15); fit2d->Draw("SURF1");
  fit2d->SetTitle("");
  fit2d->GetXaxis()->SetTitleOffset(2); data2d->GetXaxis()->SetTitle("M_{bb} [GeV/c^{2}]");
  fit2d->GetYaxis()->SetTitleOffset(2.2); data2d->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
  fit2d->GetZaxis()->SetTitleOffset(1.75);
  cv->SaveAs(Form("Plots/AllSignalBkgd/Fits/Toys/twoDimensionalFits_%d.gif",1));
  
  //Plot of massBjet and massPho projections from 2D fit
  massPho->setRange("PhoWindow",120,130);
  framex = massBjet->frame(Bins(50)); 
  framex->SetTitle(""); framex->SetXTitle("M_{bb} [GeV/c^{2}]");  framex->SetYTitle("Number of Events");
  pseudoData2D->plotOn(framex, CutRange("PhoWindow"));
  model2Dpdf->plotOn(framex, ProjectionRange("PhoWindow"), VisualizeError(*fitResult2D), FillStyle(3001));
  model2Dpdf->plotOn(framex, ProjectionRange("PhoWindow"));
  model2Dpdf->plotOn(framex, ProjectionRange("PhoWindow"), Components("sigPDFBjet"), LineStyle(kDashed), LineColor(kRed));
  model2Dpdf->plotOn(framex, ProjectionRange("PhoWindow"), Components("resPDFBjet"), LineStyle(kDashed), LineColor(kOrange));
  model2Dpdf->plotOn(framex, ProjectionRange("PhoWindow"), Components("nonresPDFBjet"), LineStyle(kDashed), LineColor(kGreen));
  
  framey = massPho->frame(Bins(50)); 
  framey->SetTitle(""); framey->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");  framey->SetYTitle("Number of Events");
  pseudoData2D->plotOn(framey);
  model2Dpdf->plotOn(framey, VisualizeError(*fitResult2D), FillStyle(3001));
  model2Dpdf->plotOn(framey);
  model2Dpdf->plotOn(framey,Components("sigPDFPho"), LineStyle(kDashed), LineColor(kRed));
  model2Dpdf->plotOn(framey,Components("resPDFPho"), LineStyle(kDashed), LineColor(kOrange));
  model2Dpdf->plotOn(framey,Components("nonresPDFPho"), LineStyle(kDashed), LineColor(kGreen));
  
  cv = new TCanvas("cv","cv",1600,600);
  cv->Divide(2);
  cv->cd(1); framex->Draw();
  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.042);
  tex->SetTextFont(42);
  tex->DrawLatex(0.52, 0.84, Form("N_{Sig} = %.2f +/- %.2f", nsig->getVal(), nsig->getPropagatedError(*fitResult2D)));
  tex->DrawLatex(0.52, 0.79, Form("N_{ResBkg} = %.2f +/- %.2f", nres->getVal(), nres->getPropagatedError(*fitResult2D)));
  tex->DrawLatex(0.52, 0.74, Form("N_{NonResBkg} = %.2f +/- %.2f", nnonres->getVal(), nnonres->getPropagatedError(*fitResult2D)));
  tex->Draw();
  cv->Update();
  cv->cd(2); framey->Draw();
  tex->DrawLatex(0.52, 0.84, Form("N_{Sig} = %.2f +/- %.2f", nsig->getVal(), nsig->getPropagatedError(*fitResult2D)));
  tex->DrawLatex(0.52, 0.79, Form("N_{ResBkg} = %.2f +/- %.2f", nres->getVal(), nres->getPropagatedError(*fitResult2D)));
  tex->DrawLatex(0.52, 0.74, Form("N_{NonResBkg} = %.2f +/- %.2f", nnonres->getVal(), nnonres->getPropagatedError(*fitResult2D)));
  tex->Draw();
  cv->Update();
  cv->SaveAs(Form("Plots/AllSignalBkgd/Fits/Toys/projectionFits_%d.gif",1));
}



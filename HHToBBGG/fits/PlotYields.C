//
//
// root -l CMSAna/HHToBBGG/fits/PlotYields.C+'("/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/MassFitResults/ResolutionAnalysis/MassFitTwoD_ResStep0.9.root")'
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
#include <TGraph.h>									// plotting class
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

TCanvas *cv = 0;
TGraphErrors *gr = 0;
TLine *line = 0;
Float_t ResAnalysisNum = 20.0;

void plotFit(TH1F *hist, TH1F *histPercent, TH1F *histErr, TF1 *func, const string outputName);
void plotErrVsRes();

void PlotYields(const string inputfile = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/MassFitResults/ResolutionAnalysis/MassFitTwoD.root") {
	
	plotErrVsRes();
	/*
  //Import the yield data
	TFile *file = new TFile(inputfile.c_str(), "READ");
	TTree *nTree = (TTree*)file->Get("resultTree");
  Float_t nsigOut, nresOut, nnonresOut;
  Float_t nsigStd, nresStd, nnonresStd;
  
  nTree->SetBranchAddress("nsigOut",&nsigOut);
  nTree->SetBranchAddress("nresOut",&nresOut);
  nTree->SetBranchAddress("nnonresOut",&nnonresOut);
  nTree->SetBranchAddress("nsigStd",&nsigStd);
  nTree->SetBranchAddress("nresStd",&nresStd);
  nTree->SetBranchAddress("nnonresStd",&nnonresStd);
  
  //Create histograms for the yield pull distributions
  TH1F *sigPull = new TH1F("sigPull", ";Signal Pull;Number of Events", 50, -8,8);
  sigPull->SetLineColor(kRed); sigPull->SetLineWidth(4);
  TH1F *resPull = new TH1F("resPull", ";Resonant Background Pull;Number of Events", 50, -8,8);
  resPull->SetLineColor(kOrange); resPull->SetLineWidth(4);
  TH1F *nonresPull = new TH1F("nonresPull", ";Non-Resonant Background Pull;Number of Events", 50, -8,8);
  nonresPull->SetLineColor(kGreen); nonresPull->SetLineWidth(4);
  //Create histograms for percentage of pull means
  TH1F *sigPercent = new TH1F("sigPercent", ";Signal Percent;Number of Events", 50, -1,1);
  sigPercent->SetLineColor(kRed); sigPercent->SetLineWidth(4);
  TH1F *resPercent = new TH1F("resPercent", ";Resonant Background Percent;Number of Events", 50, -1,1);
  resPercent->SetLineColor(kOrange); resPercent->SetLineWidth(4);
  TH1F *nonresPercent = new TH1F("nonresPercent", ";Non-Resonant Background Percent;Number of Events", 50, -1,1);
  nonresPercent->SetLineColor(kGreen); nonresPercent->SetLineWidth(4);
  //Create histograms for relative error
  TH1F *sigErr = new TH1F("sigErr", ";Signal Error;Number of Events", 50, -1,1);
  sigErr->SetLineColor(kRed); sigErr->SetLineWidth(4);
  TH1F *resErr = new TH1F("resErr", ";Resonant Background Error;Number of Events", 50, -1,1);
  resErr->SetLineColor(kOrange); resErr->SetLineWidth(4);
  TH1F *nonresErr = new TH1F("nonresErr", ";Non-Resonant Background Error;Number of Events", 50, -1,1);
  nonresErr->SetLineColor(kGreen); nonresErr->SetLineWidth(4);
  
  //Fill pull histograms
  for (int i=0; i < nTree->GetEntries(); i++) {
    nTree->GetEntry(i);
  	Float_t signal = (nsigOut - 16.5) / nsigStd; //15.85, 12.3, 16.5
  	Float_t resbkg = (nresOut - 27.5) / nresStd; //26.2, 17.6, 27.5
  	Float_t nonresbkg = (nnonresOut - 284.4) / nnonresStd; //53.7, 95.6, 284.4
  	Float_t signalPercent = (nsigOut - 16.5) / 16.5;
  	Float_t resbkgPercent = (nresOut - 27.5) / 27.5;
  	Float_t nonresbkgPercent = (nnonresOut - 284.4) / 284.4;
  	
  	sigPull->Fill(signal);
  	resPull->Fill(resbkg);
  	nonresPull->Fill(nonresbkg);
  	sigPercent->Fill(signalPercent);
  	resPercent->Fill(resbkgPercent);
  	nonresPercent->Fill(nonresbkgPercent);
  	sigErr->Fill(nsigStd/16.5);
  	resErr->Fill(nresStd/27.5);
  	nonresErr->Fill(nnonresStd/284.4); 
  }
  
  //Store pull histograms
  TFile *outfile = new TFile (Form("CMSAna/HHToBBGG/data/MassFitResults/ResolutionAnalysis/YieldsTwoD_ResStep%.1f.root", ResAnalysisNum), "RECREATE");
  outfile->cd();
  outfile->WriteTObject(sigPull, sigPull->GetName(), "WriteDelete");
  outfile->WriteTObject(resPull, resPull->GetName(), "WriteDelete");
  outfile->WriteTObject(nonresPull, nonresPull->GetName(), "WriteDelete");
  outfile->WriteTObject(sigPercent, sigPercent->GetName(), "WriteDelete");
  outfile->WriteTObject(resPercent, resPercent->GetName(), "WriteDelete");
  outfile->WriteTObject(nonresPercent, nonresPercent->GetName(), "WriteDelete");
  outfile->WriteTObject(sigErr, sigErr->GetName(), "WriteDelete");
  outfile->WriteTObject(resErr, resErr->GetName(), "WriteDelete");
  outfile->WriteTObject(nonresErr, nonresErr->GetName(), "WriteDelete");
  outfile->Close();
  delete outfile;
  
  //Fit histograms to gaussians
  TF1 *sigFit = new TF1("N (Sig) Fit","gaus", -8, 8);
  TF1 *resFit = new TF1("N (ResBkg) Fit","gaus", -8, 8);
  TF1 *nonresFit = new TF1("N (NonResBkg) Fit","gaus", -8, 8);
  sigPull->Fit(sigFit,"R");
  resPull->Fit(resFit,"R");
  nonresPull->Fit(nonresFit,"R");
  
  plotFit(sigPull,sigPercent,sigErr,sigFit,"sigPullFit");
  plotFit(resPull,resPercent,resErr,resFit,"resPullFit");
  plotFit(nonresPull,nonresPercent,nonresErr,nonresFit,"nonresPullFit");
  */
}

void plotFit(TH1F *hist, TH1F *histPercent, TH1F *histErr, TF1 *func, const string outputName) {
  cv = new TCanvas("cv","cv", 800,600);
  cv->SetFillStyle(4000);
  hist->Draw("HIST"); 
  hist->SetStats(kFALSE);
  hist->GetYaxis()->SetTitleOffset(1.5);
  std::cout << hist->Integral() << std::endl;
  func->SetLineColor(kGray+3);
  func->SetLineWidth(3);
  func->Draw("SAME");
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.032);
  tex->SetTextFont(2);
  tex->DrawLatex(0.69, 0.86, Form("#bar{N} = %.3f +/- %.3f", func->GetParameter(1), func->GetParError(1)));
  tex->DrawLatex(0.69, 0.82, Form("#bar{N}/#mu_{N} = %.3f +/- %.3f", histPercent->GetMean(), histPercent->GetRMS()));
  tex->DrawLatex(0.69, 0.78, Form("#sigma_{N} = %.3f +/- %.3f", func->GetParameter(2), func->GetParError(2)));
  tex->DrawLatex(0.69, 0.74, Form("#sigma_{N}/#mu_{N} = %.3f +/- %.3f", histErr->GetMean(), histErr->GetRMS()));
  tex->Draw();
  cv->Update();
  cv->SaveAs(("Plots/AllSignalBkgd/Fits/PullPlots/"+outputName+"TwoD_ResStep20.0.gif").c_str());
}

void plotErrVsRes() {

	Float_t xpho[11] = {0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9};
	Float_t ypho[11] = {.430, .436, .442, .447, .453, .459, .464, .470, .475, .480, .485};
	Float_t expho[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	Float_t eypho[11] = {.043, .044, .044, .044, .044, .045, .045, .045, .045, .045, .045};
	
	Float_t xbjet[11] = {10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20.};
	Float_t ybjet[11] = {.376, .390, .404, .419, .434, .451, .470, .488, .507, .528, .550};
	Float_t exbjet[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	Float_t eybjet[11] = {.040, .039, .040, .042, .042, .043, .045, .046, .047, .049, .049};
	
	//Draw
	cv = new TCanvas("cv","cv", 800,600);
	gr = new TGraphErrors(11,xpho,ypho,expho,eypho);
	gr->SetTitle("Relative Error vs. diPhoton Resolution (Signal)");
	gr->SetMinimum(0.0); gr->SetMaximum(1.0);
	gr->GetXaxis()->SetTitle("M_{#gamma#gamma} Signal Width");
	gr->GetYaxis()->SetTitle("#sigma_{N}/#mu_{N}");
	gr->SetMarkerColor(2);
	gr->SetMarkerStyle(3);
	gr->Draw("ALP");
	TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.024);
  tex->DrawLatex(0.54, 0.82, "#splitline{Current resolution}{from MC histogram fit}");
  tex->Draw();
	cv->Update();
	line = new TLine(1.43,0,1.43,1);
	line->SetLineWidth(2); line->SetLineColor(kBlue);
	line->Draw();
	cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/ResolutionPhotonSteps.gif");
	
	cv = new TCanvas("cv","cv", 800,600);
	gr = new TGraphErrors(11,xbjet,ybjet,exbjet,eybjet);
	gr->SetTitle("Relative Error vs. di-Bjet Resolution (Signal)");
	gr->SetMinimum(0.0); gr->SetMaximum(1.0);
	gr->GetXaxis()->SetTitle("M_{bb} Signal Width");
	gr->GetYaxis()->SetTitle("#sigma_{N}/#mu_{N}");
	gr->SetMarkerColor(2);
	gr->SetMarkerStyle(3);
	gr->Draw("ALP");
	TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.024);
  tex->DrawLatex(0.55, 0.82, "#splitline{Current resolution}{from MC histogram fit}");
  tex->Draw();
	cv->Update();
	line = new TLine(15.49,0,15.49,1);
	line->SetLineWidth(2); line->SetLineColor(kBlue);
	line->Draw();
	cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/ResolutionJetEnergySteps.gif");
}




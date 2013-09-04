//
//
// root -l CMSAna/HHToBBGG/fits/PlotScenarioYields.C+'("pho")'
//
// scanOption can be = "pho", "jet", "lum"
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

TCanvas *cv = 0;
TGraphErrors *gr = 0;
TLine *line = 0;
TLatex *tex = 0;
int numSteps = 0;
Float_t luminosity[13] = {0.5/3., 1./3., 1.5/3., 2./3., 2.5/3., 3./3., 4./3., 5./3., 6./3., 7./3., 8./3., 9./3., 10./3.}; 

void plotErrVsRes(const string scanOption, Float_t relativeError[], Float_t sigmaOfRelativeError[]);

void PlotScenarioYields(const string scanOption = "pho") {
	
  if (scanOption == "lum") numSteps = 13;
  else numSteps = 13;
  Float_t relativeError[numSteps];
  Float_t sigmaOfRelativeError[numSteps];
  for (int s = 0; s < numSteps; s++) {
    
    //Import the yield data
    TFile *file = new TFile(Form(("CMSAna/HHToBBGG/data/MassFitResults/ResolutionAnalysis/"+scanOption+"FitScenario%d.root").c_str(), s), "READ");
    TTree *nTree = (TTree*)file->Get("resultTree");
    Float_t nsigOut, nresOut, nnonresOut;
    Float_t nsigStd, nresStd, nnonresStd;
    
    nTree->SetBranchAddress("nsigOut",&nsigOut);
    nTree->SetBranchAddress("nresOut",&nresOut);
    nTree->SetBranchAddress("nnonresOut",&nnonresOut);
    nTree->SetBranchAddress("nsigStd",&nsigStd);
    nTree->SetBranchAddress("nresStd",&nresStd);
    nTree->SetBranchAddress("nnonresStd",&nnonresStd);
    
    //Create histogram for relative error
    TH1F *sigErr = new TH1F("sigErr", ";Signal Error;Number of Events", 100, -1,1);
    sigErr->SetLineColor(kRed); sigErr->SetLineWidth(4);
    
    Float_t nsigActual;
    if (scanOption == "lum") nsigActual = 16.3*luminosity[s];
    else nsigActual = 16.3;
    //Fill pull histograms
    for (int i=0; i < nTree->GetEntries(); i++) {
      nTree->GetEntry(i);
      sigErr->Fill(nsigStd/nsigActual);
    }
    
    relativeError[s] = sigErr->GetMean()*100;
    sigmaOfRelativeError[s] = sigErr->GetRMS()*100;
    file->Close();
  }
	
  for (int s = 0; s < 11; s++) {
    cout << relativeError[s] << " +/- " << sigmaOfRelativeError[s] << endl;
  }
  plotErrVsRes(scanOption, relativeError, sigmaOfRelativeError);
}

void plotErrVsRes(const string scanOption, Float_t relativeError[], Float_t sigmaOfRelativeError[]) {

  if (scanOption == "pho") {
    //diPhoton points for graph
    Float_t xpho[13] = {0.9, 1.15, 1.4, 1.65, 1.9, 2.15, 2.4, 2.65, 2.9, 3.15, 3.4, 3.65, 3.9};
    Float_t expho[13] = {0};
    
    //Draw diPhoton Resolution Steps
    cv = new TCanvas("cv","cv", 800,600);
    gr = new TGraphErrors(11, xpho,relativeError, expho,sigmaOfRelativeError);
    gr->SetTitle(""); gr->SetMinimum(0.0); gr->SetMaximum(100.);
    gr->GetXaxis()->SetTitle("M_{#gamma#gamma} Signal Width [GeV/c^{2}]");
    gr->GetYaxis()->SetTitle("Relative Uncertainty on Fitted Signal Yield");
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->SetMarkerColor(2); gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
    gr->Draw("ALP");
    tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.040);
    tex->SetTextFont(42);
    tex->DrawLatex(0.37, 0.70, "Nominal Resolution");
    tex->Draw();
    cv->Update();
    line = new TLine(1.43,0,1.43,100.);
    line->SetLineWidth(3); line->SetLineColor(kBlue);
    line->Draw();
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/ResolutionPhotonSteps.gif");
  }
	
  else if (scanOption == "jet") {
    //diBjet points for graph
    Float_t xbjet[13] = {10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0};
    Float_t exbjet[13] = {0};
    
    //Draw diBjet Resolution Steps
    cv = new TCanvas("cv","cv", 800,600);
    gr = new TGraphErrors(11, xbjet,relativeError, exbjet,sigmaOfRelativeError);
    gr->SetTitle(""); gr->SetMinimum(0.0); gr->SetMaximum(100.);
    gr->GetXaxis()->SetTitle("M_{bb} Signal Width [GeV/c^{2}]");
    gr->GetYaxis()->SetTitle("Relative Uncertainty on Fitted Signal Yield");
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->SetMarkerColor(2); gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
    gr->Draw("ALP");
    tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.040);
    tex->SetTextFont(42);
    tex->DrawLatex(0.5, 0.85, "Resolution at 140 Pileup");
    tex->Draw();
    cv->Update();
    line = new TLine(20.0,0,20.0,100.);
    line->SetLineWidth(3); line->SetLineColor(kBlue);
    line->Draw();
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/ResolutionJetEnergySteps.gif");
  }
  
  else {
    //Luminosity points for graph
    Float_t x[13] = {0.5, 1., 1.5, 2., 2.5, 3., 4., 5., 6., 7., 8., 9., 10.};
    Float_t ex[13] = {0};
    
    //Draw Luminosity Steps
    cv = new TCanvas("cv","cv", 800,600);
    gr = new TGraphErrors(13, x,relativeError, ex,sigmaOfRelativeError);
    gr->SetTitle(""); gr->SetMinimum(0.0); gr->SetMaximum(100.);
    gr->GetXaxis()->SetTitle("Integrated Luminosity [10^{3} fb^{-1}]");
    gr->GetYaxis()->SetTitle("Relative Uncertainty on Fitted Signal Yield");
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->SetMarkerColor(2); gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
    gr->Draw("ALP");
    tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.040);
    tex->SetTextFont(42);
    tex->DrawLatex(0.34, 0.82, "Nominal Luminosity");
    tex->Draw();
    cv->Update();
    line = new TLine(3.,0,3.,100.);
    line->SetLineWidth(3); line->SetLineColor(kBlue);
    line->Draw();
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/ResolutionLumiSteps.gif");
  }

}




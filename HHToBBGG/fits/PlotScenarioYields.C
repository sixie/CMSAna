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
#include <TGraphErrors.h>						                // plotting error class
#include <TLine.h>									// plotting lines
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TLegend.h>  
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
Float_t luminosity[13] = {1.5/3., 2./3., 2.5/3., 3./3., 4./3., 5./3., 6./3., 7./3., 8./3., 9./3., 10./3., 11./3., 12./3.}; 

void plotErrVsRes(const string scanOption, Float_t relativeError1[], Float_t sigmaOfRelativeError1[], Float_t relativeError2[], Float_t sigmaOfRelativeError2[]);

void PlotScenarioYields(const string scanOption = "pho") {
	
  if (scanOption == "lum") numSteps = 13;
  else if (scanOption == "pho") numSteps = 13;
  else if (scanOption == "jet") numSteps = 13;
  else numSteps = 13;
  Float_t relativeError1[numSteps];
  Float_t sigmaOfRelativeError1[numSteps];
  Float_t relativeError2[numSteps];
  Float_t sigmaOfRelativeError2[numSteps];

  for (int s = 0; s < numSteps; s++) {
    
    //Import the yield data
    TFile *file = new TFile(Form(("CMSAna/HHToBBGG/data/MassFitResults/ResolutionAnalysis/OptionB_HighPileupEfficiency/"+scanOption+"FitScenario%d.root").c_str(), s), "READ");
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
    TH1F *sigErr = new TH1F("sigErr", ";Signal Error;Number of Events", 100, 0,5);
    sigErr->SetLineColor(kRed); sigErr->SetLineWidth(4);
    
    Float_t nsigActual;
    if (scanOption == "lum") nsigActual = 4.93*luminosity[s];
    else nsigActual = 4.93;
    //Fill pull histograms
    for (int i=0; i < nTree->GetEntries(); i++) {
      nTree->GetEntry(i);
      sigErr->Fill(nsigStd/nsigActual);
    }
    
    relativeError1[s] = sigErr->GetMean()*100;
    sigmaOfRelativeError1[s] = sigErr->GetRMS()*100;
    file->Close();
  }
	
  for (int s = 0; s < numSteps; s++) {
      
    //Import the yield data
    TFile *file = new TFile(Form(("CMSAna/HHToBBGG/data/MassFitResults/ResolutionAnalysis/OptionB_LumiTimes2/"+scanOption+"FitScenario%d.root").c_str(), s), "READ");
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
    TH1F *sigErr = new TH1F("sigErr", ";Signal Error;Number of Events", 100, 0,5);
    sigErr->SetLineColor(kRed); sigErr->SetLineWidth(4);
    
    Float_t nsigActual;
    if (scanOption == "lum") nsigActual = 9.86*luminosity[s];
    else nsigActual = 9.86;
    //Fill pull histograms
    for (int i=0; i < nTree->GetEntries(); i++) {
      nTree->GetEntry(i);
      sigErr->Fill(nsigStd/nsigActual);
    }
    
    relativeError2[s] = sigErr->GetMean()*100;
    sigmaOfRelativeError2[s] = sigErr->GetRMS()*100;
    file->Close();
  }
	
  for (int s = 0; s < numSteps; s++) {
    cout << relativeError1[s] << " +/- " << sigmaOfRelativeError1[s] << endl;
  }
  for (int s = 0; s < numSteps; s++) {
    cout << relativeError2[s] << " +/- " << sigmaOfRelativeError2[s] << endl;
  }

  plotErrVsRes(scanOption, relativeError1, sigmaOfRelativeError1, relativeError2, sigmaOfRelativeError2);
}

void plotErrVsRes(const string scanOption, Float_t relativeError1[], Float_t sigmaOfRelativeError1[], Float_t relativeError2[], Float_t sigmaOfRelativeError2[]) {

  if (scanOption == "pho") {
    //diPhoton points for graph
    Float_t xpho[13] = {0.9, 1.15, 1.4, 1.65, 1.9, 2.15, 2.4, 2.65, 2.9, 3.15, 3.4, 3.65, 3.9};
    Float_t expho[13] = {0};
    
    //Draw diPhoton Resolution Steps
    cv = new TCanvas("cv","cv", 800,600);
    gr = new TGraphErrors(11, xpho,relativeError1, expho,sigmaOfRelativeError1);
    gr->SetTitle(""); gr->SetMinimum(0.0); gr->SetMaximum(175.);
    gr->GetXaxis()->SetTitle("M_{#gamma#gamma} Resolution Parameter [GeV/c^{2}]");
    gr->GetYaxis()->SetTitle("Relative Uncertainty on Fitted Signal Yield (%)");
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->SetMarkerColor(kRed); gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
    gr->SetLineColor(kRed); gr->SetLineWidth(2);
    gr->SetFillStyle(3004); gr->SetFillColor(kRed);
    gr->Draw("ALPE3");

    TGraphErrors *gr2 = new TGraphErrors(11, xpho,relativeError2, expho,sigmaOfRelativeError2);
    gr2->SetTitle(""); gr2->SetMinimum(0.0); gr2->SetMaximum(175.);
    gr2->GetXaxis()->SetTitle("M_{#gamma#gamma} Resolution Parameter [GeV/c^{2}]");
    gr2->GetYaxis()->SetTitle("Relative Uncertainty on Fitted Signal Yield (%)");
    gr2->GetYaxis()->SetTitleSize(0.04);
    gr2->SetMarkerColor(kGreen+2); gr2->SetMarkerStyle(20); gr2->SetMarkerSize(1.2);
    gr2->SetLineColor(kGreen+2); gr2->SetLineWidth(2);
    gr2->SetFillStyle(3004); gr2->SetFillColor(kGreen+2);
    gr2->Draw("LPE3,same");

    TLegend *legend = new TLegend(0.45,0.20,0.80,0.40);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(gr,"Nominal Signal Efficiency", "LP");
    legend->AddEntry(gr2,"x2 Signal Efficiency", "LP");
    legend->Draw();

    tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.040);
    tex->SetTextFont(42);
    tex->SetTextColor(kBlue);
    tex->DrawLatex(0.37, 0.80, "Nominal Resolution");
    tex->Draw();
    cv->Update();
    line = new TLine(1.43,0,1.43,175.);
    line->SetLineWidth(3); line->SetLineColor(kBlue);
    line->Draw();
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/ResolutionPhotonSteps.gif");
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/ResolutionPhotonSteps.pdf");
  }
	
  else if (scanOption == "jet") {
    //diBjet points for graph
    Float_t xbjet[13] = {10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0};
    Float_t exbjet[13] = {0};
    
    //Draw diBjet Resolution Steps
    cv = new TCanvas("cv","cv", 800,600);
    gr = new TGraphErrors(13, xbjet,relativeError1, exbjet,sigmaOfRelativeError1);
    gr->SetTitle(""); gr->SetMinimum(0.0); gr->SetMaximum(250.);
    gr->GetXaxis()->SetTitle("M_{bb} Resolution Parameter [GeV/c^{2}]");
    gr->GetYaxis()->SetTitle("Relative Uncertainty on Fitted Signal Yield (%)");
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->SetMarkerColor(kRed); gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
    gr->SetLineColor(kRed); gr->SetLineWidth(2);
    gr->SetFillStyle(3004); gr->SetFillColor(kRed);
    gr->Draw("ALPE3");

    TGraphErrors *gr2 = new TGraphErrors(13, xbjet,relativeError2, exbjet,sigmaOfRelativeError2);
    gr2->SetTitle(""); gr2->SetMinimum(0.0); gr2->SetMaximum(250.);
    gr2->GetXaxis()->SetTitle("M_{bb} Resolution Parameter [GeV/c^{2}]");
    gr2->GetYaxis()->SetTitle("Relative Uncertainty on Fitted Signal Yield (%)");
    gr2->GetYaxis()->SetTitleSize(0.04);
    gr2->SetMarkerColor(kGreen+2); gr2->SetMarkerStyle(20); gr2->SetMarkerSize(1.2);
    gr2->SetLineColor(kGreen+2); gr2->SetLineWidth(2);
    gr2->SetFillStyle(3004); gr2->SetFillColor(kGreen+2);
    gr2->Draw("LPE3,same");

    TLegend *legend = new TLegend(0.55,0.20,0.75,0.35);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(gr,"Nominal Signal Efficiency", "LP");
    legend->AddEntry(gr2,"x2 Signal Efficiency", "LP");
    legend->Draw();

    tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.040);
    tex->SetTextFont(42);
    tex->SetTextColor(kBlack);
    tex->DrawLatex(0.45, 0.85, "Resolution");
    tex->DrawLatex(0.45, 0.80, "(140 Pileup)");
    tex->Draw();
    cv->Update();
    line = new TLine(20.0,0,20.0,250.);
    line->SetLineWidth(3); line->SetLineColor(kBlack);
    line->Draw();
    TLatex *tex2 = new TLatex();
    tex2->SetNDC();
    tex2->SetTextSize(0.040);
    tex2->SetTextFont(42);
    tex2->SetTextColor(kBlue);
    tex2->DrawLatex(0.16, 0.85, "Resolution");
    tex2->DrawLatex(0.16, 0.80, "(LHC Run1)");
    tex2->Draw();
    TLine *line2 = new TLine(14.29,0,14.29,250.);
    line2->SetLineWidth(3); line2->SetLineColor(kBlue);
    line2->Draw();
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/ResolutionJetEnergySteps.gif");
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/ResolutionJetEnergySteps.pdf");
  }
  
  else {
    //Luminosity points for graph
    Float_t x[13] = { 1.5, 2., 2.5, 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.};
    Float_t ex[13] = {0};
    
    //Draw Luminosity Steps
    cv = new TCanvas("cv","cv", 800,600);
    gr = new TGraphErrors(13, x,relativeError1, ex,sigmaOfRelativeError1);
    gr->SetTitle(""); gr->SetMinimum(0.0); gr->SetMaximum(200.);
    gr->GetXaxis()->SetTitle("Integrated Luminosity [10^{3} fb^{-1}]");
    gr->GetYaxis()->SetTitle("Relative Uncertainty on Fitted Signal Yield (%)");
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->SetMarkerColor(kRed); gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
    gr->SetLineColor(kRed); gr->SetLineWidth(2);
    gr->SetFillStyle(3004); gr->SetFillColor(kRed);
    gr->Draw("ALPE3");

    TGraphErrors *gr2 = new TGraphErrors(13, x,relativeError2, ex,sigmaOfRelativeError2);
    gr2->SetTitle(""); gr2->SetMinimum(0.0); gr2->SetMaximum(200.);
    gr2->GetXaxis()->SetTitle("M_{bb} Resolution Parameter [GeV/c^{2}]");
    gr2->GetYaxis()->SetTitle("Relative Uncertainty on Fitted Signal Yield (%)");
    gr2->GetYaxis()->SetTitleSize(0.04);
    gr2->SetMarkerColor(kGreen+2); gr2->SetMarkerStyle(20); gr2->SetMarkerSize(1.2);
    gr2->SetLineColor(kGreen+2); gr2->SetLineWidth(2);
    gr2->SetFillStyle(3004); gr2->SetFillColor(kGreen+2);
    gr2->Draw("LPE3,same");

    TLegend *legend = new TLegend(0.55,0.64,0.80,0.84);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(gr,"Nominal Signal Efficiency", "LP");
    legend->AddEntry(gr2,"x2 Signal Efficiency", "LP");
    legend->Draw();


    tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.040);
    tex->SetTextFont(42);
    tex->SetTextColor(kBlue);
    tex->DrawLatex(0.33, 0.22, "Nominal Luminosity");
    tex->Draw();
    cv->Update();
    line = new TLine(3.,0,3.,200.);
    line->SetLineWidth(3); line->SetLineColor(kBlue);
    line->Draw();
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/ResolutionLumiSteps.gif");
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/ResolutionLumiSteps.pdf");
  }

}




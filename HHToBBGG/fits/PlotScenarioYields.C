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
Float_t phoEffRatio[6] = {0.9, 1.0, 1.1, 1.2, 1.3, 1.4};
Float_t btagEffRatio[6] = {0.9, 1.0, 1.1, 1.2, 1.3, 1.4};


void plotErrVsRes(const string scanOption, Float_t relativeError1[], Float_t sigmaOfRelativeError1[], Float_t relativeError2[], Float_t sigmaOfRelativeError2[]);

void PlotScenarioYields(const string scanOption = "pho") {
	
  if (scanOption == "lum") numSteps = 13;
  else if (scanOption == "pho") numSteps = 13;
  else if (scanOption == "jet") numSteps = 13;
  else if (scanOption == "endcapPhotonFakerate") numSteps = 6;
  else if (scanOption == "photonFakerate") numSteps = 6;
  else if (scanOption == "phoEff") numSteps = 6;
  else if (scanOption == "btagEff") numSteps = 6;
  else numSteps = 13;
  Float_t relativeError1[numSteps];
  Float_t sigmaOfRelativeError1[numSteps];
  Float_t relativeError2[numSteps];
  Float_t sigmaOfRelativeError2[numSteps];

  double myPhoEffRatio = 1.0;
  double myBtagEffRatio = 1.0;

  for (int s = 0; s < numSteps; s++) {

    //Import the yield data
    TFile *file = 0;
    if (scanOption == "endcapPhotonFakerate") {
      file = new TFile(Form(("CMSAna/HHToBBGG/data/MassFitResults/ResolutionAnalysis/OptionB_Categories/baseline140PU/"+scanOption+"FitScenario%d.root").c_str(), s), "READ");     
    } else {
      file = new TFile(Form(("CMSAna/HHToBBGG/data/MassFitResults/ResolutionAnalysis/OptionB_Categories/baseline140PU/"+scanOption+"FitScenario%d.root").c_str(), s), "READ");
    }
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
    
    if (scanOption == "phoEff") {
      myPhoEffRatio = phoEffRatio[s];
      myBtagEffRatio = 1.25;
    }
    if (scanOption == "btagEff") {
      myBtagEffRatio = btagEffRatio[s];
      myPhoEffRatio = 1.25; 
    }

    Float_t nsigActual;
    if (scanOption == "lum") nsigActual = 6.11*luminosity[s]*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio;
    else nsigActual = 6.11*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio;
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
    
    if (scanOption == "phoEff") {
      myBtagEffRatio  = 1.25;
    }
    if (scanOption == "btagEff") {
      myPhoEffRatio = 1.25;
    }
    

    //Import the yield data
    TFile *file = 0;
    if (scanOption == "photonFakerate") {
      file = new TFile(Form(("CMSAna/HHToBBGG/data/MassFitResults/ResolutionAnalysis/OptionB_Categories/backup/"+scanOption+"FitScenario%d.root").c_str(), s), "READ");
    } else {
      file = new TFile(Form(("CMSAna/HHToBBGG/data/MassFitResults/ResolutionAnalysis/OptionB_Categories/improved/"+scanOption+"FitScenario%d.root").c_str(), s), "READ");
    }
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
    
    myPhoEffRatio = 1.25;
    myBtagEffRatio = 1.25;

    Float_t nsigActual;
    if (scanOption == "lum") nsigActual = 6.11*luminosity[s]*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio;
    else nsigActual = 6.11*myPhoEffRatio*myPhoEffRatio*myBtagEffRatio*myBtagEffRatio;
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
    gr->GetYaxis()->SetTitle("Relative Uncertainty on Fitted Signal Yield [%]");
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->SetMarkerColor(kRed); gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
    gr->SetLineColor(kRed); gr->SetLineWidth(2);
    gr->SetFillStyle(3004); gr->SetFillColor(kRed);
    gr->Draw("ALPE3");

    TGraphErrors *gr2 = new TGraphErrors(11, xpho,relativeError2, expho,sigmaOfRelativeError2);
    gr2->SetTitle(""); gr2->SetMinimum(0.0); gr2->SetMaximum(175.);
    gr2->GetXaxis()->SetTitle("M_{#gamma#gamma} Resolution Parameter [GeV/c^{2}]");
    gr2->GetYaxis()->SetTitle("Relative Uncertainty on Fitted Signal Yield [%]");
    gr2->GetYaxis()->SetTitleSize(0.04);
    gr2->SetMarkerColor(kGreen+2); gr2->SetMarkerStyle(20); gr2->SetMarkerSize(1.2);
    gr2->SetLineColor(kGreen+2); gr2->SetLineWidth(2);
    gr2->SetFillStyle(3004); gr2->SetFillColor(kGreen+2);
    gr2->Draw("LPE3,same");

    TLegend *legend = new TLegend(0.45,0.18,0.80,0.35);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(gr,"Nominal Performance", "LP");
    legend->AddEntry(gr2,"Improved Performance", "LP");
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
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSUncertaintyVsPhotonResolution.gif");
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSUncertaintyVsPhotonResolution.pdf");
  }
	
  else if (scanOption == "jet") {
    //diBjet points for graph
    Float_t xbjet[13] = {10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0};
    Float_t exbjet[13] = {0};
    
    //Draw diBjet Resolution Steps
    cv = new TCanvas("cv","cv", 800,600);
    gr = new TGraphErrors(13, xbjet,relativeError1, exbjet,sigmaOfRelativeError1);
    gr->SetTitle(""); gr->SetMinimum(0.0); gr->SetMaximum(220.);
    gr->GetXaxis()->SetTitle("M_{bb} Resolution Parameter [GeV/c^{2}]");
    gr->GetYaxis()->SetTitle("Relative Uncertainty on Fitted Signal Yield [%]");
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->SetMarkerColor(kRed); gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
    gr->SetLineColor(kRed); gr->SetLineWidth(2);
    gr->SetFillStyle(3004); gr->SetFillColor(kRed);
    gr->Draw("ALPE3");

    TGraphErrors *gr2 = new TGraphErrors(13, xbjet,relativeError2, exbjet,sigmaOfRelativeError2);
    gr2->SetTitle(""); gr2->SetMinimum(0.0); gr2->SetMaximum(220.);
    gr2->GetXaxis()->SetTitle("M_{bb} Resolution Parameter [GeV/c^{2}]");
    gr2->GetYaxis()->SetTitle("Relative Uncertainty on Fitted Signal Yield [%]");
    gr2->GetYaxis()->SetTitleSize(0.04);
    gr2->SetMarkerColor(kGreen+2); gr2->SetMarkerStyle(20); gr2->SetMarkerSize(1.2);
    gr2->SetLineColor(kGreen+2); gr2->SetLineWidth(2);
    gr2->SetFillStyle(3004); gr2->SetFillColor(kGreen+2);
    gr2->Draw("LPE3,same");

    TLegend *legend = new TLegend(0.45,0.20,0.75,0.35);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(gr,"Nominal Performance", "LP");
    legend->AddEntry(gr2,"Improved Performance", "LP");
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
    line = new TLine(20.0,0,20.0,220.);
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
    TLine *line2 = new TLine(14.29,0,14.29,220.);
    line2->SetLineWidth(3); line2->SetLineColor(kBlue);
    line2->Draw();
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSUncertaintyVsJetEnergyResolution.gif");
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSUncertaintyVsJetEnergyResolution.pdf");
  }
  else if (scanOption == "endcapPhotonFakerate") {
    //diBjet points for graph
    Float_t x[6] = {0, 0.0001, 0.0003, 0.001, 0.003, 0.01}; 
    Float_t ex[6] = {0};
    
    //Draw diBjet Resolution Steps
    cv = new TCanvas("cv","cv", 800,600);
    gr = new TGraphErrors(6, x,relativeError1, ex,sigmaOfRelativeError1);
    gr->SetTitle(""); gr->SetMinimum(0.0); gr->SetMaximum(150.);
    gr->GetXaxis()->SetTitle("Photon Fake Rate In Endcap");
    gr->GetYaxis()->SetTitle("Relative Uncertainty on Fitted Signal Yield [%]");
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->SetMarkerColor(kRed); gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
    gr->SetLineColor(kRed); gr->SetLineWidth(2);
    gr->SetFillStyle(3004); gr->SetFillColor(kRed);
    gr->Draw("ALPE3");

    TGraphErrors *gr2 = new TGraphErrors(6, x,relativeError2, ex,sigmaOfRelativeError2);
    gr2->SetTitle(""); gr2->SetMinimum(0.0); gr2->SetMaximum(150.);
    gr2->GetXaxis()->SetTitle("Photon Fake Rate In Endcap");
    gr2->GetYaxis()->SetTitle("Relative Uncertainty on Fitted Signal Yield [%]");
    gr2->GetYaxis()->SetTitleSize(0.04);
    gr2->SetMarkerColor(kGreen+2); gr2->SetMarkerStyle(20); gr2->SetMarkerSize(1.2);
    gr2->SetLineColor(kGreen+2); gr2->SetLineWidth(2);
    gr2->SetFillStyle(3004); gr2->SetFillColor(kGreen+2);
    gr2->Draw("LPE3,same");

    TLegend *legend = new TLegend(0.20,0.20,0.45,0.35);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(gr,"Nominal High Pileup Performance", "LP");
    legend->AddEntry(gr2,"All Improvements Applied", "LP");
    legend->Draw();

    tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.040);
    tex->SetTextFont(42);
    tex->SetTextColor(kBlack);
    tex->DrawLatex(0.60, 0.85, "Nominal Endcap");
    tex->DrawLatex(0.60, 0.80, "Fake Rate");
    tex->Draw();
    cv->Update();
    line = new TLine(0.003,0,0.003,150.);
    line->SetLineWidth(3); line->SetLineColor(kBlack);
    line->Draw();
//     TLatex *tex2 = new TLatex();
//     tex2->SetNDC();
//     tex2->SetTextSize(0.040);
//     tex2->SetTextFont(42);
//     tex2->SetTextColor(kBlue);
//     tex2->DrawLatex(0.16, 0.85, "Resolution");
//     tex2->DrawLatex(0.16, 0.80, "(LHC Run1)");
//     tex2->Draw();
//     TLine *line2 = new TLine(14.29,0,14.29,250.);
//     line2->SetLineWidth(3); line2->SetLineColor(kBlue);
//     line2->Draw();
    cv->SetLogx(true);
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSUncertaintyVsPhotonFakerate.gif");
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSUncertaintyVsPhotonFakerate.pdf");
  }
   else if (scanOption == "phoEff") {
    //diBjet points for graph
    Float_t x[6] = {-10, 0, 10, 20, 30, 40};
    Float_t ex[6] = {0};
    
    //Draw diBjet Resolution Steps
    cv = new TCanvas("cv","cv", 800,600);
    gr = new TGraphErrors(6, x,relativeError1, ex,sigmaOfRelativeError1);
    gr->SetTitle(""); gr->SetMinimum(0.0); gr->SetMaximum(100.);
    gr->GetXaxis()->SetTitle("Relative Improvement In Photon Selection Efficiency [%]");
    gr->GetYaxis()->SetTitle("Relative Uncertainty on Fitted Signal Yield [%]");
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->SetMarkerColor(kRed); gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
    gr->SetLineColor(kRed); gr->SetLineWidth(2);
    gr->SetFillStyle(3004); gr->SetFillColor(kRed);
    gr->Draw("ALPE3");

    TLegend *legend = new TLegend(0.20,0.25,0.45,0.40);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(gr,"All Other Improvements Applied", "LP");
    legend->Draw();

    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSUncertaintyVsPhotonEffRatio.gif");
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSUncertaintyVsPhotonEffRatio.pdf");
  }
   else if (scanOption == "btagEff") {
    //diBjet points for graph
    Float_t x[6] = {-10, 0, 10, 20, 30, 40};
    Float_t ex[6] = {0};
    
    //Draw diBjet Resolution Steps
    cv = new TCanvas("cv","cv", 800,600);
    gr = new TGraphErrors(6, x,relativeError1, ex,sigmaOfRelativeError1);
    gr->SetTitle(""); gr->SetMinimum(0.0); gr->SetMaximum(100.);
    gr->GetXaxis()->SetTitle("Relative Improvement In B-Tagging Efficiency [%]");
    gr->GetYaxis()->SetTitle("Relative Uncertainty on Fitted Signal Yield [%]");
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->SetMarkerColor(kRed); gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
    gr->SetLineColor(kRed); gr->SetLineWidth(2);
    gr->SetFillStyle(3004); gr->SetFillColor(kRed);
    gr->Draw("ALPE3");

    TLegend *legend = new TLegend(0.20,0.25,0.45,0.40);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(gr,"All Other Improvements Applied", "LP");
    legend->Draw();

    cv->Update();
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSUncertaintyVsBtagEffRatio.gif");
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSUncertaintyVsBtagEffRatio.pdf");
  }
  
  else if (scanOption == "lum"){
    //Luminosity points for graph
    Float_t x[13] = { 1.5, 2., 2.5, 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.};
    Float_t ex[13] = {0};
    
    //Draw Luminosity Steps
    cv = new TCanvas("cv","cv", 800,600);

    gr = new TGraphErrors(13, x,relativeError1, ex,sigmaOfRelativeError1);
    gr->SetTitle(""); gr->SetMinimum(0.0); gr->SetMaximum(200.);
    gr->GetXaxis()->SetTitle("Integrated Luminosity [10^{3} fb^{-1}]");
    gr->GetYaxis()->SetTitle("Relative Uncertainty on Fitted Signal Yield [%]");
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->SetMarkerColor(kRed); gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
    gr->SetLineColor(kRed); gr->SetLineWidth(2);
    gr->SetFillStyle(3004); gr->SetFillColor(kRed);
    gr->Draw("ALPE3");

    TGraphErrors *gr2 = new TGraphErrors(13, x,relativeError2, ex,sigmaOfRelativeError2);
    gr2->SetTitle(""); gr2->SetMinimum(0.0); gr2->SetMaximum(200.);
    gr2->GetXaxis()->SetTitle("M_{bb} Resolution Parameter [GeV/c^{2}]");
    gr2->GetYaxis()->SetTitle("Relative Uncertainty on Fitted Signal Yield [%]");
    gr2->GetYaxis()->SetTitleSize(0.04);
    gr2->SetMarkerColor(kGreen+2); gr2->SetMarkerStyle(20); gr2->SetMarkerSize(1.2);
    gr2->SetLineColor(kGreen+2); gr2->SetLineWidth(2);
    gr2->SetFillStyle(3004); gr2->SetFillColor(kGreen+2);
    gr2->Draw("LPE3,same");

    TLegend *legend = new TLegend(0.55,0.64,0.80,0.84);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(gr,"Nominal Performance", "LP");
    legend->AddEntry(gr2,"Improved Performance", "LP");
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
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSUncertaintyVsLumi.gif");
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSUncertaintyVsLumi.pdf");
  }

}




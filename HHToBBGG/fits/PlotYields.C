//
//
// root -l CMSAna/HHToBBGG/fits/PlotYields.C+'("/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/MassFitResults/ResolutionAnalysis/DoubleLumi_ResStep0.9.root")'
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
TLatex *tex = 0;

void plotFit(TH1F *hist, TH1F *histPercent, TH1F *histErr, TF1 *func, const string outputName);

void PlotYields(const string inputfile = "/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/MassFitResults/MassFitTwoD.root") {
	
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
  TH1F *sigPercent = new TH1F("sigPercent", ";Signal Percent;Number of Events", 50, -5,5);
  sigPercent->SetLineColor(kRed); sigPercent->SetLineWidth(4);
  TH1F *resPercent = new TH1F("resPercent", ";Resonant Background Percent;Number of Events", 50, -1,1);
  resPercent->SetLineColor(kOrange); resPercent->SetLineWidth(4);
  TH1F *nonresPercent = new TH1F("nonresPercent", ";Non-Resonant Background Percent;Number of Events", 50, -1,1);
  nonresPercent->SetLineColor(kGreen); nonresPercent->SetLineWidth(4);
  
  //Create histograms for relative error
  TH1F *sigErr = new TH1F("sigErr", ";Signal Error;Number of Events", 100, 0,5);
  sigErr->SetLineColor(kRed); sigErr->SetLineWidth(4);
  TH1F *resErr = new TH1F("resErr", ";Resonant Background Error;Number of Events", 50, 0,5);
  resErr->SetLineColor(kOrange); resErr->SetLineWidth(4);
  TH1F *nonresErr = new TH1F("nonresErr", ";Non-Resonant Background Error;Number of Events", 50, 0,5);
  nonresErr->SetLineColor(kGreen); nonresErr->SetLineWidth(4);
  
  Float_t nsigActual = 4.93;
  //Float_t nsigActual = 9.86;
  //Fill pull histograms
  for (int i=0; i < nTree->GetEntries(); i++) {
    nTree->GetEntry(i);
    //This gives the Z score (how many std. dev.'s the measurement is from the actual)
  	Float_t signal = (nsigOut - nsigActual) / nsigStd; //15.85, 12.3, 16.5
  	//Float_t resbkg = (nresOut - 27.5) / nresStd; //26.2, 17.6, 27.5
  	//Float_t nonresbkg = (nnonresOut - 284.4) / nnonresStd; //53.7, 95.6, 284.4
  	
  	//This is relative error (how many % off from the mean)
//   	Float_t signalPercent = (nsigOut - nsigActual) / nsigActual;
  	Float_t signalPercent = (nsigOut - 4.93) / 4.93;

  	//Float_t resbkgPercent = (nresOut - 27.5) / 27.5;
  	//Float_t nonresbkgPercent = (nnonresOut - 284.4) / 284.4;
  	
  	sigPull->Fill(signal);
  	//resPull->Fill(resbkg);
  	//nonresPull->Fill(nonresbkg);
  	sigPercent->Fill(signalPercent);
  	//resPercent->Fill(resbkgPercent);
  	//nonresPercent->Fill(nonresbkgPercent);
  	
  	//This is the relative uncertainty
  	sigErr->Fill(nsigStd/nsigActual);
  	//resErr->Fill(nresStd/27.5);
  	//nonresErr->Fill(nnonresStd/284.4); 
  	//plotErrVsRes();
  }
  
  /*
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
  */
  
  //Fit histograms to gaussians
  TF1 *sigFit = new TF1("N (Sig) Fit","gaus", -1.5, 8);
//   TF1 *sigFit = new TF1("N (Sig) Fit","gaus", -5, 8);
//   TF1 *resFit = new TF1("N (ResBkg) Fit","gaus", -5, 8);
  //TF1 *nonresFit = new TF1("N (NonResBkg) Fit","gaus", -8, 8);
  sigPull->Fit(sigFit,"R");
  //resPull->Fit(resFit,"R");
  //nonresPull->Fit(nonresFit,"R");
  
  plotFit(sigPull,sigPercent,sigErr,sigFit,"sigPullFit");
  //plotFit(resPull,resPercent,resErr,resFit,"resPullFit");
  //plotFit(nonresPull,nonresPercent,nonresErr,nonresFit,"nonresPullFit");
  
  
}

void plotFit(TH1F *hist, TH1F *histPercent, TH1F *histErr, TF1 *func, const string outputName) {
  cv = new TCanvas("cv","cv", 800,600);
  cv->SetFillStyle(4000);
  hist->SetFillStyle(0);
  hist->Draw("HIST");
  hist->SetStats(kFALSE);
  hist->GetYaxis()->SetTitle("Number of Toy Experiments"); hist->GetYaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitle("(N^{signal}_{fit} - N^{signal}_{input}) / #sigma_{fit}"); hist->GetXaxis()->SetRangeUser(-6,8);
  std::cout << hist->Integral() << std::endl;
  std::cout << "Pull Mean: " << hist->GetMean() << "\n";
  std::cout << "Pull RMS: " << hist->GetRMS() << "\n";
  std::cout << "Bias: " << histPercent->GetMean() << "\n";

   func->SetLineColor(kGray+3);
   func->SetLineWidth(3);
   func->Draw("SAME");

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.042);
  tex->SetTextFont(42);
  tex->DrawLatex(0.55, 0.84, Form("Pull Mean = %.3f +/- %.3f", func->GetParameter(1), func->GetParError(1)));
  tex->DrawLatex(0.55, 0.79, Form("Pull Width = %.3f +/- %.3f", func->GetParameter(2), func->GetParError(2)));
  tex->DrawLatex(0.55, 0.74, Form("Avg Bias = %.1f%%", histPercent->GetMean()*100.));
  tex->Draw();
  cv->Update();
  cv->SaveAs(("Plots/AllSignalBkgd/Fits/PullPlots/"+outputName+".gif").c_str());
  cv->SaveAs(("Plots/AllSignalBkgd/Fits/PullPlots/"+outputName+".pdf").c_str());
  
  cv = new TCanvas("cv","cv", 800,600);
  cv->SetFillStyle(4000);
  histErr->Draw("HIST");
  histErr->SetStats(kFALSE);
  histErr->GetYaxis()->SetTitle("Number of Toy Experiments"); histErr->GetYaxis()->SetTitleOffset(1.3);
  histErr->GetXaxis()->SetTitle("Relative Uncertainty"); histErr->GetXaxis()->SetRangeUser(0,2);
  TLatex *tex2 = new TLatex();
  tex2->SetNDC();
  tex2->SetTextSize(0.042);
  tex2->SetTextFont(42);
  tex2->DrawLatex(0.17, 0.79, Form("#splitline{Mean of Relative}{         Uncertainty}  = %.0f%%", histErr->GetMean()*100.));
  tex2->Draw();
  cv->Update();
  cv->SaveAs(("Plots/AllSignalBkgd/Fits/PullPlots/"+outputName+"_RelErr.gif").c_str());
  cv->SaveAs(("Plots/AllSignalBkgd/Fits/PullPlots/"+outputName+"_RelErr.pdf").c_str());
}




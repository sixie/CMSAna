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


Float_t expLimitVsLuminosityBaseline[13] = {5.79688,
4.10938,
2.64844,
2.27344,
1.44922,
0.839844,
0.716797,
0.548828,
0.417969,
0.383789,
0.325195,
0.264648,
0.239258,
}; 

Float_t expLimitSigmaVsLuminosityBaseline[13] = {2.95764,
2.09666,
1.35127,
1.15994,
0.739411,
0.4285,
0.365719,
0.28002,
0.213253,
0.195814,
0.165919,
0.135027,
0.122073,
}; 

Float_t expLimitVsPhotonResolutionBaseline[13] = {1.96094,
2.05469,
2.17969,
2.55469,
2.42969,
2.39844,
2.50781,
2.64844,
2.69531,
2.74219,
2.89844,
3.03906,
2.86719,
}; 
Float_t expLimitSigmaVsPhotonResolutionBaseline[13] = {1.0005,
1.04833,
1.11211,
1.30344,
1.23966,
1.22372,
1.27952,
1.35127,
1.37518,
1.3991,
1.47882,
1.55057,
1.46288,
}; 

Float_t expLimitVsJetEnergyResolutionBaseline[13] = {1.50391,
1.66406,
1.89844,
2.08594,
2.27344,
2.38281,
2.60938,
2.86719,
2.96094,
3.17188,
3.32812,
3.48438,
3.60938,                                   
}; 
Float_t expLimitSigmaVsJetEnergyResolutionBaseline[13] = {0.767313,
0.849027,
0.968608,
1.06427,
1.15994,
1.21574,
1.33134,
1.46288,
1.51071,
1.61833,
1.69805,
1.77778,
1.84155,                           
}; 

Float_t expLimitVsLuminosityImproved[13] = {3.11719,
1.88281,
1.13672,
1.01172,
0.638672,
0.427734,
0.354492,
0.282227,
0.233398,
0.179199,
0.157715,
0.138184,
0.121582,
};
Float_t expLimitSigmaVsLuminosityImproved[13] = {1.59043,
0.960636,
0.579969,
0.516193,
0.325859,
0.218236,
0.180867,
0.143996,
0.119083,
0.0914299,
0.0804682,
0.0705031,
0.0620328,
};

Float_t expLimitVsPhotonResolutionImproved[13] = {0.832031,
0.933594,
1.07422,
1.09766,
1.13672,
1.08203,
1.14453,
1.19922,
1.20703,
1.26953,
1.20703,
1.33984,
1.31641,
}; 
Float_t expLimitSigmaVsPhotonResolutionImproved[13] = {0.424514,
0.476332,
0.548081,
0.560039,
0.579969,
0.552067,
0.583955,
0.611858,
0.615844,
0.647732,
0.615844,
0.683606,
0.671648,
}; 

Float_t expLimitVsJetEnergyResolutionImproved[13] = {0.808594,
0.910156,
1.05078,
1.05859,
1.22266,
1.25391,
1.37109,
1.55078,
1.61719,
1.69531,
1.89844,
1.99219,
1.86719,
}; 
Float_t expLimitSigmaVsJetEnergyResolutionImproved[13] = {0.412555,
0.464374,
0.536123,
0.540109,
0.623816,
0.63976,
0.699551,
0.79123,
0.825111,
0.864971,
0.968608,
1.01644,
0.952664,
}; 



Float_t expLimitVsEndcapPhotonFakerateBaseline[6] = {1.92969,
1.85156,
1.85156,
1.75781,
2.27344,
1.97656,
};
Float_t expLimitSigmaVsEndcapPhotonFakerateBaseline[6] = {0.984553,
0.944692,
0.944692,
0.89686,
1.15994,
1.00847,
};
Float_t expLimitVsEndcapPhotonFakerateImproved[6] = {0.871094,
0.878906,
0.878906,
0.917969,
1.01172,
1.01172,
};
Float_t expLimitSigmaVsEndcapPhotonFakerateImproved[6] = {0.444444,
0.44843,
0.44843,
0.46836,
0.516193,
0.516193,
};
Float_t expLimitVsPhotonEffRatio[6] = {1.58594,
1.28516,
1.08984,
1.01172,
0.902344,
0.847656,
}; 
Float_t expLimitSigmaVsPhotonEffRatio[6] = {0.809167,
0.655704,
0.556053,
0.516193,
0.460388,
0.432486,
}; 
Float_t expLimitVsBTagEffRatio[6] = {1.41016,
1.28516,
1.05859,
1,
0.964844,
0.855469,
};
Float_t expLimitSigmaVsBTagEffRatio[6] = {0.719481,
0.655704,
0.540109,
0.510213,
0.492276,
0.436472,
};



void plotErrVsRes(const string scanOption, Float_t relativeError1[], Float_t sigmaOfRelativeError1[], Float_t relativeError2[], Float_t sigmaOfRelativeError2[]);

void PlotScenarioLimits(const string scanOption = "pho") {
	
  const double HHXS = 33.9;
  if (scanOption == "lum") numSteps = 13;
  else if (scanOption == "pho") numSteps = 13;
  else if (scanOption == "jet") numSteps = 13;
  else if (scanOption == "endcapPhotonFakerate") numSteps = 6;
  else if (scanOption == "photonFakerate") numSteps = 6;
  else if (scanOption == "phoEff") numSteps = 6;
  else if (scanOption == "btagEff") numSteps = 6;
  else numSteps = 13;
  Float_t expectedLimit1[numSteps];
  Float_t sigmaOfExpectedLimit1[numSteps];
  Float_t expectedLimit2[numSteps];
  Float_t sigmaOfExpectedLimit2[numSteps];

  for (int s = 0; s < numSteps; s++) {   
    if (scanOption == "lum") {
      expectedLimit1[s] = expLimitVsLuminosityBaseline[s] * HHXS;
      sigmaOfExpectedLimit1[s] = expLimitSigmaVsLuminosityBaseline[s] * HHXS;
    } else if (scanOption == "pho") {
      expectedLimit1[s] = expLimitVsPhotonResolutionBaseline[s] * HHXS;
      sigmaOfExpectedLimit1[s] = expLimitSigmaVsPhotonResolutionBaseline[s] * HHXS;
    } else if (scanOption == "jet") {
      expectedLimit1[s] = expLimitVsJetEnergyResolutionBaseline[s] * HHXS;
      sigmaOfExpectedLimit1[s] = expLimitSigmaVsJetEnergyResolutionBaseline [s] * HHXS;
    }  else if (scanOption == "endcapPhotonFakerate") {
      expectedLimit1[s] = expLimitVsEndcapPhotonFakerateBaseline[s] * HHXS;
      sigmaOfExpectedLimit1[s] = expLimitSigmaVsEndcapPhotonFakerateBaseline[s] * HHXS;
    } else if (scanOption == "phoEff") {
      expectedLimit1[s] = expLimitVsPhotonEffRatio[s] * HHXS;
      sigmaOfExpectedLimit1[s] = expLimitSigmaVsPhotonEffRatio[s] * HHXS;
    } else if (scanOption == "btagEff") {
      expectedLimit1[s] = expLimitVsBTagEffRatio[s] * HHXS;
      sigmaOfExpectedLimit1[s] = expLimitSigmaVsBTagEffRatio[s] * HHXS;
    } 
  }
	
  for (int s = 0; s < numSteps; s++) {      
    if (scanOption == "lum") {
      expectedLimit2[s] = expLimitVsLuminosityImproved[s] * HHXS;
      sigmaOfExpectedLimit2[s] = expLimitSigmaVsLuminosityImproved[s] * HHXS;
    } else if (scanOption == "pho") {
      expectedLimit2[s] = expLimitVsPhotonResolutionImproved[s] * HHXS;
      sigmaOfExpectedLimit2[s] = expLimitSigmaVsPhotonResolutionImproved[s] * HHXS;
    } else if (scanOption == "jet") {
      expectedLimit2[s] = expLimitVsJetEnergyResolutionImproved[s] * HHXS;
      sigmaOfExpectedLimit2[s] = expLimitSigmaVsJetEnergyResolutionImproved[s] * HHXS;
    }  else if (scanOption == "endcapPhotonFakerate") {
      expectedLimit2[s] = expLimitVsEndcapPhotonFakerateImproved[s] * HHXS;
      sigmaOfExpectedLimit2[s] = expLimitSigmaVsEndcapPhotonFakerateImproved[s] * HHXS;
    } 

  }
	
  for (int s = 0; s < numSteps; s++) {
    cout << expectedLimit1[s] << " +/- " << sigmaOfExpectedLimit1[s] << endl;
  }
  for (int s = 0; s < numSteps; s++) {
    cout << expectedLimit2[s] << " +/- " << sigmaOfExpectedLimit2[s] << endl;
  }

  plotErrVsRes(scanOption, expectedLimit1, sigmaOfExpectedLimit1, expectedLimit2, sigmaOfExpectedLimit2);
}

void plotErrVsRes(const string scanOption, Float_t relativeError1[], Float_t sigmaOfRelativeError1[], Float_t relativeError2[], Float_t sigmaOfRelativeError2[]) {

  if (scanOption == "pho") {
    //diPhoton points for graph
    Float_t xpho[13] = {0.9, 1.15, 1.4, 1.65, 1.9, 2.15, 2.4, 2.65, 2.9, 3.15, 3.4, 3.65, 3.9};
    Float_t expho[13] = {0};
    
    //Draw diPhoton Resolution Steps
    cv = new TCanvas("cv","cv", 800,600);
    cv->SetLogy();
    gr = new TGraphErrors(13, xpho,relativeError1, expho,sigmaOfRelativeError1);
    gr->SetTitle(""); gr->SetMinimum(0.0); gr->SetMaximum(1000.);
    gr->GetXaxis()->SetTitle("M_{#gamma#gamma} Resolution Parameter [GeV/c^{2}]");
    gr->GetYaxis()->SetTitle("Aymptotic CLs 95% Cross Section Limit [fb]");
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->SetMarkerColor(kRed); gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
    gr->SetLineColor(kRed); gr->SetLineWidth(2);
    gr->SetFillStyle(3004); gr->SetFillColor(kRed);
    gr->Draw("ALPE3");

    TGraphErrors *gr2 = new TGraphErrors(13, xpho,relativeError2, expho,sigmaOfRelativeError2);
    gr2->SetTitle(""); gr2->SetMinimum(0.0); gr2->SetMaximum(1000.);
    gr2->GetXaxis()->SetTitle("M_{#gamma#gamma} Resolution Parameter [GeV/c^{2}]");
    gr2->GetYaxis()->SetTitle("Aymptotic CLs 95% Cross Section Limit [fb]");
    gr2->GetYaxis()->SetTitleSize(0.04);
    gr2->SetMarkerColor(kGreen+2); gr2->SetMarkerStyle(20); gr2->SetMarkerSize(1.2);
    gr2->SetLineColor(kGreen+2); gr2->SetLineWidth(2);
    gr2->SetFillStyle(3004); gr2->SetFillColor(kGreen+2);
    gr2->Draw("LPE3,same");

    TLegend *legend = new TLegend(0.35,0.20,0.80,0.40);
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
    tex->DrawLatex(0.37, 0.85, "Nominal Resolution");
    tex->Draw();
    cv->Update();
    line = new TLine(1.43,0,1.43,1000.);
    line->SetLineWidth(3); line->SetLineColor(kBlue);
    line->Draw();

    TLine *line2 = new TLine(0.6,33.9,4.2,33.9);
    line2->SetLineWidth(4); line2->SetLineColor(kBlack);
    line2->Draw();
    TLatex *tex2 = new TLatex();
    tex2->SetNDC();
    tex2->SetTextSize(0.05);
    tex2->SetTextFont(42);
    tex2->SetTextColor(kBlack);
    tex2->DrawLatex(0.67, 0.45, "SM Cross Section");
    tex2->Draw();

    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSLimitVsPhotonResolution.gif");
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSLimitVsPhotonResolution.pdf");
  }
	
  else if (scanOption == "jet") {
    //diBjet points for graph
    Float_t xbjet[13] = {10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0};
    Float_t exbjet[13] = {0};
    
    //Draw diBjet Resolution Steps
    cv = new TCanvas("cv","cv", 800,600);
    cv->SetLogy();
    gr = new TGraphErrors(13, xbjet,relativeError1, exbjet,sigmaOfRelativeError1);
    gr->SetTitle(""); gr->SetMinimum(0.0); gr->SetMaximum(1000.);
    gr->GetXaxis()->SetTitle("M_{bb} Resolution Parameter [GeV/c^{2}]");
    gr->GetYaxis()->SetTitle("Aymptotic CLs 95% Cross Section Limit [fb]");
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->SetMarkerColor(kRed); gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
    gr->SetLineColor(kRed); gr->SetLineWidth(2);
    gr->SetFillStyle(3004); gr->SetFillColor(kRed);
    gr->Draw("ALPE3");

    TGraphErrors *gr2 = new TGraphErrors(13, xbjet,relativeError2, exbjet,sigmaOfRelativeError2);
    gr2->SetTitle(""); gr2->SetMinimum(0.0); gr2->SetMaximum(1000.);
    gr2->GetXaxis()->SetTitle("M_{bb} Resolution Parameter [GeV/c^{2}]");
    gr2->GetYaxis()->SetTitle("Aymptotic CLs 95% Cross Section Limit [fb]");
    gr2->GetYaxis()->SetTitleSize(0.04);
    gr2->SetMarkerColor(kGreen+2); gr2->SetMarkerStyle(20); gr2->SetMarkerSize(1.2);
    gr2->SetLineColor(kGreen+2); gr2->SetLineWidth(2);
    gr2->SetFillStyle(3004); gr2->SetFillColor(kGreen+2);
    gr2->Draw("LPE3,same");

    TLegend *legend = new TLegend(0.60,0.20,0.85,0.40);
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
    tex->SetTextColor(kMagenta+2);
    tex->DrawLatex(0.45, 0.85, "Resolution");
    tex->DrawLatex(0.45, 0.80, "(140 Pileup)");
    tex->Draw();
    cv->Update();
    line = new TLine(20.0,0,20.0,1000.);
    line->SetLineWidth(3); line->SetLineColor(kMagenta+2);
    line->Draw();
    TLatex *tex2 = new TLatex();
    tex2->SetNDC();
    tex2->SetTextSize(0.040);
    tex2->SetTextFont(42);
    tex2->SetTextColor(kBlue);
    tex2->DrawLatex(0.16, 0.85, "Resolution");
    tex2->DrawLatex(0.16, 0.80, "(LHC Run1)");
    tex2->Draw();
    TLine *line2 = new TLine(14.29,0,14.29,1000.);
    line2->SetLineWidth(3); line2->SetLineColor(kBlue);
    line2->Draw();

    TLine *line3 = new TLine(7,33.9,43,33.9);
    line3->SetLineWidth(4); line3->SetLineColor(kBlack);
    line3->Draw();
    TLatex *tex3 = new TLatex();
    tex3->SetNDC();
    tex3->SetTextSize(0.05);
    tex3->SetTextFont(42);
    tex3->SetTextColor(kBlack);
    tex3->DrawLatex(0.67, 0.45, "SM Cross Section");
    tex3->Draw();

    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSLimitVsJetEnergyResolution.gif");
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSLimitVsJetEnergyResolution.pdf");
  }
  else if (scanOption == "endcapPhotonFakerate") {
    //diBjet points for graph
    Float_t x[6] = {0, 0.0001, 0.0003, 0.001, 0.003, 0.01}; 
    Float_t ex[6] = {0};
    
    //Draw diBjet Resolution Steps
    cv = new TCanvas("cv","cv", 800,600);
    gr = new TGraphErrors(6, x,relativeError1, ex,sigmaOfRelativeError1);
    gr->SetTitle(""); gr->SetMinimum(0.0); gr->SetMaximum(200.);
    gr->GetXaxis()->SetTitle("Photon Fake Rate In Endcap");
    gr->GetYaxis()->SetTitle("Aymptotic CLs 95% Cross Section Limit [fb]");
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->SetMarkerColor(kRed); gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
    gr->SetLineColor(kRed); gr->SetLineWidth(2);
    gr->SetFillStyle(3004); gr->SetFillColor(kRed);
    gr->Draw("ALPE3");

    TGraphErrors *gr2 = new TGraphErrors(6, x,relativeError2, ex,sigmaOfRelativeError2);
    gr2->SetTitle(""); gr2->SetMinimum(0.0); gr2->SetMaximum(200.);
    gr2->GetXaxis()->SetTitle("Photon Fake Rate In Endcap");
    gr2->GetYaxis()->SetTitle("Aymptotic CLs 95% Cross Section Limit [fb]");
    gr2->GetYaxis()->SetTitleSize(0.04);
    gr2->SetMarkerColor(kGreen+2); gr2->SetMarkerStyle(20); gr2->SetMarkerSize(1.2);
    gr2->SetLineColor(kGreen+2); gr2->SetLineWidth(2);
    gr2->SetFillStyle(3004); gr2->SetFillColor(kGreen+2);
    gr2->Draw("LPE3,same");

    TLegend *legend = new TLegend(0.20,0.70,0.45,0.85);
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
    tex->DrawLatex(0.60, 0.85, "Nominal Endcap");
    tex->DrawLatex(0.60, 0.80, "Fake Rate");
    tex->Draw();
    cv->Update();
    line = new TLine(0.003,0,0.003,200.);
    line->SetLineWidth(3); line->SetLineColor(kBlue);
    line->Draw();

    TLine *line3 = new TLine(0,33.9,0.015,33.9);
    line3->SetLineWidth(4); line3->SetLineColor(kBlack);
    line3->Draw();
    TLatex *tex3 = new TLatex();
    tex3->SetNDC();
    tex3->SetTextSize(0.05);
    tex3->SetTextFont(42);
    tex3->SetTextColor(kBlack);
    tex3->DrawLatex(0.67, 0.30, "SM Cross Section");
    tex3->Draw();

   cv->SetLogx(true);
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSLimitVsPhotonFakerate.gif");
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSLimitVsPhotonFakerate.pdf");
  }
   else if (scanOption == "phoEff") {
    //diBjet points for graph
    Float_t x[6] = {-10, 0, 10, 20, 30, 40};
    Float_t ex[6] = {0};
    
    //Draw diBjet Resolution Steps
    cv = new TCanvas("cv","cv", 800,600);
    gr = new TGraphErrors(6, x,relativeError1, ex,sigmaOfRelativeError1);
    gr->SetTitle(""); gr->SetMinimum(0.0); gr->SetMaximum(150.);
    gr->GetXaxis()->SetTitle("Relative Improvement In Photon Selection Efficiency [%]");
    gr->GetYaxis()->SetTitle("Aymptotic CLs 95% Cross Section Limit [fb]");
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->SetMarkerColor(kRed); gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
    gr->SetLineColor(kRed); gr->SetLineWidth(2);
    gr->SetFillStyle(3004); gr->SetFillColor(kRed);
    gr->Draw("ALPE3");


    TLegend *legend = new TLegend(0.50,0.65,0.75,0.80);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(gr,"All Other Improvements Applied", "LP");
    legend->Draw();


    TLine *line2 = new TLine(-15,33.9,45,33.9);
    line2->SetLineWidth(4); line2->SetLineColor(kBlack);
    line2->Draw();
    TLatex *tex2 = new TLatex();
    tex2->SetNDC();
    tex2->SetTextSize(0.05);
    tex2->SetTextFont(42);
    tex2->SetTextColor(kBlack);
    tex2->DrawLatex(0.67, 0.35, "SM Cross Section");
    tex2->Draw();

    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSLimitVsPhotonEffRatio.gif");
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSLimitVsPhotonEffRatio.pdf");
  }
   else if (scanOption == "btagEff") {
    //diBjet points for graph
    Float_t x[6] = {-10, 0, 10, 20, 30, 40};
    Float_t ex[6] = {0};
    
    //Draw diBjet Resolution Steps
    cv = new TCanvas("cv","cv", 800,600);
    gr = new TGraphErrors(6, x,relativeError1, ex,sigmaOfRelativeError1);
    gr->SetTitle(""); gr->SetMinimum(0.0); gr->SetMaximum(150.);
    gr->GetXaxis()->SetTitle("Relative Improvement In B-Tagging Efficiency [%]");
    gr->GetYaxis()->SetTitle("Aymptotic CLs 95% Cross Section Limit [fb]");
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->SetMarkerColor(kRed); gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
    gr->SetLineColor(kRed); gr->SetLineWidth(2);
    gr->SetFillStyle(3004); gr->SetFillColor(kRed);
    gr->Draw("ALPE3");


    TLegend *legend = new TLegend(0.50,0.65,0.75,0.80);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(gr,"All Other Improvements Applied", "LP");
    legend->Draw();


    TLine *line2 = new TLine(-15,33.9,45,33.9);
    line2->SetLineWidth(4); line2->SetLineColor(kBlack);
    line2->Draw();
    TLatex *tex2 = new TLatex();
    tex2->SetNDC();
    tex2->SetTextSize(0.05);
    tex2->SetTextFont(42);
    tex2->SetTextColor(kBlack);
    tex2->DrawLatex(0.67, 0.35, "SM Cross Section");
    tex2->Draw();

    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSLimitVsBtagEffRatio.gif");
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSLimitVsBtagEffRatio.pdf");
  }
  
   else if (scanOption == "lum") {
    //Luminosity points for graph
    Float_t x[13] = { 1.5, 2., 2.5, 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.};
    Float_t ex[13] = {0};
    
    //Draw Luminosity Steps
    cv = new TCanvas("cv","cv", 800,600);
    cv->SetLogy();

    gr = new TGraphErrors(13, x,relativeError1, ex,sigmaOfRelativeError1);
    gr->SetTitle(""); gr->SetMinimum(0.0); gr->SetMaximum(1000.);
    gr->GetXaxis()->SetTitle("Integrated Luminosity [10^{3} fb^{-1}]");
    gr->GetYaxis()->SetTitle("Aymptotic CLs 95% Cross Section Limit [fb]");
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->SetMarkerColor(kRed); gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
    gr->SetLineColor(kRed); gr->SetLineWidth(2);
    gr->SetFillStyle(3004); gr->SetFillColor(kRed);
    gr->Draw("ALPE3");

    TGraphErrors *gr2 = new TGraphErrors(13, x,relativeError2, ex,sigmaOfRelativeError2);
    gr2->SetTitle(""); gr2->SetMinimum(0.0); gr2->SetMaximum(1000.);
    gr2->GetXaxis()->SetTitle("M_{bb} Resolution Parameter [GeV/c^{2}]");
    gr2->GetYaxis()->SetTitle("Aymptotic CLs 95% Cross Section Limit [fb]");
    gr2->GetYaxis()->SetTitleSize(0.04);
    gr2->SetMarkerColor(kGreen+2); gr2->SetMarkerStyle(20); gr2->SetMarkerSize(1.2);
    gr2->SetLineColor(kGreen+2); gr2->SetLineWidth(2);
    gr2->SetFillStyle(3004); gr2->SetFillColor(kGreen+2);
    gr2->Draw("LPE3,same");

    TLegend *legend = new TLegend(0.55,0.70,0.80,0.90);
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
    tex->DrawLatex(0.33, 0.20, "Nominal Luminosity");
    tex->Draw();
    cv->Update();
    line = new TLine(3.,0,3.,1000.);
    line->SetLineWidth(3); line->SetLineColor(kBlue);
    line->Draw();

    TLine *line2 = new TLine(0.45,33.9,13.05,33.9);
    line2->SetLineWidth(4); line2->SetLineColor(kBlack);
    line2->Draw();
    TLatex *tex2 = new TLatex();
    tex2->SetNDC();
    tex2->SetTextSize(0.05);
    tex2->SetTextFont(42);
    tex2->SetTextColor(kBlack);
    tex2->DrawLatex(0.67, 0.57, "SM Cross Section");
    tex2->Draw();

    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSLimitVsLumi.gif");
    cv->SaveAs("Plots/AllSignalBkgd/Fits/PullPlots/XSLimitVsLumi.pdf");
  }

}




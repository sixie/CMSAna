//
//
//
// root -l CMSAna/macros/OverlayEfficiencyPlots.C+'(0 or 4 or 5,"energyRange")'
//
//
// maybe change location of legend
//
//

#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <vector>                   // STL vector class
#include <map>
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <TH1F.h>                
#include <TCanvas.h>                
#include <TGraphAsymmErrors.h>                
#include <THStack.h>
#include <TLegend.h>



TGraphAsymmErrors* loadGraph( string filename, string graphname) {
  
  TFile *file = new TFile(filename.c_str(), "READ");
  TGraphAsymmErrors* graphclone = (TGraphAsymmErrors*)( (TGraphAsymmErrors*)file->Get(graphname.c_str()) )->Clone(graphname.c_str());
  file->Close(); delete file;
  return graphclone;
  
}

void PlotBTagEfficiency() {

  Int_t colors[6] = { kBlue, kRed, kMagenta, kCyan, kGreen+2 , kBlack };

  TGraphAsymmErrors* BTagEff_type5_Pt_ttbar = loadGraph ( "BJetMistagRate_ttbar_type5.root", "MistagRate_CSVMedium_Pt");
  TGraphAsymmErrors* BTagEff_type5_Pt_qcd = loadGraph ( "BJetMistagRate_qcdPt80To120_type5.root", "MistagRate_CSVMedium_Pt");
  TGraphAsymmErrors* BTagEff_type5_Pt_diphotonjets = loadGraph ( "BJetMistagRate_DiphotonJets_type5.root", "MistagRate_CSVMedium_Pt");
  TGraphAsymmErrors* BTagEff_type5_Pt_diphotonjetsSherpa = loadGraph ( "BJetMistagRate_DiphotonJetsSherpa_type5.root", "MistagRate_CSVMedium_Pt");


  TGraphAsymmErrors* BTagEff_type4_Pt_ttbar = loadGraph ( "BJetMistagRate_ttbar_type4.root", "MistagRate_CSVMedium_Pt");
  TGraphAsymmErrors* BTagEff_type4_Pt_qcd = loadGraph ( "BJetMistagRate_qcdPt80To120_type4.root", "MistagRate_CSVMedium_Pt");
  TGraphAsymmErrors* BTagEff_type4_Pt_diphotonjets = loadGraph ( "BJetMistagRate_DiphotonJets_type4.root", "MistagRate_CSVMedium_Pt");
  TGraphAsymmErrors* BTagEff_type4_Pt_diphotonjetsSherpa = loadGraph ( "BJetMistagRate_DiphotonJetsSherpa_type4.root", "MistagRate_CSVMedium_Pt");


  TGraphAsymmErrors* BTagEff_type0_Pt_ttbar = loadGraph ( "BJetMistagRate_ttbar_type0.root", "MistagRate_CSVMedium_Pt");
  TGraphAsymmErrors* BTagEff_type0_Pt_qcd = loadGraph ( "BJetMistagRate_qcdPt80To120_type0.root", "MistagRate_CSVMedium_Pt");
  TGraphAsymmErrors* BTagEff_type0_Pt_diphotonjets = loadGraph ( "BJetMistagRate_DiphotonJets_type0.root", "MistagRate_CSVMedium_Pt");
  TGraphAsymmErrors* BTagEff_type0_Pt_diphotonjetsSherpa = loadGraph ( "BJetMistagRate_DiphotonJetsSherpa_type0.root", "MistagRate_CSVMedium_Pt");




  TCanvas *cv = 0;
  TLegend *legend = 0;




  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.2,0.7,0.5,0.9);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(BTagEff_type5_Pt_ttbar, "ttbar (madgraph)", "LP");
  legend->AddEntry(BTagEff_type5_Pt_qcd, "qcd 80-120 (pythia)", "LP");
  legend->AddEntry(BTagEff_type5_Pt_diphotonjets, "Diphoton+jets (Madgraph)", "LP");
  legend->AddEntry(BTagEff_type5_Pt_diphotonjetsSherpa, "Diphoton+jets (Sherpa)", "LP");
  BTagEff_type5_Pt_ttbar->GetXaxis()->SetRangeUser(30,100);
  BTagEff_type5_Pt_ttbar->GetYaxis()->SetRangeUser(0,1.0);
  BTagEff_type5_Pt_ttbar->GetXaxis()->SetTitle("jet p_{T} [GeV/c]");

  BTagEff_type5_Pt_ttbar->SetLineColor(colors[0]);
  BTagEff_type5_Pt_qcd->SetLineColor(colors[1]);
  BTagEff_type5_Pt_diphotonjets->SetLineColor(colors[2]);
  BTagEff_type5_Pt_diphotonjetsSherpa->SetLineColor(colors[3]);
  BTagEff_type5_Pt_ttbar->SetMarkerColor(colors[0]);
  BTagEff_type5_Pt_qcd->SetMarkerColor(colors[1]);
  BTagEff_type5_Pt_diphotonjets->SetMarkerColor(colors[2]);
  BTagEff_type5_Pt_diphotonjetsSherpa->SetMarkerColor(colors[3]);

  BTagEff_type5_Pt_ttbar->Draw("AP");
  BTagEff_type5_Pt_qcd->Draw("P");
  BTagEff_type5_Pt_diphotonjets->Draw("P");
  BTagEff_type5_Pt_diphotonjetsSherpa->Draw("P");
  legend->Draw();

  cv->SaveAs("BTagEfficiency_bquark.gif");
    



  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.2,0.7,0.5,0.9);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(BTagEff_type4_Pt_ttbar, "ttbar (madgraph)", "LP");
  legend->AddEntry(BTagEff_type4_Pt_qcd, "qcd 80-120 (pythia)", "LP");
  legend->AddEntry(BTagEff_type4_Pt_diphotonjets, "Diphoton+jets (Madgraph)", "LP");
  legend->AddEntry(BTagEff_type4_Pt_diphotonjetsSherpa, "Diphoton+jets (Sherpa)", "LP");
  BTagEff_type4_Pt_ttbar->GetXaxis()->SetRangeUser(30,100);
  BTagEff_type4_Pt_ttbar->GetYaxis()->SetRangeUser(0,0.3);
  BTagEff_type4_Pt_ttbar->GetXaxis()->SetTitle("jet p_{T} [GeV/c]");

  BTagEff_type4_Pt_ttbar->SetLineColor(colors[0]);
  BTagEff_type4_Pt_qcd->SetLineColor(colors[1]);
  BTagEff_type4_Pt_diphotonjets->SetLineColor(colors[2]);
  BTagEff_type4_Pt_diphotonjetsSherpa->SetLineColor(colors[3]);
  BTagEff_type4_Pt_ttbar->SetMarkerColor(colors[0]);
  BTagEff_type4_Pt_qcd->SetMarkerColor(colors[1]);
  BTagEff_type4_Pt_diphotonjets->SetMarkerColor(colors[2]);
  BTagEff_type4_Pt_diphotonjetsSherpa->SetMarkerColor(colors[3]);

  BTagEff_type4_Pt_ttbar->Draw("AP");
  BTagEff_type4_Pt_qcd->Draw("P");
  BTagEff_type4_Pt_diphotonjets->Draw("P");
  BTagEff_type4_Pt_diphotonjetsSherpa->Draw("P");
  legend->Draw();

  cv->SaveAs("BTagEfficiency_cquark.gif");
    


  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.2,0.7,0.5,0.9);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(BTagEff_type0_Pt_ttbar, "ttbar (madgraph)", "LP");
  legend->AddEntry(BTagEff_type0_Pt_qcd, "qcd 80-120 (pythia)", "LP");
  legend->AddEntry(BTagEff_type0_Pt_diphotonjets, "Diphoton+jets (Madgraph)", "LP");
  legend->AddEntry(BTagEff_type0_Pt_diphotonjetsSherpa, "Diphoton+jets (Sherpa)", "LP");
  BTagEff_type0_Pt_ttbar->GetXaxis()->SetRangeUser(30,100);
  BTagEff_type0_Pt_ttbar->GetYaxis()->SetRangeUser(0,0.05);
  BTagEff_type0_Pt_ttbar->GetXaxis()->SetTitle("jet p_{T} [GeV/c]");

  BTagEff_type0_Pt_ttbar->SetLineColor(colors[0]);
  BTagEff_type0_Pt_qcd->SetLineColor(colors[1]);
  BTagEff_type0_Pt_diphotonjets->SetLineColor(colors[2]);
  BTagEff_type0_Pt_diphotonjetsSherpa->SetLineColor(colors[3]);
  BTagEff_type0_Pt_ttbar->SetMarkerColor(colors[0]);
  BTagEff_type0_Pt_qcd->SetMarkerColor(colors[1]);
  BTagEff_type0_Pt_diphotonjets->SetMarkerColor(colors[2]);
  BTagEff_type0_Pt_diphotonjetsSherpa->SetMarkerColor(colors[3]);

  BTagEff_type0_Pt_ttbar->Draw("AP");
  BTagEff_type0_Pt_qcd->Draw("P");
  BTagEff_type0_Pt_diphotonjets->Draw("P");
  BTagEff_type0_Pt_diphotonjetsSherpa->Draw("P");
  legend->Draw();

  cv->SaveAs("BTagEfficiency_mistags.gif");
    


}

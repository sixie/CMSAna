//================================================================================================
//
// Simple Example
//
//root -l CMSAna/ObjectStudies/BTagging/ComputeBTaggingEfficiency.C+'("/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/BJetEfficiencyRate/ntuples/BJetEfficiencyNtuple.diphjets-START53_V7A.root",0,-1,"type0")'
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <TH1F.h>                
#include <TCanvas.h>                
#include <TGraphAsymmErrors.h>                

#include "CMSAna/Utils/EfficiencyUtils.hh"
#include "CMSAna/ObjectStudies/interface/EfficiencyTree.h"

#endif

bool passSelection( double pt, double eta, int option) {

    bool passSelection = false;
    if (option == -1) {
      passSelection = true;
    }
    else if (option == 1) {
      if (pt > 20 && pt <= 30) passSelection = true;
    }
    else if (option == 2) {
      if (pt > 30 && pt <= 40) passSelection = true;
    } 
    else if (option == 3) {
      if (pt > 40 && pt <= 50) passSelection = true;
    }
    else if (option == 4) {
      if (pt > 50 && pt <= 60) passSelection = true;
    } 
    else if (option == 5) {
      if (pt > 60 && pt <= 70) passSelection = true;
    }
    else if (option == 6) {
      if (pt > 70 && pt <= 80) passSelection = true;
    }
    else if (option == 7) {
      if (pt > 80 && pt <= 90) passSelection = true;
    } 
    else if (option == 8) {
      if (pt > 90 && pt <= 100) passSelection = true;
    }

    return passSelection;

}



//=== MAIN MACRO ================================================================================================= 

void ComputeBTaggingEfficiency(const string inputfile, int type = 0, int option = -1, string label = "") {
  
  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;


  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  TH1F *histDenominatorPt = new TH1F ("histDenominatorPt",";Electron p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorPt = new TH1F ("histNumeratorPt",";Electron p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorEta = new TH1F ("histDenominatorEta",";Electron Eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histNumeratorEta = new TH1F ("histNumeratorEta",";Electron Eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histDenominatorPhi = new TH1F ("histDenominatorPhi",";Electron Phi; Number of Events", 50, 0 , 3.2);
  TH1F *histNumeratorPhi = new TH1F ("histNumeratorPhi",";Electron Phi; Number of Events", 50, 0 , 3.2);
  TH1F *histDenominatorRho = new TH1F ("histDenominatorRho",";Electron Rho; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorRho = new TH1F ("histNumeratorRho",";Electron Rho; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpv = new TH1F ("histDenominatorNpv",";Electron Npv; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpv = new TH1F ("histNumeratorNpv",";Electron Npv; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpu = new TH1F ("histDenominatorNpu",";Electron Npu; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpu = new TH1F ("histNumeratorNpu",";Electron Npu; Number of Events", 50, 0 , 100);

  TH2F *histDenominatorPtEta = new TH2F ("histDenominatorPtEta",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", 50, 0 , 200, 50, -3.0, 3.0);
  TH2F *histNumeratorPtEta = new TH2F ("histNumeratorPtEta",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", 50, 0 , 200, 50, -3.0, 3.0);

  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  cmsana::EfficiencyTree efftree;
  efftree.LoadTree(inputfile.c_str());
  efftree.InitTree();

  for (int n=0;n<efftree.tree_->GetEntries();n++) { 
    efftree.tree_->GetEntry(n);

    if (type == 0) {
      if (!(abs(efftree.matchedPdgId) != 5 && abs(efftree.matchedPdgId) != 4 )) continue;
    } else if ( (type >= 1 && type <= 5) || type == 21) {
      if (!(abs(efftree.matchedPdgId) == type)) continue;
    }

    if (!passSelection(efftree.pt, efftree.eta, option)) continue;

    //**** PT - ETA ****
    histDenominatorPtEta->Fill(efftree.pt,efftree.eta);
    if(efftree.pass) {
      histNumeratorPtEta->Fill(efftree.pt,efftree.eta);
    }

    //**** PT ****
    if (fabs(efftree.eta) < 2.4) {
      histDenominatorPt->Fill(efftree.pt);

      //Numerator
      if(efftree.pass) {
        histNumeratorPt->Fill(efftree.pt);        
      }

    }

    //**** Eta ****
    if (fabs(efftree.pt) > 30) {
      histDenominatorEta->Fill(efftree.eta);

      //Numerator
      if(efftree.pass) {
        histNumeratorEta->Fill(efftree.eta);        
      }

    }

    //**** Phi ****
    if (fabs(efftree.eta) < 2.4) {
      histDenominatorPhi->Fill(efftree.phi);

      //Numerator
      if(efftree.pass) {
        histNumeratorPhi->Fill(efftree.phi);        
      }

    }

    //**** Rho ****
    if (fabs(efftree.eta) < 2.4) {
      histDenominatorRho->Fill(efftree.rho);

      //Numerator
      if(efftree.pass) {
        histNumeratorRho->Fill(efftree.rho);        
      }

    }
    //**** Npv ****
    if (fabs(efftree.eta) < 2.4) {
      histDenominatorNpv->Fill(efftree.npv);

      //Numerator
      if(efftree.pass) {
        histNumeratorNpv->Fill(efftree.npv);        
      }

    }
    //**** Npu ****
    if (fabs(efftree.eta) < 2.4) {
      histDenominatorNpu->Fill(efftree.npu);

      //Numerator
      if(efftree.pass) {
        histNumeratorNpu->Fill(efftree.npu);        
      }

    }

  }


  //--------------------------------------------------------------------------------------------------------------
  // Make Efficiency Plots
  //==============================================================================================================

  TGraphAsymmErrors *efficiency_pt = cmsana::createEfficiencyGraph(histNumeratorPt, histDenominatorPt, "Efficiency_Pt" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_eta = cmsana::createEfficiencyGraph(histNumeratorEta, histDenominatorEta, "Efficiency_Eta" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_phi = cmsana::createEfficiencyGraph(histNumeratorPhi, histDenominatorPhi, "Efficiency_Phi" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_rho = cmsana::createEfficiencyGraph(histNumeratorRho, histDenominatorRho, "Efficiency_Rho" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_npv = cmsana::createEfficiencyGraph(histNumeratorNpv, histDenominatorNpv, "Efficiency_Npv" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_npu = cmsana::createEfficiencyGraph(histNumeratorNpu, histDenominatorNpu, "Efficiency_Npu" , vector<double>() ,  -99, -99, 0, 1);  
  TH2F *efficiency_pteta = cmsana::createEfficiencyHist2D(histNumeratorPtEta, histDenominatorPtEta, "Efficiency_PtEta" , vector<double>() ,vector<double>());  


  //--------------------------------------------------------------------------------------------------------------
  // Draw
  //==============================================================================================================
  TCanvas *cv =0;

  cv = new TCanvas("cv","cv",800,600);
  efficiency_pt->Draw("AP");
  efficiency_pt->SetTitle("");
  efficiency_pt->GetYaxis()->SetRangeUser(0.0,0.1);
  cv->SaveAs("Efficiency_Pt.gif");

  cv = new TCanvas("cv","cv",800,600);
  efficiency_eta->Draw("AP");
  efficiency_eta->SetTitle("");
  efficiency_eta->GetYaxis()->SetRangeUser(0.0,0.1);
  cv->SaveAs("Efficiency_Eta.gif");

  cv = new TCanvas("cv","cv",800,600);
  efficiency_phi->Draw("AP");
  efficiency_phi->SetTitle("");
  efficiency_phi->GetYaxis()->SetRangeUser(0.0,0.1);
  cv->SaveAs("Efficiency_Phi.gif");

  cv = new TCanvas("cv","cv",800,600);
  efficiency_rho->Draw("AP");
  efficiency_rho->SetTitle("");
  efficiency_rho->GetYaxis()->SetRangeUser(0.0,0.1);
  cv->SaveAs("Efficiency_Rho.gif");

  cv = new TCanvas("cv","cv",800,600);
  efficiency_npv->Draw("AP");
  efficiency_npv->SetTitle("");
  efficiency_npv->GetYaxis()->SetRangeUser(0.0,0.1);
  cv->SaveAs("Efficiency_Npv.gif");

  cv = new TCanvas("cv","cv",800,600);
  efficiency_npu->Draw("AP");
  efficiency_npu->SetTitle("");
  efficiency_npu->GetYaxis()->SetRangeUser(0.0,0.1);
  cv->SaveAs("Efficiency_Npu.gif");


  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("Efficiency"+Label+".root").c_str(), "UPDATE");
  file->cd();
  file->WriteTObject(efficiency_pt, efficiency_pt->GetName(), "WriteDelete");
  file->WriteTObject(efficiency_eta, efficiency_eta->GetName(), "WriteDelete");
  file->WriteTObject(efficiency_phi, efficiency_phi->GetName(), "WriteDelete");
  file->WriteTObject(efficiency_rho, efficiency_rho->GetName(), "WriteDelete");
  file->WriteTObject(efficiency_npv, efficiency_npv->GetName(), "WriteDelete");
  file->WriteTObject(efficiency_npu, efficiency_npu->GetName(), "WriteDelete");
  file->WriteTObject(efficiency_pteta, efficiency_pteta->GetName(), "WriteDelete");

  file->Close();
  delete file;       

}

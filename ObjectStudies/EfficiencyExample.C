//================================================================================================
//
// Simple Example
//
//root -l CMSAna/ObjectStudies/EfficiencyExample.C+'("/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/BJetMistagRate/ntuples/BJetMistagNtuple.qcd_pt30To50-START53_V7A.root")'
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

void EfficiencyExample(const string inputfile, int option = -1) {
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;


  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  TH1F *histDenominatorPt = new TH1F ("histDenominatorPt","; p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorPt = new TH1F ("histNumeratorPt","; p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorEta = new TH1F ("histDenominatorEta","; Eta; Number of Events", 50, -3 , 3);
  TH1F *histNumeratorEta = new TH1F ("histNumeratorEta","; Eta; Number of Events", 50, -3 , -3);
  TH1F *histDenominatorPhi = new TH1F ("histDenominatorPhi","; Phi; Number of Events", 50, 0 , 3.2);
  TH1F *histNumeratorPhi = new TH1F ("histNumeratorPhi","; Phi; Number of Events", 50, 0 , 3.2);
  TH1F *histDenominatorRho = new TH1F ("histDenominatorRho","; Rho; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorRho = new TH1F ("histNumeratorRho","; Rho; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpv = new TH1F ("histDenominatorNpv","; Number of Primary Vertices; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpv = new TH1F ("histNumeratorNpv","; Number of Primary Vertices; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpu = new TH1F ("histDenominatorNpu","; Number of Pileup Events; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpu = new TH1F ("histNumeratorNpu","; Number of Pileup Events; Number of Events", 50, 0 , 100);

  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  cmsana::EfficiencyTree efftree;
  efftree.LoadTree(inputfile.c_str());
  efftree.InitTree();

  for (int n=0;n<efftree.tree_->GetEntries();n++) { 
    efftree.tree_->GetEntry(n);

    if (!passSelection(efftree.pt, efftree.eta, option)) continue;

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

  TGraphAsymmErrors *efficiency_pt = cmsana::createEfficiencyGraph(histNumeratorPt, histDenominatorPt, "MistagRate_CSVMedium_Pt" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_eta = cmsana::createEfficiencyGraph(histNumeratorEta, histDenominatorEta, "MistagRate_CSVMedium_Eta" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_phi = cmsana::createEfficiencyGraph(histNumeratorPhi, histDenominatorPhi, "MistagRate_CSVMedium_Phi" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_rho = cmsana::createEfficiencyGraph(histNumeratorRho, histDenominatorRho, "MistagRate_CSVMedium_Rho" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_npv = cmsana::createEfficiencyGraph(histNumeratorNpv, histDenominatorNpv, "MistagRate_CSVMedium_Npv" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_npu = cmsana::createEfficiencyGraph(histNumeratorNpu, histDenominatorNpu, "MistagRate_CSVMedium_Npu" , vector<double>() ,  -99, -99, 0, 1);  


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
  TFile *file = TFile::Open("MistagRate.root", "UPDATE");
  file->cd();
  file->WriteTObject(efficiency_pt, efficiency_pt->GetName(), "WriteDelete");
  file->Close();
  delete file;       
  
}

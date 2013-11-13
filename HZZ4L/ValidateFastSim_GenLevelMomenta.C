// root -l HiggsAna/HZZ4l/HLL/ValidateFastSimExample.C+'("HZZEventNtuple_ZZ_0000.root","HZZEventNtuple_FastSim_ZZ_0000.root","ZZ")'

//================================================================================================
//
//
//
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <iostream>                 // standard I/O
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TPad.h"

#include "CMSAna/HZZ4L/interface/HZZEventTree.h"
#include "CMSAna/Utils/CommonDefs.hh"
#include "CMSAna/DataTree/interface/Types.h"

// #include "TLorentzVector.h"

#endif

//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
void NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return;
}

//*************************************************************************************************
//Divide Histograms
//*************************************************************************************************
TH1F* MakeRelative(TH1F *hist, TH1F *reference) {

  TH1F *relative = (TH1F*)hist->Clone((string(hist->GetName())+"_relative").c_str());
  relative->SetTitle("");
  assert(hist->GetXaxis()->GetNbins() == reference->GetXaxis()->GetNbins());

  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    if (reference->GetBinContent(b) > 0) {
      relative->SetBinContent(b,100 * ( hist->GetBinContent(b) - reference->GetBinContent(b) )/reference->GetBinContent(b));
    } else {
      relative->SetBinContent(b,0);
    }
  }

  return relative;
}

void DrawComparison( TH1F* hist1, TH1F* hist2, string label1, string label2, string plotname, bool doNormalize = false) {

  TH1F* relative = MakeRelative(hist2,hist1);

  TCanvas *cv = new TCanvas("cv","cv", 800,600);

  TPad *pad1 = new TPad("pad1","pad1", 0,0.2,1,1);
  pad1->SetBottomMargin(0.125);
  pad1->Draw();
  pad1->cd();

  TLegend *legend = new TLegend(0.54,0.14,0.94,0.44);
  legend = new TLegend(0.74,0.70,0.94,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(hist1,label1.c_str(), "LP");
  legend->AddEntry(hist2,label2.c_str(), "LP");

  if (doNormalize) {
    NormalizeHist(hist1);
    NormalizeHist(hist2);
  }
  hist1->SetMaximum(1.2*fmax(hist1->GetMaximum(), hist2->GetMaximum()));

  hist1->SetLineColor(kRed);
  hist1->Draw("hist");
  hist2->SetLineColor(kBlue);
  hist2->Draw("hist,same");
  legend->Draw();

  cv->cd();
  TPad *pad2 = new TPad("pad2","pad2", 0,0,1,0.2);
  pad1->SetTopMargin(0.01);
  pad2->Draw();
  pad2->cd();

  relative->GetYaxis()->SetTitle("% Difference");
  relative->GetYaxis()->SetNdivisions(306);
  relative->GetYaxis()->SetTitleSize(0.15);
  relative->GetYaxis()->SetTitleOffset(0.3);
  relative->GetYaxis()->SetRangeUser(-10,10);
  relative->GetYaxis()->SetLabelSize(0.15);
  relative->GetXaxis()->SetLabelSize(0.0);
  relative->SetLineColor(kBlack);
  relative->SetMarkerColor(kBlack);
  relative->Draw("hist");

  cv->SaveAs(plotname.c_str());

}


//*************************************************************************************************
//Main Function
//*************************************************************************************************
void ValidateFastSim_GenLevelMomenta(string fullsimFilename, string fastsimFilename,
                           const string Label = "ZZ") {

  string label = Label;
  if (Label != "") label = "_" + Label;

  //*****************************************************************************************
  // Load Input
  //*****************************************************************************************
  HZZEventTree fullsimHZZEventTree;
  fullsimHZZEventTree.LoadTree(fullsimFilename.c_str());
  fullsimHZZEventTree.InitTree();
  HZZEventTree fastsimHZZEventTree;
  fastsimHZZEventTree.LoadTree(fastsimFilename.c_str());
  fastsimHZZEventTree.InitTree();
  cout << "Events in the ntuple: " << fastsimHZZEventTree.tree_->GetEntries() << endl;

  //*************************************************************************************************
  //Histograms
  //*************************************************************************************************
  TH1::AddDirectory(kFALSE);
  TH1F *histZ1Mass_ee_fullsim = new TH1F( "histZ1Mass_ee_fullsim", ";Z1 Mass [GeV/c^{2}]; Number of Events", 100, 0,200);
  TH1F *histZ1Mass_ee_fastsim = new TH1F( "histZ1Mass_ee_fastsim", ";Z1 Mass [GeV/c^{2}]; Number of Events", 100, 0,200);
  TH1F *histZ1Mass_mm_fullsim = new TH1F( "histZ1Mass_mm_fullsim", ";Z1 Mass [GeV/c^{2}]; Number of Events", 100, 0,200);
  TH1F *histZ1Mass_mm_fastsim = new TH1F( "histZ1Mass_mm_fastsim", ";Z1 Mass [GeV/c^{2}]; Number of Events", 100, 0,200);
  TH1F *histZ2Mass_ee_fullsim = new TH1F( "histZ2Mass_ee_fullsim", ";Z2 Mass [GeV/c^{2}]; Number of Events", 100, 0,200);
  TH1F *histZ2Mass_ee_fastsim = new TH1F( "histZ2Mass_ee_fastsim", ";Z2 Mass [GeV/c^{2}]; Number of Events", 100, 0,200);
  TH1F *histZ2Mass_mm_fullsim = new TH1F( "histZ2Mass_mm_fullsim", ";Z2 Mass [GeV/c^{2}]; Number of Events", 100, 0,200);
  TH1F *histZ2Mass_mm_fastsim = new TH1F( "histZ2Mass_mm_fastsim", ";Z2 Mass [GeV/c^{2}]; Number of Events", 100, 0,200);
  TH1F *histZZMass_4e_fullsim = new TH1F( "histZZMass_4e_fullsim", ";ZZ Mass [GeV/c^{2}]; Number of Events", 300, 50,350);
  TH1F *histZZMass_4e_fastsim = new TH1F( "histZZMass_4e_fastsim", ";ZZ Mass [GeV/c^{2}]; Number of Events", 300, 50,350);
  TH1F *histZZMass_4m_fullsim = new TH1F( "histZZMass_4m_fullsim", ";ZZ Mass [GeV/c^{2}]; Number of Events", 300, 50,350);
  TH1F *histZZMass_4m_fastsim = new TH1F( "histZZMass_4m_fastsim", ";ZZ Mass [GeV/c^{2}]; Number of Events", 300, 50,350);
  TH1F *histZZMass_2e2m_fullsim = new TH1F( "histZZMass_2e2m_fullsim", ";ZZ Mass [GeV/c^{2}]; Number of Events", 300, 50,350);
  TH1F *histZZMass_2e2m_fastsim = new TH1F( "histZZMass_2e2m_fastsim", ";ZZ Mass [GeV/c^{2}]; Number of Events", 300, 50,350);
  TH1F *histZZMass_2m2e_fullsim = new TH1F( "histZZMass_2m2e_fullsim", ";ZZ Mass [GeV/c^{2}]; Number of Events", 300, 50,350);
  TH1F *histZZMass_2m2e_fastsim = new TH1F( "histZZMass_2m2e_fastsim", ";ZZ Mass [GeV/c^{2}]; Number of Events", 300, 50,350);

  TH1F *histLep1Pt_e_fullsim = new TH1F( "histLep1Pt_e_fullsim", ";lep1 pT [GeV/c]; Number of Events", 100, 0,200);
  TH1F *histLep1Pt_e_fastsim = new TH1F( "histLep1Pt_e_fastsim", ";lep1 pT [GeV/c]; Number of Events", 100, 0,200);
  TH1F *histLep1Pt_m_fullsim = new TH1F( "histLep1Pt_m_fullsim", ";lep1 pT [GeV/c]; Number of Events", 100, 0,200);
  TH1F *histLep1Pt_m_fastsim = new TH1F( "histLep1Pt_m_fastsim", ";lep1 pT [GeV/c]; Number of Events", 100, 0,200);
  TH1F *histLep2Pt_e_fullsim = new TH1F( "histLep2Pt_e_fullsim", ";lep2 pT [GeV/c]; Number of Events", 100, 0,100);
  TH1F *histLep2Pt_e_fastsim = new TH1F( "histLep2Pt_e_fastsim", ";lep2 pT [GeV/c]; Number of Events", 100, 0,100);
  TH1F *histLep2Pt_m_fullsim = new TH1F( "histLep2Pt_m_fullsim", ";lep2 pT [GeV/c]; Number of Events", 100, 0,100);
  TH1F *histLep2Pt_m_fastsim = new TH1F( "histLep2Pt_m_fastsim", ";lep2 pT [GeV/c]; Number of Events", 100, 0,100);
  TH1F *histLep3Pt_e_fullsim = new TH1F( "histLep3Pt_e_fullsim", ";lep3 pT [GeV/c]; Number of Events", 100, 0,100);
  TH1F *histLep3Pt_e_fastsim = new TH1F( "histLep3Pt_e_fastsim", ";lep3 pT [GeV/c]; Number of Events", 100, 0,100);
  TH1F *histLep3Pt_m_fullsim = new TH1F( "histLep3Pt_m_fullsim", ";lep3 pT [GeV/c]; Number of Events", 100, 0,100);
  TH1F *histLep3Pt_m_fastsim = new TH1F( "histLep3Pt_m_fastsim", ";lep3 pT [GeV/c]; Number of Events", 100, 0,100);
  TH1F *histLep4Pt_e_fullsim = new TH1F( "histLep4Pt_e_fullsim", ";lep4 pT [GeV/c]; Number of Events", 100, 0,50);
  TH1F *histLep4Pt_e_fastsim = new TH1F( "histLep4Pt_e_fastsim", ";lep4 pT [GeV/c]; Number of Events", 100, 0,50);
  TH1F *histLep4Pt_m_fullsim = new TH1F( "histLep4Pt_m_fullsim", ";lep4 pT [GeV/c]; Number of Events", 100, 0,50);
  TH1F *histLep4Pt_m_fastsim = new TH1F( "histLep4Pt_m_fastsim", ";lep4 pT [GeV/c]; Number of Events", 100, 0,50);

  TH1F *histLep1Eta_e_fullsim = new TH1F( "histLep1Eta_e_fullsim", ";lep1 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep1Eta_e_fastsim = new TH1F( "histLep1Eta_e_fastsim", ";lep1 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep1Eta_m_fullsim = new TH1F( "histLep1Eta_m_fullsim", ";lep1 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep1Eta_m_fastsim = new TH1F( "histLep1Eta_m_fastsim", ";lep1 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep2Eta_e_fullsim = new TH1F( "histLep2Eta_e_fullsim", ";lep2 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep2Eta_e_fastsim = new TH1F( "histLep2Eta_e_fastsim", ";lep2 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep2Eta_m_fullsim = new TH1F( "histLep2Eta_m_fullsim", ";lep2 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep2Eta_m_fastsim = new TH1F( "histLep2Eta_m_fastsim", ";lep2 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep3Eta_e_fullsim = new TH1F( "histLep3Eta_e_fullsim", ";lep3 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep3Eta_e_fastsim = new TH1F( "histLep3Eta_e_fastsim", ";lep3 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep3Eta_m_fullsim = new TH1F( "histLep3Eta_m_fullsim", ";lep3 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep3Eta_m_fastsim = new TH1F( "histLep3Eta_m_fastsim", ";lep3 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep4Eta_e_fullsim = new TH1F( "histLep4Eta_e_fullsim", ";lep4 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep4Eta_e_fastsim = new TH1F( "histLep4Eta_e_fastsim", ";lep4 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep4Eta_m_fullsim = new TH1F( "histLep4Eta_m_fullsim", ";lep4 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep4Eta_m_fastsim = new TH1F( "histLep4Eta_m_fastsim", ";lep4 #eta; Number of Events", 100, -2.5,2.5);

  //*************************************************************************************************
  //Yields
  //*************************************************************************************************
  double fullsimYield = 0;
  double fastsimYield = 0;

  //*****************************************************************************************
  // Loop over fullsim events
  //*****************************************************************************************
  for(UInt_t ientry=0; ientry < fullsimHZZEventTree.tree_->GetEntries(); ientry++) {       	
    fullsimHZZEventTree.tree_->GetEntry(ientry);
    if (ientry % 10000 == 0) cout << "Event " << ientry << endl;

    //how to get lepton four momenta

    cmsana::FourVectorM lepton1FourVector;
    lepton1FourVector.SetPt(fullsimHZZEventTree.l1pt);
    lepton1FourVector.SetEta(fullsimHZZEventTree.l1eta);
    lepton1FourVector.SetPhi(fullsimHZZEventTree.l1phi);
    lepton1FourVector.SetM( (abs(fullsimHZZEventTree.l1id) == 11) ? ELECTRONMASS : MUONMASS);

    
    //*****************************************************************************************
    // leading and subleading lepton
    //*****************************************************************************************
    double leadingLepPt = fullsimHZZEventTree.l1pt;
    double subleadingLepPt = fullsimHZZEventTree.l2pt;
    if (fullsimHZZEventTree.l2pt > fullsimHZZEventTree.l1pt) {
      leadingLepPt = fullsimHZZEventTree.l2pt;
      subleadingLepPt = fullsimHZZEventTree.l1pt;
    }
    if (fullsimHZZEventTree.l3pt > leadingLepPt) {
      subleadingLepPt = leadingLepPt;
      leadingLepPt = fullsimHZZEventTree.l3pt;      
    } else if (fullsimHZZEventTree.l3pt > subleadingLepPt) {
      subleadingLepPt = fullsimHZZEventTree.l3pt;
    }
    if (fullsimHZZEventTree.l4pt > leadingLepPt) {
      subleadingLepPt = leadingLepPt;
      leadingLepPt = fullsimHZZEventTree.l4pt;      
    } else if (fullsimHZZEventTree.l4pt > subleadingLepPt) {
      subleadingLepPt = fullsimHZZEventTree.l4pt;
    }


    //fill events passing selection
    if (fullsimHZZEventTree.l1pt > (abs(fullsimHZZEventTree.l1id) == 11 ? 7 : 5)
        && fullsimHZZEventTree.l2pt > (abs(fullsimHZZEventTree.l2id) == 11 ? 7 : 5)
        && fullsimHZZEventTree.l3pt > (abs(fullsimHZZEventTree.l3id) == 11 ? 7 : 5)
        && fullsimHZZEventTree.l4pt > (abs(fullsimHZZEventTree.l4id) == 11 ? 7 : 5)
        && abs(fullsimHZZEventTree.l1eta) < (abs(fullsimHZZEventTree.l1id) == 11 ? 2.5 : 2.4)
        && abs(fullsimHZZEventTree.l2eta) < (abs(fullsimHZZEventTree.l2id) == 11 ? 2.5 : 2.4)
        && abs(fullsimHZZEventTree.l3eta) < (abs(fullsimHZZEventTree.l3id) == 11 ? 2.5 : 2.4)
        && abs(fullsimHZZEventTree.l4eta) < (abs(fullsimHZZEventTree.l4id) == 11 ? 2.5 : 2.4)
        && leadingLepPt > 20
        && subleadingLepPt > 10
      ) {

      fullsimYield++;

      //******************************************
      //ZMass
      //******************************************
      if (abs(fullsimHZZEventTree.l1id) == 11) {
        histZ1Mass_ee_fullsim->Fill(fullsimHZZEventTree.genz1mass, fullsimHZZEventTree.weight);
      } else {
        histZ1Mass_mm_fullsim->Fill(fullsimHZZEventTree.genz1mass, fullsimHZZEventTree.weight);
      }
      if (abs(fullsimHZZEventTree.l3id) == 11) {
        histZ2Mass_ee_fullsim->Fill(fullsimHZZEventTree.genz2mass, fullsimHZZEventTree.weight);
      } else {
        histZ2Mass_mm_fullsim->Fill(fullsimHZZEventTree.genz2mass, fullsimHZZEventTree.weight);
      }
      if (fullsimHZZEventTree.channel == 0) histZZMass_4e_fullsim->Fill(fullsimHZZEventTree.genzzmass, fullsimHZZEventTree.weight);
      if (fullsimHZZEventTree.channel == 1) histZZMass_4m_fullsim->Fill(fullsimHZZEventTree.genzzmass, fullsimHZZEventTree.weight);
      if (fullsimHZZEventTree.channel == 2) histZZMass_2e2m_fullsim->Fill(fullsimHZZEventTree.genzzmass, fullsimHZZEventTree.weight);
      if (fullsimHZZEventTree.channel == 3) histZZMass_2m2e_fullsim->Fill(fullsimHZZEventTree.genzzmass, fullsimHZZEventTree.weight);

      //******************************************
      //Lepton Pt, Eta
      //******************************************
      if (abs(fullsimHZZEventTree.genl1id) == 11) {
        if (fullsimHZZEventTree.l1pt > 7) {
          histLep1Pt_e_fullsim->Fill(fullsimHZZEventTree.genl1pt, fullsimHZZEventTree.weight);
          histLep1Eta_e_fullsim->Fill(fullsimHZZEventTree.genl1eta, fullsimHZZEventTree.weight);
        }
      } else {
        if (fullsimHZZEventTree.genl1pt > 5) {
          histLep1Pt_m_fullsim->Fill(fullsimHZZEventTree.genl1pt, fullsimHZZEventTree.weight);
          histLep1Eta_m_fullsim->Fill(fullsimHZZEventTree.genl1eta, fullsimHZZEventTree.weight);
        }
      }
      if (abs(fullsimHZZEventTree.genl2id) == 11) {
        if (fullsimHZZEventTree.genl2pt > 7) {
          histLep2Pt_e_fullsim->Fill(fullsimHZZEventTree.genl2pt, fullsimHZZEventTree.weight);
          histLep2Eta_e_fullsim->Fill(fullsimHZZEventTree.genl2eta, fullsimHZZEventTree.weight);
        }
      } else {
        if (fullsimHZZEventTree.genl2pt > 5) {
          histLep2Pt_m_fullsim->Fill(fullsimHZZEventTree.genl2pt, fullsimHZZEventTree.weight);
          histLep2Eta_m_fullsim->Fill(fullsimHZZEventTree.genl2eta, fullsimHZZEventTree.weight);
        }
      }

      if (abs(fullsimHZZEventTree.genl3id) == 11) {
        if (fullsimHZZEventTree.genl3pt > 7) {
          histLep3Pt_e_fullsim->Fill(fullsimHZZEventTree.genl3pt, fullsimHZZEventTree.weight);
          histLep3Eta_e_fullsim->Fill(fullsimHZZEventTree.genl3eta, fullsimHZZEventTree.weight);
        }
      } else {
        if (fullsimHZZEventTree.genl3pt > 5) {
          histLep3Pt_m_fullsim->Fill(fullsimHZZEventTree.genl3pt, fullsimHZZEventTree.weight);
          histLep3Eta_m_fullsim->Fill(fullsimHZZEventTree.genl3eta, fullsimHZZEventTree.weight);
        }
      }

      if (abs(fullsimHZZEventTree.genl4id) == 11) {
        if (fullsimHZZEventTree.genl4pt > 7) {
          histLep4Pt_e_fullsim->Fill(fullsimHZZEventTree.genl4pt, fullsimHZZEventTree.weight);
          histLep4Eta_e_fullsim->Fill(fullsimHZZEventTree.genl4eta, fullsimHZZEventTree.weight);
        }
      } else {
        if (fullsimHZZEventTree.genl4pt > 5) {
          histLep4Pt_m_fullsim->Fill(fullsimHZZEventTree.genl4pt, fullsimHZZEventTree.weight);
          histLep4Eta_m_fullsim->Fill(fullsimHZZEventTree.genl4eta, fullsimHZZEventTree.weight);
        }
      }

      
    }//if pass event selection
   
  } //loop over full sim events

  //*****************************************************************************************
  // Loop over fastsim events
  //*****************************************************************************************
  for(UInt_t ientry=0; ientry < fastsimHZZEventTree.tree_->GetEntries(); ientry++) {       	
    fastsimHZZEventTree.tree_->GetEntry(ientry);
    if (ientry % 10000 == 0) cout << "Event " << ientry << endl;

    //how to get lepton four momenta

    cmsana::FourVectorM lepton1FourVector;
    lepton1FourVector.SetPt(fastsimHZZEventTree.l1pt);
    lepton1FourVector.SetEta(fastsimHZZEventTree.l1eta);
    lepton1FourVector.SetPhi(fastsimHZZEventTree.l1phi);
    lepton1FourVector.SetM( (abs(fastsimHZZEventTree.l1id) == 11) ? ELECTRONMASS : MUONMASS);

    //*****************************************************************************************
    // leading and subleading lepton
    //*****************************************************************************************
    double leadingLepPt = fastsimHZZEventTree.l1pt;
    double subleadingLepPt = fastsimHZZEventTree.l2pt;
    if (fastsimHZZEventTree.l2pt > fastsimHZZEventTree.l1pt) {
      leadingLepPt = fastsimHZZEventTree.l2pt;
      subleadingLepPt = fastsimHZZEventTree.l1pt;
    }
    if (fastsimHZZEventTree.l3pt > leadingLepPt) {
      subleadingLepPt = leadingLepPt;
      leadingLepPt = fastsimHZZEventTree.l3pt;      
    } else if (fastsimHZZEventTree.l3pt > subleadingLepPt) {
      subleadingLepPt = fastsimHZZEventTree.l3pt;
    }
    if (fastsimHZZEventTree.l4pt > leadingLepPt) {
      subleadingLepPt = leadingLepPt;
      leadingLepPt = fastsimHZZEventTree.l4pt;      
    } else if (fastsimHZZEventTree.l4pt > subleadingLepPt) {
      subleadingLepPt = fastsimHZZEventTree.l4pt;
    }


    //********************************
    //weight modifier
    //********************************
    double weightModifier = 1.0;
    if (abs(fastsimHZZEventTree.l1id) == 11 && 
        fastsimHZZEventTree.genl1pt < 8) weightModifier *= 1.15;
    if (abs(fastsimHZZEventTree.l1id) == 11 && 
        fastsimHZZEventTree.genl1pt >= 8 && fastsimHZZEventTree.genl1pt < 10) weightModifier *= 1.05;
    if (abs(fastsimHZZEventTree.l1id) == 11 && 
        fastsimHZZEventTree.genl1pt >= 10 &&  fastsimHZZEventTree.genl1pt < 16) weightModifier *= 1.03;
    if (abs(fastsimHZZEventTree.l1id) == 11 && 
        fastsimHZZEventTree.genl1pt >= 16 &&  fastsimHZZEventTree.genl1pt < 20) weightModifier *= 1.01;
    if (abs(fastsimHZZEventTree.l1id) == 11 && 
        fastsimHZZEventTree.genl1pt >= 20 &&  fastsimHZZEventTree.genl1pt < 40) weightModifier *= 0.975;
    if (abs(fastsimHZZEventTree.l1id) == 11 && 
        fastsimHZZEventTree.genl1pt >= 40 ) weightModifier *= 0.97;

    if (abs(fastsimHZZEventTree.l2id) == 11 && 
        fastsimHZZEventTree.genl2pt < 8) weightModifier *= 1.15;
    if (abs(fastsimHZZEventTree.l2id) == 11 && 
        fastsimHZZEventTree.genl2pt >= 8 && fastsimHZZEventTree.genl2pt < 10) weightModifier *= 1.05;
    if (abs(fastsimHZZEventTree.l1id) == 11 && 
        fastsimHZZEventTree.genl2pt >= 10 &&  fastsimHZZEventTree.genl2pt < 16) weightModifier *= 1.03;
    if (abs(fastsimHZZEventTree.l1id) == 11 && 
        fastsimHZZEventTree.genl2pt >= 16 &&  fastsimHZZEventTree.genl2pt < 20) weightModifier *= 1.01;
    if (abs(fastsimHZZEventTree.l2id) == 11 && 
        fastsimHZZEventTree.genl2pt >= 20 &&  fastsimHZZEventTree.genl2pt < 40) weightModifier *= 0.975;
    if (abs(fastsimHZZEventTree.l2id) == 11 && 
        fastsimHZZEventTree.genl2pt >= 40 ) weightModifier *= 0.97;

    if (abs(fastsimHZZEventTree.l3id) == 11 && 
        fastsimHZZEventTree.genl3pt < 8) weightModifier *= 1.15;
    if (abs(fastsimHZZEventTree.l3id) == 11 && 
        fastsimHZZEventTree.genl3pt >= 8 && fastsimHZZEventTree.genl3pt < 10) weightModifier *= 1.05;
    if (abs(fastsimHZZEventTree.l1id) == 11 && 
        fastsimHZZEventTree.genl3pt >= 10 &&  fastsimHZZEventTree.genl3pt < 16) weightModifier *= 1.03;
    if (abs(fastsimHZZEventTree.l3id) == 11 && 
        fastsimHZZEventTree.genl3pt >= 16 &&  fastsimHZZEventTree.genl3pt < 20) weightModifier *= 1.01;
    if (abs(fastsimHZZEventTree.l3id) == 11 && 
        fastsimHZZEventTree.genl3pt >= 20 &&  fastsimHZZEventTree.genl3pt < 40) weightModifier *= 0.975;
    if (abs(fastsimHZZEventTree.l3id) == 11 && 
        fastsimHZZEventTree.genl3pt >= 40 ) weightModifier *= 0.97;

    if (abs(fastsimHZZEventTree.l4id) == 11 && 
        fastsimHZZEventTree.genl4pt < 8) weightModifier *= 1.15;
    if (abs(fastsimHZZEventTree.l4id) == 11 && 
        fastsimHZZEventTree.genl4pt >= 8 && fastsimHZZEventTree.genl4pt < 10) weightModifier *= 1.05;
    if (abs(fastsimHZZEventTree.l1id) == 11 && 
        fastsimHZZEventTree.genl4pt >= 10 &&  fastsimHZZEventTree.genl4pt < 16) weightModifier *= 1.03;
    if (abs(fastsimHZZEventTree.l1id) == 11 && 
        fastsimHZZEventTree.genl4pt >= 16 &&  fastsimHZZEventTree.genl4pt < 20) weightModifier *= 1.01;
    if (abs(fastsimHZZEventTree.l4id) == 11 && 
        fastsimHZZEventTree.genl4pt >= 20 &&  fastsimHZZEventTree.genl4pt < 40) weightModifier *= 0.975;
    if (abs(fastsimHZZEventTree.l4id) == 11 && 
        fastsimHZZEventTree.genl4pt >= 40) weightModifier *= 0.97;

    if (abs(fastsimHZZEventTree.l1id) == 13 && fastsimHZZEventTree.genl1pt < 6) weightModifier *= 1.10;
    if (abs(fastsimHZZEventTree.l1id) == 13 && 
        fastsimHZZEventTree.genl1pt >= 6 && fastsimHZZEventTree.genl1pt < 7) weightModifier *= 1.05;
    if (abs(fastsimHZZEventTree.l1id) == 13 && 
        fastsimHZZEventTree.genl1pt >= 7 && fastsimHZZEventTree.genl1pt < 10) weightModifier *= 1.015;
    if (abs(fastsimHZZEventTree.l1id) == 13 && 
        fastsimHZZEventTree.genl1pt >= 10 &&  fastsimHZZEventTree.genl1pt < 20) weightModifier *= 1.015;
    if (abs(fastsimHZZEventTree.l1id) == 13 && 
        fastsimHZZEventTree.genl1pt >= 20 &&  fastsimHZZEventTree.genl1pt < 40) weightModifier *= 1.00;
    if (abs(fastsimHZZEventTree.l1id) == 13 &&  
        fastsimHZZEventTree.genl1pt >= 40 ) weightModifier *= 0.99;

    if (abs(fastsimHZZEventTree.l2id) == 13 && fastsimHZZEventTree.genl2pt < 6) weightModifier *= 1.10;
    if (abs(fastsimHZZEventTree.l2id) == 13 && 
        fastsimHZZEventTree.genl2pt >= 6 && fastsimHZZEventTree.genl2pt < 7) weightModifier *= 1.05;
    if (abs(fastsimHZZEventTree.l2id) == 13 && 
        fastsimHZZEventTree.genl2pt >= 7 && fastsimHZZEventTree.genl2pt < 10) weightModifier *= 1.015;
    if (abs(fastsimHZZEventTree.l1id) == 13 && 
        fastsimHZZEventTree.genl2pt >= 10 &&  fastsimHZZEventTree.genl2pt < 20) weightModifier *= 1.015;
    if (abs(fastsimHZZEventTree.l2id) == 13 && 
        fastsimHZZEventTree.genl2pt >= 20 &&  fastsimHZZEventTree.genl2pt < 40) weightModifier *= 1.00;
    if (abs(fastsimHZZEventTree.l2id) == 13 && 
        fastsimHZZEventTree.genl2pt >= 40 ) weightModifier *= 0.99;

    if (abs(fastsimHZZEventTree.l3id) == 13 && fastsimHZZEventTree.genl3pt < 6) weightModifier *= 1.10;
    if (abs(fastsimHZZEventTree.l3id) == 13 && 
        fastsimHZZEventTree.genl3pt >= 6 && fastsimHZZEventTree.genl3pt < 7) weightModifier *= 1.05;
    if (abs(fastsimHZZEventTree.l3id) == 13 && 
        fastsimHZZEventTree.genl3pt >= 7 && fastsimHZZEventTree.genl3pt < 10) weightModifier *= 1.015;
    if (abs(fastsimHZZEventTree.l1id) == 13 && 
        fastsimHZZEventTree.genl3pt >= 10 &&  fastsimHZZEventTree.genl3pt < 20) weightModifier *= 1.015;
    if (abs(fastsimHZZEventTree.l3id) == 13 && 
        fastsimHZZEventTree.genl3pt >= 20 &&  fastsimHZZEventTree.genl3pt < 40) weightModifier *= 1.00;
    if (abs(fastsimHZZEventTree.l3id) == 13 && 
        fastsimHZZEventTree.genl3pt >= 40 ) weightModifier *= 0.99;

    if (abs(fastsimHZZEventTree.l4id) == 13 && fastsimHZZEventTree.genl4pt < 6) weightModifier *= 1.10;
    if (abs(fastsimHZZEventTree.l4id) == 13 && 
        fastsimHZZEventTree.genl4pt >= 6 && fastsimHZZEventTree.genl4pt < 7) weightModifier *= 1.05;
    if (abs(fastsimHZZEventTree.l4id) == 13 && 
        fastsimHZZEventTree.genl4pt >= 7 && fastsimHZZEventTree.genl4pt < 10) weightModifier *= 1.015;
    if (abs(fastsimHZZEventTree.l1id) == 13 && 
        fastsimHZZEventTree.genl4pt >= 10 &&  fastsimHZZEventTree.genl4pt < 20) weightModifier *= 1.015;
    if (abs(fastsimHZZEventTree.l4id) == 13 && 
        fastsimHZZEventTree.genl4pt >= 20 &&  fastsimHZZEventTree.genl4pt < 40) weightModifier *= 1.00;
    if (abs(fastsimHZZEventTree.l4id) == 13 && 
        fastsimHZZEventTree.genl4pt >= 40) weightModifier *= 0.99;





    //fill events passing selection
    if (fastsimHZZEventTree.l1pt > (abs(fastsimHZZEventTree.l1id) == 11 ? 7 : 5)
        && fastsimHZZEventTree.l2pt > (abs(fastsimHZZEventTree.l2id) == 11 ? 7 : 5)
        && fastsimHZZEventTree.l3pt > (abs(fastsimHZZEventTree.l3id) == 11 ? 7 : 5)
        && fastsimHZZEventTree.l4pt > (abs(fastsimHZZEventTree.l4id) == 11 ? 7 : 5)
        && abs(fastsimHZZEventTree.l1eta) < (abs(fastsimHZZEventTree.l1id) == 11 ? 2.5 : 2.4)
        && abs(fastsimHZZEventTree.l2eta) < (abs(fastsimHZZEventTree.l2id) == 11 ? 2.5 : 2.4)
        && abs(fastsimHZZEventTree.l3eta) < (abs(fastsimHZZEventTree.l3id) == 11 ? 2.5 : 2.4)
        && abs(fastsimHZZEventTree.l4eta) < (abs(fastsimHZZEventTree.l4id) == 11 ? 2.5 : 2.4)
        && leadingLepPt > 20
        && subleadingLepPt > 10
      ) {

      fastsimYield += weightModifier*fastsimHZZEventTree.weight;

      //******************************************
      //ZMass
      //******************************************
      if (abs(fastsimHZZEventTree.l1id) == 11) {
        histZ1Mass_ee_fastsim->Fill(fastsimHZZEventTree.genz1mass, weightModifier*fastsimHZZEventTree.weight);
      } else {
        histZ1Mass_mm_fastsim->Fill(fastsimHZZEventTree.genz1mass, weightModifier*fastsimHZZEventTree.weight);
      }
      if (abs(fastsimHZZEventTree.l3id) == 11) {
        histZ2Mass_ee_fastsim->Fill(fastsimHZZEventTree.genz2mass, weightModifier*fastsimHZZEventTree.weight);
      } else {
        histZ2Mass_mm_fastsim->Fill(fastsimHZZEventTree.genz2mass, weightModifier*fastsimHZZEventTree.weight);
      }
      if (fastsimHZZEventTree.channel == 0) histZZMass_4e_fastsim->Fill(fastsimHZZEventTree.genzzmass, weightModifier*fastsimHZZEventTree.weight);
      if (fastsimHZZEventTree.channel == 1) histZZMass_4m_fastsim->Fill(fastsimHZZEventTree.genzzmass, weightModifier*fastsimHZZEventTree.weight);
      if (fastsimHZZEventTree.channel == 2) histZZMass_2e2m_fastsim->Fill(fastsimHZZEventTree.genzzmass, weightModifier*fastsimHZZEventTree.weight);
      if (fastsimHZZEventTree.channel == 3) histZZMass_2m2e_fastsim->Fill(fastsimHZZEventTree.genzzmass, weightModifier*fastsimHZZEventTree.weight);

      //******************************************
      //Lepton Pt
      //******************************************
      if (abs(fastsimHZZEventTree.genl1id) == 11) {
        if (fastsimHZZEventTree.genl1pt > 7) {
          histLep1Pt_e_fastsim->Fill(fastsimHZZEventTree.genl1pt, weightModifier*fastsimHZZEventTree.weight);
          histLep1Eta_e_fastsim->Fill(fastsimHZZEventTree.genl1eta, weightModifier*fastsimHZZEventTree.weight);
        }
      } else {
        if (fastsimHZZEventTree.genl1pt > 5) {
          histLep1Pt_m_fastsim->Fill(fastsimHZZEventTree.genl1pt, weightModifier*fastsimHZZEventTree.weight);
          histLep1Eta_m_fastsim->Fill(fastsimHZZEventTree.genl1eta, weightModifier*fastsimHZZEventTree.weight);
        }
      }
      if (abs(fastsimHZZEventTree.genl2id) == 11) {
        if (fastsimHZZEventTree.genl2pt > 7) {
          histLep2Pt_e_fastsim->Fill(fastsimHZZEventTree.genl2pt, weightModifier*fastsimHZZEventTree.weight);
          histLep2Eta_e_fastsim->Fill(fastsimHZZEventTree.genl2eta, weightModifier*fastsimHZZEventTree.weight);
        }
      } else {
        if (fastsimHZZEventTree.genl2pt > 5) {
          histLep2Pt_m_fastsim->Fill(fastsimHZZEventTree.genl2pt, weightModifier*fastsimHZZEventTree.weight);
          histLep2Eta_m_fastsim->Fill(fastsimHZZEventTree.genl2eta, weightModifier*fastsimHZZEventTree.weight);
        }
      }

      if (abs(fastsimHZZEventTree.genl3id) == 11) {
        if (fastsimHZZEventTree.genl3pt > 7) {
          histLep3Pt_e_fastsim->Fill(fastsimHZZEventTree.genl3pt, weightModifier*fastsimHZZEventTree.weight);
          histLep3Eta_e_fastsim->Fill(fastsimHZZEventTree.genl3eta, weightModifier*fastsimHZZEventTree.weight);
        }
      } else {
        if (fastsimHZZEventTree.genl3pt > 5) {
          histLep3Pt_m_fastsim->Fill(fastsimHZZEventTree.genl3pt, weightModifier*fastsimHZZEventTree.weight);
          histLep3Eta_m_fastsim->Fill(fastsimHZZEventTree.genl3eta, weightModifier*fastsimHZZEventTree.weight);
        }
      }

      if (abs(fastsimHZZEventTree.genl4id) == 11) {
        if (fastsimHZZEventTree.genl4pt > 7) {
          histLep4Pt_e_fastsim->Fill(fastsimHZZEventTree.genl4pt, weightModifier*fastsimHZZEventTree.weight);
          histLep4Eta_e_fastsim->Fill(fastsimHZZEventTree.genl4eta, weightModifier*fastsimHZZEventTree.weight);
        }
      } else {
        if (fastsimHZZEventTree.genl4pt > 5) {
          histLep4Pt_m_fastsim->Fill(fastsimHZZEventTree.genl4pt, weightModifier*fastsimHZZEventTree.weight);
          histLep4Eta_m_fastsim->Fill(fastsimHZZEventTree.genl4eta, weightModifier*fastsimHZZEventTree.weight);
        }
      }


    } //if pass event selection

  } //loop over fast sim events



  //*****************************************************************************************
  // Plot
  //*****************************************************************************************
  TCanvas *cv = 0;
  TLegend *legend = 0;


  //*****************************************************************************************
  // Lepton Pt
  //*****************************************************************************************
  DrawComparison(histLep1Pt_e_fullsim, histLep1Pt_e_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_Lep1Pt_e" + label + ".gif");
  DrawComparison(histLep1Pt_m_fullsim, histLep1Pt_m_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_Lep1Pt_m" + label + ".gif");
  DrawComparison(histLep2Pt_e_fullsim, histLep2Pt_e_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_Lep2Pt_e" + label + ".gif");
  DrawComparison(histLep2Pt_m_fullsim, histLep2Pt_m_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_Lep2Pt_m" + label + ".gif");
  DrawComparison(histLep3Pt_e_fullsim, histLep3Pt_e_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_Lep3Pt_e" + label + ".gif");
  DrawComparison(histLep3Pt_m_fullsim, histLep3Pt_m_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_Lep3Pt_m" + label + ".gif");
  DrawComparison(histLep4Pt_e_fullsim, histLep4Pt_e_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_Lep4Pt_e" + label + ".gif");
  DrawComparison(histLep4Pt_m_fullsim, histLep4Pt_m_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_Lep4Pt_m" + label + ".gif");

  //*****************************************************************************************
  // Lepton Eta
  //*****************************************************************************************
  DrawComparison(histLep1Eta_e_fullsim, histLep1Eta_e_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_Lep1Eta_e" + label + ".gif");
  DrawComparison(histLep1Eta_m_fullsim, histLep1Eta_m_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_Lep1Eta_m" + label + ".gif");
  DrawComparison(histLep2Eta_e_fullsim, histLep2Eta_e_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_Lep2Eta_e" + label + ".gif");
  DrawComparison(histLep2Eta_m_fullsim, histLep2Eta_m_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_Lep2Eta_m" + label + ".gif");
  DrawComparison(histLep3Eta_e_fullsim, histLep3Eta_e_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_Lep3Eta_e" + label + ".gif");
  DrawComparison(histLep3Eta_m_fullsim, histLep3Eta_m_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_Lep3Eta_m" + label + ".gif");
  DrawComparison(histLep4Eta_e_fullsim, histLep4Eta_e_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_Lep4Eta_e" + label + ".gif");
  DrawComparison(histLep4Eta_m_fullsim, histLep4Eta_m_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_Lep4Eta_m" + label + ".gif");

  //*****************************************************************************************
  // Z1 Mass
  //*****************************************************************************************
  DrawComparison(histZ1Mass_ee_fullsim, histZ1Mass_ee_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_Z1Mass_ee" + label + ".gif");
  DrawComparison(histZ1Mass_mm_fullsim, histZ1Mass_mm_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_Z1Mass_mm" + label + ".gif");
  DrawComparison(histZ2Mass_ee_fullsim, histZ2Mass_ee_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_Z2Mass_ee" + label + ".gif");
  DrawComparison(histZ2Mass_mm_fullsim, histZ2Mass_mm_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_Z2Mass_mm" + label + ".gif");

  DrawComparison(histZZMass_4e_fullsim, histZZMass_4e_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_ZZMass_4e" + label + ".gif");
  DrawComparison(histZZMass_4m_fullsim, histZZMass_4m_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_ZZMass_4m" + label + ".gif");
  DrawComparison(histZZMass_2e2m_fullsim, histZZMass_2e2m_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_ZZMass_2e2m" + label + ".gif");
  DrawComparison(histZZMass_2m2e_fullsim, histZZMass_2m2e_fastsim, "FullSim", "FastSim", 
                 "FastSimValidation_ZZMass_2m2e" + label + ".gif");

  //*****************************************************************************************
  // Compare Yield
  //*****************************************************************************************
  cout << "FullSim Yield : " << fullsimYield << "\n";
  cout << "FastSim Yield : " << fastsimYield << "\n";
  cout << "FastSim Yield / FullSim Yield = " << fastsimYield/fullsimYield << " \n";

}

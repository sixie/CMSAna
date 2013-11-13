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
#include "CMSAna/HZZ4L/interface/TauHelperFunctions2.h"
 
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

  for (UInt_t b=0; b < UInt_t(hist->GetXaxis()->GetNbins()+2); ++b) {
    if (reference->GetBinContent(b) > 0) {
      relative->SetBinContent(b,100 * ( hist->GetBinContent(b) - reference->GetBinContent(b) )/reference->GetBinContent(b));
    } else {
      relative->SetBinContent(b,0);
    }
  }

  return relative;
}

void DrawHist( TH1F* hist1, string plotname, bool doNormalize = false) {

  TCanvas *cv = new TCanvas("cv","cv", 800,600);

  TPad *pad1 = new TPad("pad1","pad1", 0,0.2,1,1);
  pad1->SetBottomMargin(0.125);
  pad1->Draw();
  pad1->cd();

  if (doNormalize) {
    NormalizeHist(hist1);
  }

  hist1->SetLineColor(kBlack);
  hist1->Draw("E1");

  cout << plotname << " :  underflows : " << hist1->GetBinContent(0)  << "\n";
  cout << plotname << " :  overflows : " << hist1->GetBinContent(hist1->GetXaxis()->GetNbins()+1) << "\n";
  cout << "Integral : " << hist1->Integral() << "\n";

  cv->SaveAs(plotname.c_str());

}

//*************************************************************************************************
//Convert Vectors to Angles
//*************************************************************************************************
EventParameters ConvertVectorsToAnglesRoberto(FourVector L11, FourVector L12, FourVector L21, FourVector L22)
{

//    FourVector L11 = Leptons.Lepton11;
//    FourVector L12 = Leptons.Lepton12;
//    FourVector L21 = Leptons.Lepton21;
//    FourVector L22 = Leptons.Lepton22;

   L11[0] = L11.GetP();
   L12[0] = L12.GetP();
   L21[0] = L21.GetP();
   L22[0] = L22.GetP();

   double HMass = (L11 + L12 + L21 + L22).GetMass();
   double ZMass = (L11 + L12).GetMass();
   double Z2Mass = (L21 + L22).GetMass();

   double Gamma1 = HMass / (2 * ZMass) * (1 + (ZMass * ZMass - Z2Mass * Z2Mass) / (HMass * HMass));
   double Beta1 = GammaToBeta(Gamma1);
   double Gamma2 = HMass / (2 * Z2Mass) * (1 - (ZMass * ZMass - Z2Mass * Z2Mass) / (HMass * HMass));
   double Beta2 = GammaToBeta(Gamma2);

   FourVector HiggsLab = L11 + L12 + L21 + L22;
   double HiggsBoostGamma = HiggsLab[0] / HMass;
   double HiggsBoostBeta = GammaToBeta(HiggsBoostGamma);
   double PhiH = HiggsLab.GetPhi();

   L11 = L11.GammaBoost(HiggsLab, -HiggsBoostGamma);
   L12 = L12.GammaBoost(HiggsLab, -HiggsBoostGamma);
   L21 = L21.GammaBoost(HiggsLab, -HiggsBoostGamma);
   L22 = L22.GammaBoost(HiggsLab, -HiggsBoostGamma);

   FourVector Z1 = L11 + L12;
   FourVector Z2 = L21 + L22;

   double Z1BoostBeta = GammaToBeta(Z1[0] / ZMass);
   double Z2BoostBeta = GammaToBeta(Z2[0] / Z2Mass);

   double Theta1 = GetAngle(Z1, L11.GammaBoost(Z1, -Z1[0] / ZMass));
   double Theta2 = GetAngle(Z2, L21.GammaBoost(Z2, -Z2[0] / Z2Mass));

   FourVector L11Perpendicular = L11 - Z1 * (L11.SpatialDot(Z1)) / Z1.GetP2();
   FourVector L21Perpendicular = L21 - Z1 * (L21.SpatialDot(Z1)) / Z1.GetP2();
   FourVector HPerpendicular = HiggsLab - Z1 * (HiggsLab.SpatialDot(Z1)) / Z1.GetP2();

   double Phi = GetAngle(-L21Perpendicular, L11Perpendicular) - PI;
   if(Z2.SpatialDot(L21Perpendicular.SpatialCross(L11Perpendicular)) < 0)
      Phi = 2 * PI - Phi;
   while(Phi < 0)         Phi = Phi + 2 * PI;
   while(Phi >= 2 * PI)   Phi = Phi - 2 * PI;

   FourVector P1CM(HMass / 2, 0, 0, HMass / 2);
   double CosTheta0 = -(P1CM * (L11 + L12) / (0.5 * ZMass * HMass * Gamma1) - 1) / Beta1;
   double Theta0 = acos(CosTheta0);

   double CosPhi1 = (P1CM * (L11 - L12)) / (0.5 * ZMass * HMass)
      - Gamma1 * (Beta1 - cos(Theta0)) * cos(Theta1);
   CosPhi1 = CosPhi1 / (fabs(sin(Theta0) * sin(Theta1)));

   double SinPhi1 = (P1CM.SpatialCross(L11 + L12).SpatialDot(L11))
      / (0.25 * ZMass * ZMass * HMass * Beta1 * Gamma1 * fabs(sin(Theta0)) * sin(Theta1));

   double Phi1 = PI - acos(CosPhi1);
   if(SinPhi1 < 0)
      Phi1 = -Phi1;

   if(Phi1 < 0)         Phi1 = Phi1 + 2 * PI;
   if(Phi1 >= 2 * PI)   Phi1 = Phi1 - 2 * PI;
   
   FourVector Lepton11;
   Lepton11[0] = Gamma1 * (1 + Beta1 * cos(Theta1));
   Lepton11[1] = Gamma1 * (Beta1 + cos(Theta1)) * sin(Theta0) + cos(Theta0) * cos(PI-Phi1) * sin(Theta1);
   Lepton11[2] = sin(Theta1) * sin(Phi1);
   Lepton11[3] = Gamma1 * cos(Theta0) * (Beta1 + cos(Theta1)) - cos(PI-Phi1) * sin(Theta0) * sin(Theta1);
   Lepton11 = Lepton11 * ZMass / 2;
   
   EventParameters Result;

   Result.Phi0 = Phi1;
   Result.Theta0 = Theta0;
   Result.Phi = Phi;
   Result.Theta1 = Theta1;
   Result.Theta2 = Theta2;
   Result.HMass = HMass;
   Result.ZMass = ZMass;
   Result.Z2Mass = Z2Mass;
   Result.PhiH = PhiH;
   Result.PhiOffset = -(Lepton11.GetPhi() - L11.GetPhi());
   Result.PTH = HiggsLab.GetPT();
   Result.YH = HiggsLab.GetY();

   return Result;
}



//*************************************************************************************************
//Main Function
//*************************************************************************************************
void GenerateTemplatePdf(string inputFilename, 
                         const string Label = "ZJets") {

  string label = Label;
  if (Label != "") label = "_" + Label;
  bool doNormalize = true;

  //*****************************************************************************************
  // Load Input
  //*****************************************************************************************
  HZZEventTree fastsimHZZEventTree;
  fastsimHZZEventTree.LoadTree(inputFilename.c_str());
  fastsimHZZEventTree.InitTree();
  cout << "Events in the ntuple: " << fastsimHZZEventTree.tree_->GetEntries() << endl;

  //*************************************************************************************************
  //Histograms
  //*************************************************************************************************
  TH1::AddDirectory(kFALSE);
  TH1F *histZ1Mass_ee_fastsim = new TH1F( "histZ1Mass_ee_fastsim", ";Z1 Mass [GeV/c^{2}]; Number of Events", 100, 0,200);
  TH1F *histZ1Mass_mm_fastsim = new TH1F( "histZ1Mass_mm_fastsim", ";Z1 Mass [GeV/c^{2}]; Number of Events", 100, 0,200);
  TH1F *histZ2Mass_ee_fastsim = new TH1F( "histZ2Mass_ee_fastsim", ";Z2 Mass [GeV/c^{2}]; Number of Events", 100, 0,200);
  TH1F *histZ2Mass_mm_fastsim = new TH1F( "histZ2Mass_mm_fastsim", ";Z2 Mass [GeV/c^{2}]; Number of Events", 100, 0,200);

  int zzmassbins = 100;
  double zzmassmin = 100;
  double zzmassmax = 150;
  TH1F *histZZMass_4e_fastsim = new TH1F( "histZZMass_4e_fastsim", ";ZZ Mass [GeV/c^{2}]; Number of Events", zzmassbins, zzmassmin,zzmassmax);
  TH1F *histZZMass_4m_fastsim = new TH1F( "histZZMass_4m_fastsim", ";ZZ Mass [GeV/c^{2}]; Number of Events", zzmassbins, zzmassmin,zzmassmax);
  TH1F *histZZMass_2e2m_fastsim = new TH1F( "histZZMass_2e2m_fastsim", ";ZZ Mass [GeV/c^{2}]; Number of Events", zzmassbins, zzmassmin,zzmassmax);
  TH1F *histZZMass_2m2e_fastsim = new TH1F( "histZZMass_2m2e_fastsim", ";ZZ Mass [GeV/c^{2}]; Number of Events", zzmassbins, zzmassmin,zzmassmax);


  TH1F *histLep1Pt_e_fastsim = new TH1F( "histLep1Pt_e_fastsim", ";lep1 pT [GeV/c]; Number of Events", 100, 0,200);
  TH1F *histLep1Pt_m_fastsim = new TH1F( "histLep1Pt_m_fastsim", ";lep1 pT [GeV/c]; Number of Events", 100, 0,200);
  TH1F *histLep2Pt_e_fastsim = new TH1F( "histLep2Pt_e_fastsim", ";lep2 pT [GeV/c]; Number of Events", 100, 0,100);
  TH1F *histLep2Pt_m_fastsim = new TH1F( "histLep2Pt_m_fastsim", ";lep2 pT [GeV/c]; Number of Events", 100, 0,100);
  TH1F *histLep3Pt_e_fastsim = new TH1F( "histLep3Pt_e_fastsim", ";lep3 pT [GeV/c]; Number of Events", 100, 0,100);
  TH1F *histLep3Pt_m_fastsim = new TH1F( "histLep3Pt_m_fastsim", ";lep3 pT [GeV/c]; Number of Events", 100, 0,100);
  TH1F *histLep4Pt_e_fastsim = new TH1F( "histLep4Pt_e_fastsim", ";lep4 pT [GeV/c]; Number of Events", 100, 0,100);
  TH1F *histLep4Pt_m_fastsim = new TH1F( "histLep4Pt_m_fastsim", ";lep4 pT [GeV/c]; Number of Events", 100, 0,100);

  TH1F *histLep1Eta_e_fastsim = new TH1F( "histLep1Eta_e_fastsim", ";lep1 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep1Eta_m_fastsim = new TH1F( "histLep1Eta_m_fastsim", ";lep1 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep2Eta_e_fastsim = new TH1F( "histLep2Eta_e_fastsim", ";lep2 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep2Eta_m_fastsim = new TH1F( "histLep2Eta_m_fastsim", ";lep2 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep3Eta_e_fastsim = new TH1F( "histLep3Eta_e_fastsim", ";lep3 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep3Eta_m_fastsim = new TH1F( "histLep3Eta_m_fastsim", ";lep3 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep4Eta_e_fastsim = new TH1F( "histLep4Eta_e_fastsim", ";lep4 #eta; Number of Events", 100, -2.5,2.5);
  TH1F *histLep4Eta_m_fastsim = new TH1F( "histLep4Eta_m_fastsim", ";lep4 #eta; Number of Events", 100, -2.5,2.5);  


  TH1F *histPhi0_4e_fastsim = new TH1F( "histPhi0_4e_fastsim", ";#Phi_{0}; Number of Events", 100, 0,6.4);
  TH1F *histTheta0_4e_fastsim = new TH1F( "histTheta0_4e_fastsim", ";#Theta_{0}; Number of Events", 100, 0,3.2);
  TH1F *histPhi_4e_fastsim = new TH1F( "histPhi_4e_fastsim", ";#Phi; Number of Events", 100, 0,6.4);
  TH1F *histTheta1_4e_fastsim = new TH1F( "histTheta1_4e_fastsim", ";#Theta_{1}; Number of Events", 100, 0,3.2);
  TH1F *histTheta2_4e_fastsim = new TH1F( "histTheta2_4e_fastsim", ";#Theta_{2}; Number of Events", 100, 0,3.2);
  TH1F *histPhi0_4m_fastsim = new TH1F( "histPhi0_4m_fastsim", ";#Phi_{0}; Number of Events", 100, 0,6.4);
  TH1F *histTheta0_4m_fastsim = new TH1F( "histTheta0_4m_fastsim", ";#Theta_{0}; Number of Events", 100, 0,3.2);
  TH1F *histPhi_4m_fastsim = new TH1F( "histPhi_4m_fastsim", ";#Phi; Number of Events", 100, 0,6.4);
  TH1F *histTheta1_4m_fastsim = new TH1F( "histTheta1_4m_fastsim", ";#Theta_{1}; Number of Events", 100, 0,3.2);
  TH1F *histTheta2_4m_fastsim = new TH1F( "histTheta2_4m_fastsim", ";#Theta_{2}; Number of Events", 100, 0,3.2);
  TH1F *histPhi0_2e2m_fastsim = new TH1F( "histPhi0_2e2m_fastsim", ";#Phi_{0}; Number of Events", 100, 0,6.4);
  TH1F *histTheta0_2e2m_fastsim = new TH1F( "histTheta0_2e2m_fastsim", ";#Theta_{0}; Number of Events", 100, 0,3.2);
  TH1F *histPhi_2e2m_fastsim = new TH1F( "histPhi_2e2m_fastsim", ";#Phi; Number of Events", 100, 0,6.4);
  TH1F *histTheta1_2e2m_fastsim = new TH1F( "histTheta1_2e2m_fastsim", ";#Theta_{1}; Number of Events", 100, 0,3.2);
  TH1F *histTheta2_2e2m_fastsim = new TH1F( "histTheta2_2e2m_fastsim", ";#Theta_{2}; Number of Events", 100, 0,3.2);
  TH1F *histPhi0_2m2e_fastsim = new TH1F( "histPhi0_2m2e_fastsim", ";#Phi_{0}; Number of Events", 100, 0,6.4);
  TH1F *histTheta0_2m2e_fastsim = new TH1F( "histTheta0_2m2e_fastsim", ";#Theta_{0}; Number of Events", 100, 0,3.2);
  TH1F *histPhi_2m2e_fastsim = new TH1F( "histPhi_2m2e_fastsim", ";#Phi; Number of Events", 100, 0,6.4);
  TH1F *histTheta1_2m2e_fastsim = new TH1F( "histTheta1_2m2e_fastsim", ";#Theta_{1}; Number of Events", 100, 0,3.2);
  TH1F *histTheta2_2m2e_fastsim = new TH1F( "histTheta2_2m2e_fastsim", ";#Theta_{2}; Number of Events", 100, 0,3.2);


  histZ1Mass_ee_fastsim->Sumw2();
  histZ1Mass_mm_fastsim->Sumw2();
  histZ2Mass_ee_fastsim->Sumw2();
  histZ2Mass_mm_fastsim->Sumw2();
  histZZMass_4e_fastsim->Sumw2();
  histZZMass_4m_fastsim->Sumw2();
  histZZMass_2e2m_fastsim->Sumw2();
  histZZMass_2m2e_fastsim->Sumw2();
  histLep1Pt_e_fastsim->Sumw2();
  histLep1Pt_m_fastsim->Sumw2();
  histLep2Pt_e_fastsim->Sumw2();
  histLep2Pt_m_fastsim->Sumw2();
  histLep3Pt_e_fastsim->Sumw2();
  histLep3Pt_m_fastsim->Sumw2();
  histLep4Pt_e_fastsim->Sumw2();
  histLep4Pt_m_fastsim->Sumw2();
  histLep1Eta_e_fastsim->Sumw2();
  histLep1Eta_m_fastsim->Sumw2();
  histLep2Eta_e_fastsim->Sumw2();
  histLep2Eta_m_fastsim->Sumw2();
  histLep3Eta_e_fastsim->Sumw2();
  histLep3Eta_m_fastsim->Sumw2();
  histLep4Eta_e_fastsim->Sumw2();
  histLep4Eta_m_fastsim->Sumw2();
  histPhi0_4e_fastsim->Sumw2();
  histTheta0_4e_fastsim->Sumw2();
  histPhi_4e_fastsim->Sumw2();
  histTheta1_4e_fastsim->Sumw2();
  histTheta2_4e_fastsim->Sumw2();
  histPhi0_4m_fastsim->Sumw2();
  histTheta0_4m_fastsim->Sumw2();
  histPhi_4m_fastsim->Sumw2();
  histTheta1_4m_fastsim->Sumw2();
  histTheta2_4m_fastsim->Sumw2();
  histPhi0_2e2m_fastsim->Sumw2();
  histTheta0_2e2m_fastsim->Sumw2();
  histPhi_2e2m_fastsim->Sumw2();
  histTheta1_2e2m_fastsim->Sumw2();
  histTheta2_2e2m_fastsim->Sumw2();
  histPhi0_2m2e_fastsim->Sumw2();
  histTheta0_2m2e_fastsim->Sumw2();
  histPhi_2m2e_fastsim->Sumw2();
  histTheta1_2m2e_fastsim->Sumw2();
  histTheta2_2m2e_fastsim->Sumw2();




  //*************************************************************************************************
  //Yields
  //*************************************************************************************************
  double fastsimYield = 0;


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

      FourVector lepton1FourVector;
      lepton1FourVector.SetPtEtaPhiMass( fastsimHZZEventTree.l1pt, 
                                         fastsimHZZEventTree.l1eta,
                                         fastsimHZZEventTree.l1phi,
                                         (abs(fastsimHZZEventTree.l1id) == 11) ? ELECTRONMASS : MUONMASS);
      FourVector lepton2FourVector;
      lepton2FourVector.SetPtEtaPhiMass( fastsimHZZEventTree.l2pt, 
                                         fastsimHZZEventTree.l2eta,
                                         fastsimHZZEventTree.l2phi,
                                         (abs(fastsimHZZEventTree.l2id) == 11) ? ELECTRONMASS : MUONMASS);
      FourVector lepton3FourVector;
      lepton3FourVector.SetPtEtaPhiMass( fastsimHZZEventTree.l3pt, 
                                         fastsimHZZEventTree.l3eta,
                                         fastsimHZZEventTree.l3phi,
                                         (abs(fastsimHZZEventTree.l3id) == 11) ? ELECTRONMASS : MUONMASS);
      FourVector lepton4FourVector;
      lepton4FourVector.SetPtEtaPhiMass( fastsimHZZEventTree.l4pt, 
                                         fastsimHZZEventTree.l4eta,
                                         fastsimHZZEventTree.l4phi,
                                         (abs(fastsimHZZEventTree.l4id) == 11) ? ELECTRONMASS : MUONMASS);      

      EventParameters eventparameters = ConvertVectorsToAnglesRoberto( lepton1FourVector, lepton2FourVector,
                                                                       lepton3FourVector, lepton4FourVector);


      fastsimYield += fastsimHZZEventTree.weight;

      //******************************************
      //ZMass
      //******************************************
      if (abs(fastsimHZZEventTree.l1id) == 11) {
        histZ1Mass_ee_fastsim->Fill(fastsimHZZEventTree.z1mass, fastsimHZZEventTree.weight);
      } else {
        histZ1Mass_mm_fastsim->Fill(fastsimHZZEventTree.z1mass, fastsimHZZEventTree.weight);
      }
      if (abs(fastsimHZZEventTree.l3id) == 11) {
        histZ2Mass_ee_fastsim->Fill(fastsimHZZEventTree.z2mass, fastsimHZZEventTree.weight);
      } else {
        histZ2Mass_mm_fastsim->Fill(fastsimHZZEventTree.z2mass, fastsimHZZEventTree.weight);
      }
      if (fastsimHZZEventTree.channel == 0) histZZMass_4e_fastsim->Fill(fastsimHZZEventTree.zzmass, fastsimHZZEventTree.weight);
      if (fastsimHZZEventTree.channel == 1) histZZMass_4m_fastsim->Fill(fastsimHZZEventTree.zzmass, fastsimHZZEventTree.weight);
      if (fastsimHZZEventTree.channel == 2) histZZMass_2e2m_fastsim->Fill(fastsimHZZEventTree.zzmass, fastsimHZZEventTree.weight);
      if (fastsimHZZEventTree.channel == 3) histZZMass_2m2e_fastsim->Fill(fastsimHZZEventTree.zzmass, fastsimHZZEventTree.weight);

      //******************************************
      //Lepton Pt
      //******************************************
      if (abs(fastsimHZZEventTree.l1id) == 11) {
        if (fastsimHZZEventTree.l1pt > 7) {
          histLep1Pt_e_fastsim->Fill(fastsimHZZEventTree.l1pt, fastsimHZZEventTree.weight);
          histLep1Eta_e_fastsim->Fill(fastsimHZZEventTree.l1eta, fastsimHZZEventTree.weight);
        }
      } else {
        if (fastsimHZZEventTree.l1pt > 5) {
          histLep1Pt_m_fastsim->Fill(fastsimHZZEventTree.l1pt, fastsimHZZEventTree.weight);
          histLep1Eta_m_fastsim->Fill(fastsimHZZEventTree.l1eta, fastsimHZZEventTree.weight);
        }
      }
      if (abs(fastsimHZZEventTree.l2id) == 11) {
        if (fastsimHZZEventTree.l2pt > 7) {
          histLep2Pt_e_fastsim->Fill(fastsimHZZEventTree.l2pt, fastsimHZZEventTree.weight);
          histLep2Eta_e_fastsim->Fill(fastsimHZZEventTree.l2eta, fastsimHZZEventTree.weight);
        }
      } else {
        if (fastsimHZZEventTree.l2pt > 5) {
          histLep2Pt_m_fastsim->Fill(fastsimHZZEventTree.l2pt, fastsimHZZEventTree.weight);
          histLep2Eta_m_fastsim->Fill(fastsimHZZEventTree.l2eta, fastsimHZZEventTree.weight);
        }
      }

      if (abs(fastsimHZZEventTree.l3id) == 11) {
        if (fastsimHZZEventTree.l3pt > 7) {
          histLep3Pt_e_fastsim->Fill(fastsimHZZEventTree.l3pt, fastsimHZZEventTree.weight);
          histLep3Eta_e_fastsim->Fill(fastsimHZZEventTree.l3eta, fastsimHZZEventTree.weight);
        }
      } else {
        if (fastsimHZZEventTree.l3pt > 5) {
          histLep3Pt_m_fastsim->Fill(fastsimHZZEventTree.l3pt, fastsimHZZEventTree.weight);
          histLep3Eta_m_fastsim->Fill(fastsimHZZEventTree.l3eta, fastsimHZZEventTree.weight);
        }
      }

      if (abs(fastsimHZZEventTree.l4id) == 11) {
        if (fastsimHZZEventTree.l4pt > 7) {
          histLep4Pt_e_fastsim->Fill(fastsimHZZEventTree.l4pt, fastsimHZZEventTree.weight);
          histLep4Eta_e_fastsim->Fill(fastsimHZZEventTree.l4eta, fastsimHZZEventTree.weight);
        }
      } else {
        if (fastsimHZZEventTree.l4pt > 5) {
          histLep4Pt_m_fastsim->Fill(fastsimHZZEventTree.l4pt, fastsimHZZEventTree.weight);
          histLep4Eta_m_fastsim->Fill(fastsimHZZEventTree.l4eta, fastsimHZZEventTree.weight);
        }
      }

      //******************************************
      //Angles
      //******************************************
      if (fastsimHZZEventTree.channel == 0) {
        histPhi0_4e_fastsim->Fill(eventparameters.Phi0, fastsimHZZEventTree.weight);
        histTheta0_4e_fastsim->Fill(eventparameters.Theta0, fastsimHZZEventTree.weight);
        histPhi_4e_fastsim->Fill(eventparameters.Phi, fastsimHZZEventTree.weight);
        histTheta1_4e_fastsim->Fill(eventparameters.Theta1, fastsimHZZEventTree.weight);
        histTheta2_4e_fastsim->Fill(eventparameters.Theta2, fastsimHZZEventTree.weight);
      } else if (fastsimHZZEventTree.channel == 1) {
        histPhi0_4m_fastsim->Fill(eventparameters.Phi0, fastsimHZZEventTree.weight);
        histTheta0_4m_fastsim->Fill(eventparameters.Theta0, fastsimHZZEventTree.weight);
        histPhi_4m_fastsim->Fill(eventparameters.Phi, fastsimHZZEventTree.weight);
        histTheta1_4m_fastsim->Fill(eventparameters.Theta1, fastsimHZZEventTree.weight);
        histTheta2_4m_fastsim->Fill(eventparameters.Theta2, fastsimHZZEventTree.weight);
      } else if (fastsimHZZEventTree.channel == 2) {
        histPhi0_2e2m_fastsim->Fill(eventparameters.Phi0, fastsimHZZEventTree.weight);
        histTheta0_2e2m_fastsim->Fill(eventparameters.Theta0, fastsimHZZEventTree.weight);
        histPhi_2e2m_fastsim->Fill(eventparameters.Phi, fastsimHZZEventTree.weight);
        histTheta1_2e2m_fastsim->Fill(eventparameters.Theta1, fastsimHZZEventTree.weight);
        histTheta2_2e2m_fastsim->Fill(eventparameters.Theta2, fastsimHZZEventTree.weight);
      } else if (fastsimHZZEventTree.channel == 3) {
        histPhi0_2m2e_fastsim->Fill(eventparameters.Phi0, fastsimHZZEventTree.weight);
        histTheta0_2m2e_fastsim->Fill(eventparameters.Theta0, fastsimHZZEventTree.weight);
        histPhi_2m2e_fastsim->Fill(eventparameters.Phi, fastsimHZZEventTree.weight);
        histTheta1_2m2e_fastsim->Fill(eventparameters.Theta1, fastsimHZZEventTree.weight);
        histTheta2_2m2e_fastsim->Fill(eventparameters.Theta2, fastsimHZZEventTree.weight);
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
  DrawHist(histLep1Pt_e_fastsim, "ZJetsTemplateProjection_Lep1Pt_e" + label + ".gif", doNormalize);
  DrawHist(histLep1Pt_m_fastsim, "ZJetsTemplateProjection_Lep1Pt_m" + label + ".gif", doNormalize);
  DrawHist(histLep2Pt_e_fastsim, "ZJetsTemplateProjection_Lep2Pt_e" + label + ".gif", doNormalize);
  DrawHist(histLep2Pt_m_fastsim, "ZJetsTemplateProjection_Lep2Pt_m" + label + ".gif", doNormalize);
  DrawHist(histLep3Pt_e_fastsim, "ZJetsTemplateProjection_Lep3Pt_e" + label + ".gif", doNormalize);
  DrawHist(histLep3Pt_m_fastsim, "ZJetsTemplateProjection_Lep3Pt_m" + label + ".gif", doNormalize);
  DrawHist(histLep4Pt_e_fastsim, "ZJetsTemplateProjection_Lep4Pt_e" + label + ".gif", doNormalize);
  DrawHist(histLep4Pt_m_fastsim, "ZJetsTemplateProjection_Lep4Pt_m" + label + ".gif", doNormalize);

  //*****************************************************************************************
  // Lepton Eta
  //*****************************************************************************************
  DrawHist(histLep1Eta_e_fastsim, "ZJetsTemplateProjection_Lep1Eta_e" + label + ".gif", doNormalize);
  DrawHist(histLep1Eta_m_fastsim, "ZJetsTemplateProjection_Lep1Eta_m" + label + ".gif", doNormalize);
  DrawHist(histLep2Eta_e_fastsim, "ZJetsTemplateProjection_Lep2Eta_e" + label + ".gif", doNormalize);
  DrawHist(histLep2Eta_m_fastsim, "ZJetsTemplateProjection_Lep2Eta_m" + label + ".gif", doNormalize);
  DrawHist(histLep3Eta_e_fastsim, "ZJetsTemplateProjection_Lep3Eta_e" + label + ".gif", doNormalize);
  DrawHist(histLep3Eta_m_fastsim, "ZJetsTemplateProjection_Lep3Eta_m" + label + ".gif", doNormalize);
  DrawHist(histLep4Eta_e_fastsim, "ZJetsTemplateProjection_Lep4Eta_e" + label + ".gif", doNormalize);
  DrawHist(histLep4Eta_m_fastsim, "ZJetsTemplateProjection_Lep4Eta_m" + label + ".gif", doNormalize);

  //*****************************************************************************************
  // Z1 Mass
  //*****************************************************************************************
  DrawHist(histZ1Mass_ee_fastsim, "ZJetsTemplateProjection_Z1Mass_ee" + label + ".gif", doNormalize);
  DrawHist(histZ1Mass_mm_fastsim, "ZJetsTemplateProjection_Z1Mass_mm" + label + ".gif", doNormalize);
  DrawHist(histZ2Mass_ee_fastsim, "ZJetsTemplateProjection_Z2Mass_ee" + label + ".gif", doNormalize);
  DrawHist(histZ2Mass_mm_fastsim, "ZJetsTemplateProjection_Z2Mass_mm" + label + ".gif", doNormalize);

  DrawHist(histZZMass_4e_fastsim, "ZJetsTemplateProjection_ZZMass_4e" + label + ".gif", doNormalize);
  DrawHist(histZZMass_4m_fastsim, "ZJetsTemplateProjection_ZZMass_4m" + label + ".gif", doNormalize);
  DrawHist(histZZMass_2e2m_fastsim, "ZJetsTemplateProjection_ZZMass_2e2m" + label + ".gif", doNormalize);
  DrawHist(histZZMass_2m2e_fastsim, "ZJetsTemplateProjection_ZZMass_2m2e" + label + ".gif", doNormalize);


  //*****************************************************************************************
  // Angles
  //*****************************************************************************************
  DrawHist(histPhi0_4e_fastsim, "ZJetsTemplateProjection_Phi0_4e" + label + ".gif", doNormalize);
  DrawHist(histPhi0_4m_fastsim, "ZJetsTemplateProjection_Phi0_4m" + label + ".gif", doNormalize);
  DrawHist(histPhi0_2e2m_fastsim, "ZJetsTemplateProjection_Phi0_2e2m" + label + ".gif", doNormalize);
  DrawHist(histPhi0_2m2e_fastsim, "ZJetsTemplateProjection_Phi0_2m2e" + label + ".gif", doNormalize);

  DrawHist(histTheta0_4e_fastsim, "ZJetsTemplateProjection_Theta0_4e" + label + ".gif", doNormalize);
  DrawHist(histTheta0_4m_fastsim, "ZJetsTemplateProjection_Theta0_4m" + label + ".gif", doNormalize);
  DrawHist(histTheta0_2e2m_fastsim, "ZJetsTemplateProjection_Theta0_2e2m" + label + ".gif", doNormalize);
  DrawHist(histTheta0_2m2e_fastsim, "ZJetsTemplateProjection_Theta0_2m2e" + label + ".gif", doNormalize);

  DrawHist(histPhi_4e_fastsim, "ZJetsTemplateProjection_Phi_4e" + label + ".gif", doNormalize);
  DrawHist(histPhi_4m_fastsim, "ZJetsTemplateProjection_Phi_4m" + label + ".gif", doNormalize);
  DrawHist(histPhi_2e2m_fastsim, "ZJetsTemplateProjection_Phi_2e2m" + label + ".gif", doNormalize);
  DrawHist(histPhi_2m2e_fastsim, "ZJetsTemplateProjection_Phi_2m2e" + label + ".gif", doNormalize);

  DrawHist(histTheta1_4e_fastsim, "ZJetsTemplateProjection_Theta1_4e" + label + ".gif", doNormalize);
  DrawHist(histTheta1_4m_fastsim, "ZJetsTemplateProjection_Theta1_4m" + label + ".gif", doNormalize);
  DrawHist(histTheta1_2e2m_fastsim, "ZJetsTemplateProjection_Theta1_2e2m" + label + ".gif", doNormalize);
  DrawHist(histTheta1_2m2e_fastsim, "ZJetsTemplateProjection_Theta1_2m2e" + label + ".gif", doNormalize);

  DrawHist(histTheta2_4e_fastsim, "ZJetsTemplateProjection_Theta2_4e" + label + ".gif", doNormalize);
  DrawHist(histTheta2_4m_fastsim, "ZJetsTemplateProjection_Theta2_4m" + label + ".gif", doNormalize);
  DrawHist(histTheta2_2e2m_fastsim, "ZJetsTemplateProjection_Theta2_2e2m" + label + ".gif", doNormalize);
  DrawHist(histTheta2_2m2e_fastsim, "ZJetsTemplateProjection_Theta2_2m2e" + label + ".gif", doNormalize);



  //*****************************************************************************************
  // Compare Yield
  //*****************************************************************************************
  cout << "Parameterized Model Yield : " << fastsimYield << "\n";

}

//root -l CMSAna/ObjectStudies/Photons/ComparePhotonDistributions.C+'("/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonStudies/ntuples/PhotonNtuple.HHtoBBGG-14tev-START53_V7A.root","/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonStudies/ntuples/PhotonNtuple.diphjets-START53_V7A.root","example","HHToBBGG","DiPhoton+Jets")'

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include "CMSAna/ObjectStudies/interface/PhotonTree.h"
#include <vector>
#include <map>


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

bool passSelection( double pt, double eta, int option) {

    bool passSelection = false;
    if (option == 0) {
      if (pt > 20 && pt <= 40 && fabs(eta) < 1.5) passSelection = true;
    }
    else if (option == 1) {
      if (pt > 40 && pt <= 60 && fabs(eta) < 1.5) passSelection = true;
    } 
    else if (option == 2) {
      if (pt > 60 && fabs(eta) < 1.5) passSelection = true;
    }
    else if (option == 10) {
      if (pt > 20 && pt <= 40 && fabs(eta) > 1.5) passSelection = true;
    } 
    else if (option == 11) {
      if (pt > 40 && pt <= 60 && fabs(eta) > 1.5) passSelection = true;
    } 
    else if (option == 12) {
      if (pt > 60 && fabs(eta) > 1.5) passSelection = true;
    }

    return passSelection;

}

//------------------------------------------------------------------------------
// PlotHiggsRes_LP
//------------------------------------------------------------------------------
void ComparePhotonDistributions ( string filename1 = "PhotonNtuple1.root",
                                  string filename2 = "PhotonNtuple2.root",
                                  string label = "",
                                  string legendLabel1 = "sample1",
                                  string legendLabel2 = "sample2",
                                  int option = -1
  )
{

  string Label = "";
  if (label != "") Label = "_" + label;

  //*******************************************************************************************
  //Define Histograms
  //*******************************************************************************************
 
       //*** variables _1****
  TH1F *histPhoHasPixelSeed_1 =  new TH1F( "histPhoHasPixelSeed_1", "; Has Pixel Seed; Number of Events", 2, -0.5, 1.5); 
  TH1F *histPhoIsConversion_1 =  new TH1F( "histPhoIsConversion_1", "; Is Conversion; Number of Events", 2, -0.5, 1.5); 
  TH1F *histPhoPassEleVeto_1 =  new TH1F( "histPhoPassEleVeto_1", "; Pass Ele Veto; Number of Events", 2, -0.5, 1.5); 
  TH1F *histPhoHoverE_1 =  new TH1F( "histPhoHoverE_1", "; H/E; Number of Events", 100, 0, 0.2); 
  TH1F *histPhoHoverESingleTower_1 =  new TH1F( "histPhoHoverESingleTower_1", "; H/E Single Tower; Number of Events", 100, 0, 0.2); 
  
       //***shower shape _1***
  TH1F *histPhotonSigmaIEtaIEta_1 =  new TH1F( "histPhotonSigmaIEtaIEta_1", "; #sigma_{i#eta i#eta};Number of Events", 100, 0, 0.05);  
  TH1F *histPhotonSigmaIPhiIPhi_1 =  new TH1F( "histPhotonSigmaIPhiIPhi_1", "; #sigma_{i#phi i#phi};Number of Events", 100, 0, 0.05);
  TH1F *histPhotonSigmaIEtaIPhi_1 =  new TH1F( "histPhotonSigmaIEtaIPhi_1", "; #sigma_{i#eta i#phi};Number of Events", 100, -1, 1); 
  TH1F *histPhoSCEtaWidth_1 =  new TH1F( "histPhoSCEtaWidth_1", "; SC Eta Width ;Number of Events", 100, 0, 0.03); 
  TH1F *histPhoSCPhiWidth_1 =  new TH1F( "histPhoSCPhiWidth_1", "; SC Phi Width ;Number of Events", 100, 0, 0.3); 
  TH1F *histPhoR9_1 =  new TH1F( "histPhoR9_1", "; R9 ;Number of Events", 100, 0.0, 1); 
  
       //***Isolation Variables _1***
  TH1F *histPhoTrkIso03_1 =  new TH1F( "histPhoTrkIso_1", "; Track Isolation 03; Number of Events", 100, 0, 10); 
  TH1F *histPhoEMIso03_1 =  new TH1F( "histPhoEMIso03_1", "; EM Isolation 03; Number of Events", 250, 0, 50); 
  TH1F *histPhoHadIso03_1 =  new TH1F( "histPhoHadIso03_1", "; Had Isolation 03; Number of Events", 150, 0, 30); 
  TH1F *histPhoPFIso03_1 =  new TH1F( "histPhoPFIso03_1", "; PF Isolation 03 / p_{T}; Number of Events", 100, 0, 10); 
  TH1F *histPhoPFIso04_1 =  new TH1F( "histPhoPFIso04_1", "; PF Isolation 04 / p_{T}; Number of Events", 100, 0, 10); 
  TH1F *histChargedIso_DR0p0To0p1_1 =  new TH1F( "histChargesIso_DR0p0To0p1_1", "; Charged Isolation DR0.0 to 0.1 / p_{T}; Number of Events", 100, 0, 10); 
  TH1F *histChargedIso_DR0p1To0p2_1 =  new TH1F( "histChargesIso_DR0p1To0p2_1", "; Charged Isolation DR0.1 to 0.2 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histChargedIso_DR0p2To0p3_1 =  new TH1F( "histChargesIso_DR0p2To0p3_1", "; Charged Isolation DR0.2 to 0.3 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histChargedIso_DR0p3To0p4_1 =  new TH1F( "histChargesIso_DR0p3To0p4_1", "; Charged Isolation DR0.3 to 0.4 / p_{T} / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histChargedIso_DR0p4To0p5_1 =  new TH1F( "histChargesIso_DR0p4To0p5_1", "; Charged Isolation DR0.4 to 0.5; Number of Events", 100, 0, 10);
  TH1F *histGammaIso_DR0p0To0p1_1 =  new TH1F( "histGammaIso_DR0p0To0p1_1", "; Gamma Isolation DR0.0 to 0.1 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histGammaIso_DR0p1To0p2_1 =  new TH1F( "histGammaIso_DR0p1To0p2_1", "; Gamma Isolation DR0.1 to 0.2 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histGammaIso_DR0p2To0p3_1 =  new TH1F( "histGammaIso_DR0p2To0p3_1", "; Gamma Isolation DR0.2 to 0.3 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histGammaIso_DR0p3To0p4_1 =  new TH1F( "histGammaIso_DR0p3To0p4_1", "; Gamma Isolation DR0.3 to 0.4 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histGammaIso_DR0p4To0p5_1 =  new TH1F( "histGammaIso_DR0p4To0p5_1", "; Gamma Isolation DR0.4 to 0.5 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histNeutralHadronIso_DR0p0To0p1_1 =  new TH1F( "histNeutralHadronIso_DR0p0To0p1_1", "; Neutral Hadron Isolation DR0.0 to 0.1 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histNeutralHadronIso_DR0p1To0p2_1 =  new TH1F( "histNeutralHadronIso_DR0p1To0p2_1", "; Neutral Hadron Isolation DR0.1 to 0.2 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histNeutralHadronIso_DR0p2To0p3_1 =  new TH1F( "histNeutralHadronIso_DR0p2To0p3_1", "; Neutral Hadron Isolation DR0.2 to 0.3 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histNeutralHadronIso_DR0p3To0p4_1 =  new TH1F( "histNeutralHadronIso_DR0p3To0p4_1", "; Neutral Hadron Isolation DR0.3 to 0.4 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histNeutralHadronIso_DR0p4To0p5_1 =  new TH1F( "histNeutralHadronIso_DR0p4To0p5_1", "; Neutral Hadron Isolation DR0.4 to 0.5 / p_{T}; Number of Events", 100, 0, 10);

       //***Regression Variables _1***
  TH1F *histSCRawEnergy_1 =  new TH1F( "histSCRawEnergy_1", "; SC Raw Energy; Number of Events", 100, 0, 100);
  TH1F *histPreShowerOverRaw_1 =  new TH1F( "histPreShowerOverRaw_1", "; Pre-shower/Raw; Number of Events", 100, 0, 1);
  TH1F *histNClusters_1 =  new TH1F( "histNClusters_1", "; N Clusters; Number of Events", 11, -0.5, 10.5);
  TH1F *histEtaSeed_1 =  new TH1F( "histEtaSeed_1", "; Eta Seed; Number of Events", 100, -3, 3);
  TH1F *histPhiSeed_1 =  new TH1F( "histPhiSeed_1", "; Phi Seed; Number of Events", 100, -3.2, 3.2);
  TH1F *histESeed_1 =  new TH1F( "histESeed_1", "; E Seed; Number of Events", 100, 0, 1);
  TH1F *histE3x3Seed_1 =  new TH1F( "histE3x3Seed_1", "; E3x3 Seed; Number of Events", 100, 0, 1);
  TH1F *histE5x5Seed_1 =  new TH1F( "histE5x5Seed_1", "; E5x5 Seed; Number of Events", 100, 0, 1);
  TH1F *histEMaxSeed_1 =  new TH1F( "histEMaxSeed_1", "; EMax Seed; Number of Events", 100, 0, 1);
  TH1F *histE2ndSeed_1 =  new TH1F( "histE2ndSeed_1", "; E2nd Seed; Number of Events", 100, 0, 1);
  TH1F *histETopSeed_1 =  new TH1F( "histETopSeed_1", "; E Top Seed; Number of Events", 100, 0, 1);
  TH1F *histEBottomSeed_1 =  new TH1F( "histEBottomSeed_1", "; E Bottom Seed; Number of Events", 100, 0, 1);
  TH1F *histELeftSeed_1 =  new TH1F( "histELeftSeed_1", "; E Left Seed; Number of Events", 100, 0, 1);
  TH1F *histERightSeed_1 =  new TH1F( "histERightSeed_1", "; E Right Seed; Number of Events", 100, 0, 1);
  TH1F *histE2x5MaxSeed_1 =  new TH1F( "histE2x5MaxSeed_1", "; E2x5Max Seed; Number of Events", 100, 0, 1);
  TH1F *histE2x5TopSeed_1 =  new TH1F( "histE2x5TopSeed_1", "; E2x5 Top Seed; Number of Events", 100, 0, 1);
  TH1F *histE2x5BottomSeed_1 =  new TH1F( "histE2x5BottomSeed_1", "; E2x5 Bottom Seed; Number of Events", 100, 0, 1);
  TH1F *histE2x5LeftSeed_1 =  new TH1F( "histE2x5LeftSeed_1", "; E2x5 Left Seed; Number of Events", 100, 0, 1);
  TH1F *histE2x5RightSeed_1 =  new TH1F( "histE2x5RightSeed_1", "; E2x5 Right Seed; Number of Events", 100, 0, 1);
  TH1F *histIEtaSeed_1 =  new TH1F( "histIEtaSeed_1", "; IEta Seed; Number of Events", 200, -99.5, 100.5);
  TH1F *histIPhiSeed_1 =  new TH1F( "histIPhiSeed_1", "; IPhi Seed; Number of Events", 360, -0.5, 359.5);
  TH1F *histEtaCrySeed_1 =  new TH1F( "histEtaCrySeed_1", "; Eta Cry  Seed; Number of Events", 100, -1, 1);
  TH1F *histPhiCrySeed_1 =  new TH1F( "histPhiCrySeed_1", "; Phi Cry  Seed; Number of Events", 100, -1, 1);
  TH1F *histGeneratedEnergy_1 =  new TH1F( "histGeneratedEnergy_1", "; Generated Energy; Number of Events", 100, 0, 1);
  TH1F *histPhoMomentum_Regression_V0_1 =  new TH1F( "histPhoMomentum_Regression_V0_1", "; Photon Momentum_Regression_V0; Number of Events", 100, 0, 1);
  TH1F *histPhoMomentumError_Regression_V0_1 =  new TH1F( "histPhoMomentumError_Regression_V0_1", "; Photon MomentumError_Regression_V0; Number of Events", 100, 0, 1);

       //*** variables _2**** 
  TH1F *histPhoHasPixelSeed_2 =  new TH1F( "histPhoHasPixelSeed_2", "; Has Pixel Seed; Number of Events", 2, -0.5, 1.5); 
  TH1F *histPhoIsConversion_2 =  new TH1F( "histPhoIsConversion_2", "; Is Conversion; Number of Events", 2, -0.5, 1.5); 
  TH1F *histPhoPassEleVeto_2 =  new TH1F( "histPhoPassEleVeto_2", "; Pass Ele Veto; Number of Events", 2, -0.5, 1.5); 
  TH1F *histPhoHoverE_2 =  new TH1F( "histPhoHoverE_2", "; H/E; Number of Events", 100, 0, 0.2); 
  TH1F *histPhoHoverESingleTower_2 =  new TH1F( "histPhoHoverESingleTower_2", "; H/E Single Tower; Number of Events", 100, 0, 0.2); 

       //***shower shape _2***
  TH1F *histPhotonSigmaIEtaIEta_2 =  new TH1F( "histPhotonSigmaIEtaIEta_2", "; #sigma_{i#eta i#eta};Number of Events", 100, 0, 0.05);
  TH1F *histPhotonSigmaIPhiIPhi_2 =  new TH1F( "histPhotonSigmaIPhiIPhi_2", "; #sigma_{i#phi i#phi};Number of Events", 100, 0, 0.05);
  TH1F *histPhotonSigmaIEtaIPhi_2 =  new TH1F( "histPhotonSigmaIEtaIPhi_2", "; #sigma_{i#eta i#phi};Number of Events", 100, -1, 1);
  TH1F *histPhoSCEtaWidth_2 =  new TH1F( "histPhoSCEtaWidth_2", "; SC Eta Width ;Number of Events", 100, 0, 0.03);
  TH1F *histPhoSCPhiWidth_2 =  new TH1F( "histPhoSCPhiWidth_2", "; SC Phi Width ;Number of Events", 100, 0, 0.3); 
  TH1F *histPhoR9_2 =  new TH1F( "histPhoR9_2", "; R9 ;Number of Events", 100, 0.0, 1); 
 
       //***Isolation Variables _2***
  TH1F *histPhoTrkIso03_2 =  new TH1F( "histPhoTrkIso_2", "; Track Isolation 03; Number of Events", 100, 0, 10); 
  TH1F *histPhoEMIso03_2 =  new TH1F( "histPhoEMIso03_2", "; EM Isolation 03; Number of Events", 250, 0, 50); 
  TH1F *histPhoHadIso03_2 =  new TH1F( "histPhoHadIso03_2", "; Had Isolation 03; Number of Events", 150, 0, 30); 
  TH1F *histPhoPFIso03_2 =  new TH1F( "histPhoPFIso03_2", "; PF Isolation 03 / p_{T}; Number of Events", 100, 0, 10); 
  TH1F *histPhoPFIso04_2 =  new TH1F( "histPhoPFIso04_2", "; PF Isolation 04 / p_{T}; Number of Events", 100, 0, 10); 
  TH1F *histChargedIso_DR0p0To0p1_2 =  new TH1F( "histChargesIso_DR0p0To0p1_2", "; Charged Isolation DR0.0 to 0.1 / p_{T}; Number of Events", 100, 0, 10); 
  TH1F *histChargedIso_DR0p1To0p2_2 =  new TH1F( "histChargesIso_DR0p1To0p2_2", "; Charged Isolation DR0.1 to 0.2 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histChargedIso_DR0p2To0p3_2 =  new TH1F( "histChargesIso_DR0p2To0p3_2", "; Charged Isolation DR0.2 to 0.3 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histChargedIso_DR0p3To0p4_2 =  new TH1F( "histChargesIso_DR0p3To0p4_2", "; Charged Isolation DR0.3 to 0.4 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histChargedIso_DR0p4To0p5_2 =  new TH1F( "histChargesIso_DR0p4To0p5_2", "; Charged Isolation DR0.4 to 0.5 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histGammaIso_DR0p0To0p1_2 =  new TH1F( "histGammaIso_DR0p0To0p1_2", "; Gamma Isolation DR0.0 to 0.1 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histGammaIso_DR0p1To0p2_2 =  new TH1F( "histGammaIso_DR0p1To0p2_2", "; Gamma Isolation DR0.1 to 0.2 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histGammaIso_DR0p2To0p3_2 =  new TH1F( "histGammaIso_DR0p2To0p3_2", "; Gamma Isolation DR0.2 to 0.3 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histGammaIso_DR0p3To0p4_2 =  new TH1F( "histGammaIso_DR0p3To0p4_2", "; Gamma Isolation DR0.3 to 0.4 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histGammaIso_DR0p4To0p5_2 =  new TH1F( "histGammaIso_DR0p4To0p5_2", "; Gamma Isolation DR0.4 to 0.5 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histNeutralHadronIso_DR0p0To0p1_2 =  new TH1F( "histNeutralHadronIso_DR0p0To0p1_2", "; Neutral Hadron Isolation DR0.0 to 0.1 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histNeutralHadronIso_DR0p1To0p2_2 =  new TH1F( "histNeutralHadronIso_DR0p1To0p2_2", "; Neutral Hadron Isolation DR0.1 to 0.2 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histNeutralHadronIso_DR0p2To0p3_2 =  new TH1F( "histNeutralHadronIso_DR0p2To0p3_2", "; Neutral Hadron Isolation DR0.2 to 0.3 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histNeutralHadronIso_DR0p3To0p4_2 =  new TH1F( "histNeutralHadronIso_DR0p3To0p4_2", "; Neutral Hadron Isolation DR0.3 to 0.4 / p_{T}; Number of Events", 100, 0, 10);
  TH1F *histNeutralHadronIso_DR0p4To0p5_2 =  new TH1F( "histNeutralHadronIso_DR0p4To0p5_2", "; Neutral Hadron Isolation DR0.4 to 0.5 / p_{T}; Number of Events", 100, 0, 10);

       //***Regression Variables _2***
  TH1F *histSCRawEnergy_2 =  new TH1F( "histSCRawEnergy_2", "; SC Raw Energy; Number of Events", 100, 0, 100);
  TH1F *histPreShowerOverRaw_2 =  new TH1F( "histPreShowerOverRaw_2", "; Pre-shower/Raw; Number of Events", 100, 0, 1);
  TH1F *histNClusters_2 =  new TH1F( "histNClusters_2", "; N Clusters; Number of Events", 11, -0.5, 10.5);
  TH1F *histEtaSeed_2 =  new TH1F( "histEtaSeed_2", "; Eta Seed; Number of Events", 100, -3, 3);
  TH1F *histPhiSeed_2 =  new TH1F( "histPhiSeed_2", "; Phi Seed; Number of Events", 100, -3.2, 3.2);
  TH1F *histESeed_2 =  new TH1F( "histESeed_2", "; E Seed; Number of Events", 100, 0, 1);
  TH1F *histE3x3Seed_2 =  new TH1F( "histE3x3Seed_2", "; E3x3 Seed; Number of Events", 100, 0, 1);
  TH1F *histE5x5Seed_2 =  new TH1F( "histE5x5Seed_2", "; E5x5 Seed; Number of Events", 100, 0, 1);
  TH1F *histEMaxSeed_2 =  new TH1F( "histEMaxSeed_2", "; EMax Seed; Number of Events", 100, 0, 1);
  TH1F *histE2ndSeed_2 =  new TH1F( "histE2ndSeed_2", "; E2nd Seed; Number of Events", 100, 0, 1);
  TH1F *histETopSeed_2 =  new TH1F( "histETopSeed_2", "; E Top Seed; Number of Events", 100, 0, 1);
  TH1F *histEBottomSeed_2 =  new TH1F( "histEBottomSeed_2", "; E Bottom Seed; Number of Events", 100, 0, 1);
  TH1F *histELeftSeed_2 =  new TH1F( "histELeftSeed_2", "; E Left Seed; Number of Events", 100, 0, 1);
  TH1F *histERightSeed_2 =  new TH1F( "histERightSeed_2", "; E Right Seed; Number of Events", 100, 0, 1);
  TH1F *histE2x5MaxSeed_2 =  new TH1F( "histE2x5MaxSeed_2", "; E2x5Max Seed; Number of Events", 100, 0, 1);
  TH1F *histE2x5TopSeed_2 =  new TH1F( "histE2x5TopSeed_2", "; E2x5 Top Seed; Number of Events", 100, 0, 1);
  TH1F *histE2x5BottomSeed_2 =  new TH1F( "histE2x5BottomSeed_2", "; E2x5 Bottom Seed; Number of Events", 100, 0, 1);
  TH1F *histE2x5LeftSeed_2 =  new TH1F( "histE2x5LeftSeed_2", "; E2x5 Left Seed; Number of Events", 100, 0, 1);
  TH1F *histE2x5RightSeed_2 =  new TH1F( "histE2x5RightSeed_2", "; E2x5 Right Seed; Number of Events", 100, 0, 1);
  TH1F *histIEtaSeed_2 =  new TH1F( "histIEtaSeed_2", "; IEta Seed; Number of Events", 200, -99.5, 100.5);
  TH1F *histIPhiSeed_2 =  new TH1F( "histIPhiSeed_2", "; IPhi Seed; Number of Events", 360, -0.5, 359.5);
  TH1F *histEtaCrySeed_2 =  new TH1F( "histEtaCrySeed_2", "; Eta Cry  Seed; Number of Events", 100, -1, 1);
  TH1F *histPhiCrySeed_2 =  new TH1F( "histPhiCrySeed_2", "; Phi Cry  Seed; Number of Events", 100, -1, 1);
  TH1F *histGeneratedEnergy_2 =  new TH1F( "histGeneratedEnergy_2", "; Generated Energy; Number of Events", 100, 0, 1);
  TH1F *histPhoMomentum_Regression_V0_2 =  new TH1F( "histPhoMomentum_Regression_V0_2", "; Photon Momentum_Regression_V0; Number of Events", 100, 0, 1);
  TH1F *histPhoMomentumError_Regression_V0_2 =  new TH1F( "histPhoMomentumError_Regression_V0_2", "; Photon MomentumError_Regression_V0; Number of Events", 100, 0, 1);


  //*******************************************************************************************
  //Read file1
  //*******************************************************************************************                
  cmsana::PhotonTree photree1;
  photree1.LoadTree(filename1.c_str());
  photree1.InitTree();

  for (int n=0;n<photree1.tree_->GetEntries();n++) { 
    photree1.tree_->GetEntry(n);

    if (!passSelection(photree1.fPhoPt,photree1.fPhoEta,option)) continue;

      //*** variables _1****
    histPhoHasPixelSeed_1->Fill(photree1.fPhoHasPixelSeed);
    histPhoIsConversion_1->Fill(photree1.fPhoIsConversion); 
    histPhoPassEleVeto_1->Fill(photree1.fPhoPassEleVeto); 
    histPhoHoverE_1->Fill(photree1.fPhoHoverE) ; 
    histPhoHoverESingleTower_1->Fill(photree1.fPhoHoverESingleTower);    

      //***shower shape _1****
    histPhotonSigmaIEtaIEta_1->Fill(photree1.fPhoSigmaIEtaIEta);
    histPhotonSigmaIPhiIPhi_1->Fill(photree1.fPhoSigmaIPhiIPhi);
    histPhotonSigmaIEtaIPhi_1->Fill(photree1.fPhoSigmaIEtaIPhi);
    histPhoSCEtaWidth_1->Fill(photree1.fPhoSCEtaWidth);
    histPhoSCPhiWidth_1->Fill(photree1.fPhoSCPhiWidth);
    histPhoR9_1->Fill(photree1.fPhoR9);

      //****Isolation Variables _1******
    histPhoTrkIso03_1->Fill(fmin(photree1.fPhoTrkIso03,9.99));
    histPhoEMIso03_1->Fill(fmin(photree1.fPhoEMIso03,49.99));
    histPhoHadIso03_1->Fill(fmin(photree1.fPhoHadIso03,29.99));
    histPhoPFIso03_1->Fill(fmin(photree1.fPhoPFIso03/photree1.fPhoPt,9.99));
    histPhoPFIso04_1->Fill(fmin(photree1.fPhoPFIso04/photree1.fPhoPt,9.99));
    histChargedIso_DR0p0To0p1_1->Fill(fmin(photree1.fChargedIso_DR0p0To0p1/photree1.fPhoPt,9.99));
    histChargedIso_DR0p1To0p2_1->Fill(fmin(photree1.fChargedIso_DR0p1To0p2/photree1.fPhoPt,9.99));
    histChargedIso_DR0p2To0p3_1->Fill(fmin(photree1.fChargedIso_DR0p2To0p3/photree1.fPhoPt,9.99));
    histChargedIso_DR0p3To0p4_1->Fill(fmin(photree1.fChargedIso_DR0p3To0p4/photree1.fPhoPt,9.99));
    histChargedIso_DR0p4To0p5_1->Fill(fmin(photree1.fChargedIso_DR0p4To0p5/photree1.fPhoPt,9.99));
    histGammaIso_DR0p0To0p1_1->Fill(fmin(photree1.fGammaIso_DR0p0To0p1/photree1.fPhoPt,9.99));
    histGammaIso_DR0p1To0p2_1->Fill(fmin(photree1.fGammaIso_DR0p1To0p2/photree1.fPhoPt,9.99));
    histGammaIso_DR0p2To0p3_1->Fill(fmin(photree1.fGammaIso_DR0p2To0p3/photree1.fPhoPt,9.99));
    histGammaIso_DR0p3To0p4_1->Fill(fmin(photree1.fGammaIso_DR0p3To0p4/photree1.fPhoPt,9.99));
    histGammaIso_DR0p4To0p5_1->Fill(fmin(photree1.fGammaIso_DR0p4To0p5/photree1.fPhoPt,9.99));
    histNeutralHadronIso_DR0p0To0p1_1->Fill(fmin(photree1.fNeutralHadronIso_DR0p0To0p1/photree1.fPhoPt,9.99));
    histNeutralHadronIso_DR0p1To0p2_1->Fill(fmin(photree1.fNeutralHadronIso_DR0p1To0p2/photree1.fPhoPt,9.99));
    histNeutralHadronIso_DR0p2To0p3_1->Fill(fmin(photree1.fNeutralHadronIso_DR0p2To0p3/photree1.fPhoPt,9.99));
    histNeutralHadronIso_DR0p3To0p4_1->Fill(fmin(photree1.fNeutralHadronIso_DR0p3To0p4/photree1.fPhoPt,9.99));
    histNeutralHadronIso_DR0p4To0p5_1->Fill(fmin(photree1.fNeutralHadronIso_DR0p4To0p5/photree1.fPhoPt,9.99));

      //****Regression Variables _1****
    histSCRawEnergy_1->Fill(photree1.fSCRawEnergy);
    histPreShowerOverRaw_1->Fill(photree1.fPreShowerOverRaw);
    histNClusters_1->Fill(photree1.fNClusters);
    histEtaSeed_1->Fill(photree1.fEtaSeed);
    histPhiSeed_1->Fill(photree1.fPhiSeed);
    histESeed_1->Fill(photree1.fESeed/photree1.fSCRawEnergy);
    histE3x3Seed_1->Fill(photree1.fE3x3Seed/photree1.fSCRawEnergy);
    histE5x5Seed_1->Fill(photree1.fE5x5Seed/photree1.fSCRawEnergy);
    histEMaxSeed_1->Fill(photree1.fEMaxSeed/photree1.fSCRawEnergy);
    histE2ndSeed_1->Fill(photree1.fE2ndSeed/photree1.fSCRawEnergy);
    histETopSeed_1->Fill(photree1.fETopSeed/photree1.fSCRawEnergy);
    histEBottomSeed_1->Fill(photree1.fEBottomSeed/photree1.fSCRawEnergy);
    histELeftSeed_1->Fill(photree1.fELeftSeed/photree1.fSCRawEnergy);
    histERightSeed_1->Fill(photree1.fERightSeed/photree1.fSCRawEnergy);
    histE2x5MaxSeed_1->Fill(photree1.fE2x5MaxSeed/photree1.fSCRawEnergy);
    histE2x5TopSeed_1->Fill(photree1.fE2x5TopSeed/photree1.fSCRawEnergy);
    histE2x5BottomSeed_1->Fill(photree1.fE2x5BottomSeed/photree1.fSCRawEnergy);
    histE2x5LeftSeed_1->Fill(photree1.fE2x5LeftSeed/photree1.fSCRawEnergy);
    histE2x5RightSeed_1->Fill(photree1.fE2x5RightSeed/photree1.fSCRawEnergy);
    histIEtaSeed_1->Fill(photree1.fIEtaSeed);
    histIPhiSeed_1->Fill(photree1.fIPhiSeed);
    histEtaCrySeed_1->Fill(photree1.fEtaCrySeed);
    histPhiCrySeed_1->Fill(photree1.fPhiCrySeed);
    histGeneratedEnergy_1->Fill(photree1.fGeneratedEnergy);
    histPhoMomentum_Regression_V0_1->Fill(photree1.fPhoMomentum_Regression_V0);
    histPhoMomentumError_Regression_V0_1->Fill(photree1.fPhoMomentumError_Regression_V0);

  }
  
  //*******************************************************************************************
  //Read file2
  //*******************************************************************************************                
  cmsana::PhotonTree photree2;
  photree2.LoadTree(filename2.c_str());
  photree2.InitTree();

  for (int n=0;n<photree2.tree_->GetEntries();n++) { 
    photree2.tree_->GetEntry(n);

    if (!passSelection(photree2.fPhoPt,photree2.fPhoEta,option)) continue;

      //*** variables _2****
    histPhoHasPixelSeed_2->Fill(photree2.fPhoHasPixelSeed);
    histPhoIsConversion_2->Fill(photree2.fPhoIsConversion); 
    histPhoPassEleVeto_2->Fill(photree2.fPhoPassEleVeto); 
    histPhoHoverE_2->Fill(photree2.fPhoHoverE) ; 
    histPhoHoverESingleTower_2->Fill(photree2.fPhoHoverESingleTower);    

      //***shower shape _2****
    histPhotonSigmaIEtaIEta_2->Fill(photree2.fPhoSigmaIEtaIEta);
    histPhotonSigmaIPhiIPhi_2->Fill(photree2.fPhoSigmaIPhiIPhi);
    histPhotonSigmaIEtaIPhi_2->Fill(photree2.fPhoSigmaIEtaIPhi);
    histPhoSCEtaWidth_2->Fill(photree2.fPhoSCEtaWidth);
    histPhoSCPhiWidth_2->Fill(photree2.fPhoSCPhiWidth);
    histPhoR9_2->Fill(photree2.fPhoR9);

      //****Isolation Variables _2******
    histPhoTrkIso03_2->Fill(fmin(photree2.fPhoTrkIso03,9.99));
    histPhoEMIso03_2->Fill(fmin(photree2.fPhoEMIso03,49.99));
    histPhoHadIso03_2->Fill(fmin(photree2.fPhoHadIso03,29.99));
    histPhoPFIso03_2->Fill(fmin(photree2.fPhoPFIso03/photree1.fPhoPt,9.99));
    histPhoPFIso04_2->Fill(fmin(photree2.fPhoPFIso04/photree1.fPhoPt,9.99));
    histChargedIso_DR0p0To0p1_2->Fill(fmin(photree2.fChargedIso_DR0p0To0p1/photree2.fPhoPt,9.99));
    histChargedIso_DR0p1To0p2_2->Fill(fmin(photree2.fChargedIso_DR0p1To0p2/photree2.fPhoPt,9.99));
    histChargedIso_DR0p2To0p3_2->Fill(fmin(photree2.fChargedIso_DR0p2To0p3/photree2.fPhoPt,9.99));
    histChargedIso_DR0p3To0p4_2->Fill(fmin(photree2.fChargedIso_DR0p3To0p4/photree2.fPhoPt,9.99));
    histChargedIso_DR0p4To0p5_2->Fill(fmin(photree2.fChargedIso_DR0p4To0p5/photree2.fPhoPt,9.99));
    histGammaIso_DR0p0To0p1_2->Fill(fmin(photree2.fGammaIso_DR0p0To0p1/photree2.fPhoPt,9.99));
    histGammaIso_DR0p1To0p2_2->Fill(fmin(photree2.fGammaIso_DR0p1To0p2/photree2.fPhoPt,9.99));
    histGammaIso_DR0p2To0p3_2->Fill(fmin(photree2.fGammaIso_DR0p2To0p3/photree2.fPhoPt,9.99));
    histGammaIso_DR0p3To0p4_2->Fill(fmin(photree2.fGammaIso_DR0p3To0p4/photree2.fPhoPt,9.99));
    histGammaIso_DR0p4To0p5_2->Fill(fmin(photree2.fGammaIso_DR0p4To0p5/photree2.fPhoPt,9.99));
    histNeutralHadronIso_DR0p0To0p1_2->Fill(fmin(photree2.fNeutralHadronIso_DR0p0To0p1/photree2.fPhoPt,9.99));
    histNeutralHadronIso_DR0p1To0p2_2->Fill(fmin(photree2.fNeutralHadronIso_DR0p1To0p2/photree2.fPhoPt,9.99));
    histNeutralHadronIso_DR0p2To0p3_2->Fill(fmin(photree2.fNeutralHadronIso_DR0p2To0p3/photree2.fPhoPt,9.99));
    histNeutralHadronIso_DR0p3To0p4_2->Fill(fmin(photree2.fNeutralHadronIso_DR0p3To0p4/photree2.fPhoPt,9.99));
    histNeutralHadronIso_DR0p4To0p5_2->Fill(fmin(photree2.fNeutralHadronIso_DR0p4To0p5/photree2.fPhoPt,9.99));

      //****Regression Variables _2****
    histSCRawEnergy_2->Fill(photree2.fSCRawEnergy);
    histPreShowerOverRaw_2->Fill(photree2.fPreShowerOverRaw);
    histNClusters_2->Fill(photree2.fNClusters);
    histEtaSeed_2->Fill(photree2.fEtaSeed);
    histPhiSeed_2->Fill(photree2.fPhiSeed);
    histESeed_2->Fill(photree2.fESeed/photree2.fSCRawEnergy);
    histE3x3Seed_2->Fill(photree2.fE3x3Seed/photree2.fSCRawEnergy);
    histE5x5Seed_2->Fill(photree2.fE5x5Seed/photree2.fSCRawEnergy);
    histEMaxSeed_2->Fill(photree2.fEMaxSeed/photree2.fSCRawEnergy);
    histE2ndSeed_2->Fill(photree2.fE2ndSeed/photree2.fSCRawEnergy);
    histETopSeed_2->Fill(photree2.fETopSeed/photree2.fSCRawEnergy);
    histEBottomSeed_2->Fill(photree2.fEBottomSeed/photree2.fSCRawEnergy);
    histELeftSeed_2->Fill(photree2.fELeftSeed/photree2.fSCRawEnergy);
    histERightSeed_2->Fill(photree2.fERightSeed/photree2.fSCRawEnergy);
    histE2x5MaxSeed_2->Fill(photree2.fE2x5MaxSeed/photree2.fSCRawEnergy);
    histE2x5TopSeed_2->Fill(photree2.fE2x5TopSeed/photree2.fSCRawEnergy);
    histE2x5BottomSeed_2->Fill(photree2.fE2x5BottomSeed/photree2.fSCRawEnergy);
    histE2x5LeftSeed_2->Fill(photree2.fE2x5LeftSeed/photree2.fSCRawEnergy);
    histE2x5RightSeed_2->Fill(photree2.fE2x5RightSeed/photree2.fSCRawEnergy);
    histIEtaSeed_2->Fill(photree2.fIEtaSeed);
    histIPhiSeed_2->Fill(photree2.fIPhiSeed);
    histEtaCrySeed_2->Fill(photree2.fEtaCrySeed);
    histPhiCrySeed_2->Fill(photree2.fPhiCrySeed);
    histGeneratedEnergy_2->Fill(photree2.fGeneratedEnergy);
    histPhoMomentum_Regression_V0_2->Fill(photree2.fPhoMomentum_Regression_V0);
    histPhoMomentumError_Regression_V0_2->Fill(photree2.fPhoMomentumError_Regression_V0);

  }
  
  //Draw Plots
  TCanvas *cv = 0;
  TLegend *legend = 0;

  //*******************************************************************************************
  //Normalize histograms
  //*******************************************************************************************
 
     //*****Variables _1********
  NormalizeHist(histPhoHasPixelSeed_1);
  NormalizeHist(histPhoIsConversion_1);
  NormalizeHist(histPhoPassEleVeto_1);
  NormalizeHist(histPhoHoverE_1);
  NormalizeHist(histPhoHoverESingleTower_1);

     //*****shower shap _1********
  NormalizeHist(histPhotonSigmaIEtaIEta_1);
  NormalizeHist(histPhotonSigmaIPhiIPhi_1);
  NormalizeHist(histPhotonSigmaIEtaIPhi_1);
  NormalizeHist(histPhoSCEtaWidth_1);
  NormalizeHist(histPhoSCPhiWidth_1);
  NormalizeHist(histPhoR9_1);

     //****Isolation Variables _1*****
  NormalizeHist(histPhoTrkIso03_1);
  NormalizeHist(histPhoEMIso03_1);
  NormalizeHist(histPhoHadIso03_1);
  NormalizeHist(histPhoPFIso03_1);
  NormalizeHist(histPhoPFIso04_1);
  NormalizeHist(histChargedIso_DR0p0To0p1_1);
  NormalizeHist(histChargedIso_DR0p1To0p2_1);
  NormalizeHist(histChargedIso_DR0p2To0p3_1);
  NormalizeHist(histChargedIso_DR0p3To0p4_1);
  NormalizeHist(histChargedIso_DR0p4To0p5_1);
  NormalizeHist(histGammaIso_DR0p0To0p1_1);
  NormalizeHist(histGammaIso_DR0p1To0p2_1);
  NormalizeHist(histGammaIso_DR0p2To0p3_1);
  NormalizeHist(histGammaIso_DR0p3To0p4_1);
  NormalizeHist(histGammaIso_DR0p4To0p5_1);
  NormalizeHist(histNeutralHadronIso_DR0p0To0p1_1);
  NormalizeHist(histNeutralHadronIso_DR0p1To0p2_1);
  NormalizeHist(histNeutralHadronIso_DR0p2To0p3_1);
  NormalizeHist(histNeutralHadronIso_DR0p3To0p4_1);
  NormalizeHist(histNeutralHadronIso_DR0p4To0p5_1);

  //****Regression Variables _1****
  NormalizeHist(histSCRawEnergy_1);
  NormalizeHist(histPreShowerOverRaw_1); 
  NormalizeHist(histNClusters_1);
  NormalizeHist(histEtaSeed_1);
  NormalizeHist(histPhiSeed_1);
  NormalizeHist(histESeed_1);
  NormalizeHist(histE3x3Seed_1);
  NormalizeHist(histE5x5Seed_1);
  NormalizeHist(histEMaxSeed_1);
  NormalizeHist(histE2ndSeed_1);
  NormalizeHist(histETopSeed_1);
  NormalizeHist(histEBottomSeed_1);
  NormalizeHist(histELeftSeed_1);
  NormalizeHist(histERightSeed_1);
  NormalizeHist(histE2x5MaxSeed_1);
  NormalizeHist(histE2x5TopSeed_1);
  NormalizeHist(histE2x5BottomSeed_1);
  NormalizeHist(histE2x5LeftSeed_1);
  NormalizeHist(histE2x5RightSeed_1);
  NormalizeHist(histIEtaSeed_1);
  NormalizeHist(histIPhiSeed_1);
  NormalizeHist(histEtaCrySeed_1);
  NormalizeHist(histPhiCrySeed_1);
  NormalizeHist(histGeneratedEnergy_1);
  NormalizeHist(histPhoMomentum_Regression_V0_1);
  NormalizeHist(histPhoMomentumError_Regression_V0_1);

     //****Variables 2*******
  NormalizeHist(histPhoHasPixelSeed_2);
  NormalizeHist(histPhoIsConversion_2);
  NormalizeHist(histPhoPassEleVeto_2);
  NormalizeHist(histPhoHoverE_2);
  NormalizeHist(histPhoHoverESingleTower_2);

     //*****shower shap _2********
  NormalizeHist(histPhotonSigmaIEtaIEta_2);
  NormalizeHist(histPhotonSigmaIPhiIPhi_2);
  NormalizeHist(histPhotonSigmaIEtaIPhi_2);
  NormalizeHist(histPhoSCEtaWidth_2);
  NormalizeHist(histPhoSCPhiWidth_2);
  NormalizeHist(histPhoR9_2);

     //****Isolation Variables _2*****
  NormalizeHist(histPhoTrkIso03_2);
  NormalizeHist(histPhoEMIso03_2);
  NormalizeHist(histPhoHadIso03_2);
  NormalizeHist(histPhoPFIso03_2);
  NormalizeHist(histPhoPFIso04_2);
  NormalizeHist(histChargedIso_DR0p0To0p1_2);
  NormalizeHist(histChargedIso_DR0p1To0p2_2);
  NormalizeHist(histChargedIso_DR0p2To0p3_2);
  NormalizeHist(histChargedIso_DR0p3To0p4_2);
  NormalizeHist(histChargedIso_DR0p4To0p5_2);
  NormalizeHist(histGammaIso_DR0p0To0p1_2);
  NormalizeHist(histGammaIso_DR0p1To0p2_2);
  NormalizeHist(histGammaIso_DR0p2To0p3_2);
  NormalizeHist(histGammaIso_DR0p3To0p4_2);
  NormalizeHist(histGammaIso_DR0p4To0p5_2);
  NormalizeHist(histNeutralHadronIso_DR0p0To0p1_2);
  NormalizeHist(histNeutralHadronIso_DR0p1To0p2_2);
  NormalizeHist(histNeutralHadronIso_DR0p2To0p3_2);
  NormalizeHist(histNeutralHadronIso_DR0p3To0p4_2);
  NormalizeHist(histNeutralHadronIso_DR0p4To0p5_2);

  //****Regression Variables _2****
  NormalizeHist(histSCRawEnergy_2);
  NormalizeHist(histPreShowerOverRaw_2); 
  NormalizeHist(histNClusters_2);
  NormalizeHist(histEtaSeed_2);
  NormalizeHist(histPhiSeed_2);
  NormalizeHist(histESeed_2);
  NormalizeHist(histE3x3Seed_2);
  NormalizeHist(histE5x5Seed_2);
  NormalizeHist(histEMaxSeed_2);
  NormalizeHist(histE2ndSeed_2);
  NormalizeHist(histETopSeed_2);
  NormalizeHist(histEBottomSeed_2);
  NormalizeHist(histELeftSeed_2);
  NormalizeHist(histERightSeed_2);
  NormalizeHist(histE2x5MaxSeed_2);
  NormalizeHist(histE2x5TopSeed_2);
  NormalizeHist(histE2x5BottomSeed_2);
  NormalizeHist(histE2x5LeftSeed_2);
  NormalizeHist(histE2x5RightSeed_2);
  NormalizeHist(histIEtaSeed_2);
  NormalizeHist(histIPhiSeed_2);
  NormalizeHist(histEtaCrySeed_2);
  NormalizeHist(histPhiCrySeed_2);
  NormalizeHist(histGeneratedEnergy_2);
  NormalizeHist(histPhoMomentum_Regression_V0_2);
  NormalizeHist(histPhoMomentumError_Regression_V0_2);

  //*******************************************************************************************
  //Photon Has Pixel Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhoHasPixelSeed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhoHasPixelSeed_2, legendLabel2.c_str(), "L");
  histPhoHasPixelSeed_1->SetLineColor(kBlue);
  histPhoHasPixelSeed_2->SetLineColor(kRed);
  histPhoHasPixelSeed_1->Draw("hist");
  histPhoHasPixelSeed_2->Draw("hist, same");
  histPhoHasPixelSeed_1->SetMinimum(0.0);
  histPhoHasPixelSeed_1->SetMaximum(1.0);
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_HasPixelSeed"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon Is Conversion
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhoIsConversion_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhoIsConversion_2, legendLabel2.c_str(), "L");
  histPhoIsConversion_1->SetLineColor(kBlue);
  histPhoIsConversion_2->SetLineColor(kRed);
  histPhoIsConversion_1->Draw("hist");
  histPhoIsConversion_2->Draw("hist,same");
  histPhoIsConversion_1->SetMinimum(0.0);
  histPhoIsConversion_1->SetMaximum(1.0);
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_IsConversion"+Label+".gif").c_str() );

//*******************************************************************************************
  //Photon Pass EleVeto
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhoPassEleVeto_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhoPassEleVeto_2, legendLabel2.c_str(), "L");
  histPhoPassEleVeto_1->SetLineColor(kBlue);
  histPhoPassEleVeto_2->SetLineColor(kRed);
  histPhoPassEleVeto_1->Draw("hist");
  histPhoPassEleVeto_2->Draw("hist, same");
  histPhoPassEleVeto_1->SetMinimum(0.0);
  histPhoPassEleVeto_1->SetMaximum(1.0);
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_PassEleVeto"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon H over E
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhoHoverE_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhoHoverE_2, legendLabel2.c_str(), "L");
  histPhoHoverE_1->SetLineColor(kBlue);
  histPhoHoverE_2->SetLineColor(kRed);
  histPhoHoverE_1->Draw("hist");
  histPhoHoverE_2->Draw("hist, same");
  legend->Draw();
  cv->SetLogy();
  cv->SaveAs(("PhotonDistribution_HoverE"+Label+".gif").c_str() );

 //*******************************************************************************************
  //Photon H over E single tower
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhoHoverESingleTower_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhoHoverESingleTower_2, legendLabel2.c_str(), "L");
  histPhoHoverESingleTower_1->SetLineColor(kBlue);
  histPhoHoverESingleTower_2->SetLineColor(kRed);
  histPhoHoverESingleTower_1->Draw("hist");
  histPhoHoverESingleTower_2->Draw("hist, same");
  legend->Draw();
  cv->SetLogy();
  cv->SaveAs(("PhotonDistribution_HoverESingleTower"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon SigmaIEtaIEta
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhotonSigmaIEtaIEta_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhotonSigmaIEtaIEta_2, legendLabel2.c_str(), "L");
  histPhotonSigmaIEtaIEta_1->SetLineColor(kBlue);
  histPhotonSigmaIEtaIEta_2->SetLineColor(kRed);
  histPhotonSigmaIEtaIEta_1->Draw("hist");
  histPhotonSigmaIEtaIEta_2->Draw("hist, same");
  histPhotonSigmaIEtaIEta_1->SetMaximum( 1.2 * fmax(histPhotonSigmaIEtaIEta_1->GetMaximum(),histPhotonSigmaIEtaIEta_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_SigmaIEtaIEta"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon SigmaIPhiIPhi
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhotonSigmaIPhiIPhi_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhotonSigmaIPhiIPhi_2, legendLabel2.c_str(), "L");
  histPhotonSigmaIPhiIPhi_1->SetLineColor(kBlue);
  histPhotonSigmaIPhiIPhi_2->SetLineColor(kRed);
  histPhotonSigmaIPhiIPhi_1->Draw("hist");
  histPhotonSigmaIPhiIPhi_2->Draw("hist, same");
  histPhotonSigmaIPhiIPhi_1->SetMaximum( 1.2 * fmax(histPhotonSigmaIPhiIPhi_1->GetMaximum(),histPhotonSigmaIPhiIPhi_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_SigmaIPhiIPhi"+Label+".gif").c_str() );
 
  //*******************************************************************************************
  //Photon SigmaIEtaIPhi
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhotonSigmaIEtaIPhi_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhotonSigmaIEtaIPhi_2, legendLabel2.c_str(), "L");
  histPhotonSigmaIEtaIPhi_1->SetLineColor(kBlue);
  histPhotonSigmaIEtaIPhi_2->SetLineColor(kRed);
  histPhotonSigmaIEtaIPhi_1->Draw("hist");
  histPhotonSigmaIEtaIPhi_2->Draw("hist, same");
  histPhotonSigmaIEtaIPhi_1->SetMaximum( 1.2 * fmax(histPhotonSigmaIEtaIPhi_1->GetMaximum(),histPhotonSigmaIEtaIPhi_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_SigmaIEtaIPhi"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon SC Eta Width
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhoSCEtaWidth_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhoSCEtaWidth_2, legendLabel2.c_str(), "L");
  histPhoSCEtaWidth_1->SetLineColor(kBlue);
  histPhoSCEtaWidth_2->SetLineColor(kRed);
  histPhoSCEtaWidth_1->Draw("hist");
  histPhoSCEtaWidth_2->Draw("hist, same");
  histPhoSCEtaWidth_1->SetMaximum( 1.2 * fmax(histPhoSCEtaWidth_1->GetMaximum(),histPhoSCEtaWidth_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_SCEtaWidth"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon SC Phi Width
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhoSCPhiWidth_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhoSCPhiWidth_2, legendLabel2.c_str(), "L");
  histPhoSCPhiWidth_1->SetLineColor(kBlue);
  histPhoSCPhiWidth_2->SetLineColor(kRed);
  histPhoSCPhiWidth_1->Draw("hist");
  histPhoSCPhiWidth_2->Draw("hist, same");
  histPhoSCPhiWidth_1->SetMaximum( 1.2 * fmax(histPhoSCPhiWidth_1->GetMaximum(),histPhoSCPhiWidth_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_SCPhiWidth"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon R9
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhoR9_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhoR9_2, legendLabel2.c_str(), "L");
  histPhoR9_1->SetLineColor(kBlue);
  histPhoR9_2->SetLineColor(kRed);
  histPhoR9_1->Draw("hist");
  histPhoR9_2->Draw("hist, same");
  histPhoR9_1->SetMaximum( 1.2 * fmax(histPhoR9_1->GetMaximum(),histPhoR9_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_R9"+Label+".gif").c_str() );

 //*******************************************************************************************
  //Photon Track Isolation 03
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhoTrkIso03_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhoTrkIso03_2, legendLabel2.c_str(), "L");
  histPhoTrkIso03_1->SetLineColor(kBlue);
  histPhoTrkIso03_2->SetLineColor(kRed);
  histPhoTrkIso03_1->Draw("hist");
  histPhoTrkIso03_2->Draw("hist, same");
  histPhoTrkIso03_1->SetMaximum( 1.2 * fmax(histPhoTrkIso03_1->GetMaximum(),histPhoTrkIso03_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_TrkIsolation03"+Label+".gif").c_str() );

 //*******************************************************************************************
  //Photon EM Isolation 03
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhoEMIso03_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhoEMIso03_2, legendLabel2.c_str(), "L");
  histPhoEMIso03_1->SetLineColor(kBlue);
  histPhoEMIso03_2->SetLineColor(kRed);
  histPhoEMIso03_1->Draw("hist");
  histPhoEMIso03_2->Draw("hist, same");
  histPhoEMIso03_1->SetMaximum( 1.2 * fmax(histPhoEMIso03_1->GetMaximum(),histPhoEMIso03_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_EMIsolation03"+Label+".gif").c_str() );


 //*******************************************************************************************
  //Photon Had Isolation 03
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhoHadIso03_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhoHadIso03_2, legendLabel2.c_str(), "L");
  histPhoHadIso03_1->SetLineColor(kBlue);
  histPhoHadIso03_2->SetLineColor(kRed);
  histPhoHadIso03_1->Draw("hist");
  histPhoHadIso03_2->Draw("hist, same");
  histPhoHadIso03_1->SetMaximum( 1.2 * fmax(histPhoHadIso03_1->GetMaximum(),histPhoHadIso03_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_HadIsolation03"+Label+".gif").c_str() );
 //*******************************************************************************************
  //Photon PF Isolation 03
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhoPFIso03_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhoPFIso03_2, legendLabel2.c_str(), "L");
  histPhoPFIso03_1->SetLineColor(kBlue);
  histPhoPFIso03_2->SetLineColor(kRed);
  histPhoPFIso03_1->Draw("hist");
  histPhoPFIso03_2->Draw("hist, same");
  histPhoPFIso03_1->SetMaximum( 1.2 * fmax(histPhoPFIso03_1->GetMaximum(),histPhoPFIso03_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_PFIso03"+Label+".gif").c_str() );


 //*******************************************************************************************
  //Photon PF Isolation 04
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhoPFIso04_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhoPFIso04_2, legendLabel2.c_str(), "L");
  histPhoPFIso04_1->SetLineColor(kBlue);
  histPhoPFIso04_2->SetLineColor(kRed);
  histPhoPFIso04_1->Draw("hist");
  histPhoPFIso04_2->Draw("hist, same");
  histPhoPFIso04_1->SetMaximum( 1.2 * fmax(histPhoPFIso04_1->GetMaximum(),histPhoPFIso04_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_PFIso04"+Label+".gif").c_str() );


  //*******************************************************************************************
  //Photon Charged Isolation DR 0.0 to 0.1
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histChargedIso_DR0p0To0p1_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histChargedIso_DR0p0To0p1_2, legendLabel2.c_str(), "L");
  histChargedIso_DR0p0To0p1_1->SetLineColor(kBlue);
  histChargedIso_DR0p0To0p1_2->SetLineColor(kRed);
  histChargedIso_DR0p0To0p1_1->Draw("hist");
  histChargedIso_DR0p0To0p1_2->Draw("hist, same");
  histChargedIso_DR0p0To0p1_1->SetMaximum( 1.2 * fmax(histChargedIso_DR0p0To0p1_1->GetMaximum(),histChargedIso_DR0p0To0p1_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("ChargedIso_DR0p0To0p1"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon Charged Isolation DR 0.1 to 0.2
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histChargedIso_DR0p1To0p2_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histChargedIso_DR0p1To0p2_2, legendLabel2.c_str(), "L");
  histChargedIso_DR0p1To0p2_1->SetLineColor(kBlue);
  histChargedIso_DR0p1To0p2_2->SetLineColor(kRed);
  histChargedIso_DR0p1To0p2_1->Draw("hist");
  histChargedIso_DR0p1To0p2_2->Draw("hist, same");
  histChargedIso_DR0p1To0p2_1->SetMaximum( 1.2 * fmax(histChargedIso_DR0p1To0p2_1->GetMaximum(),histChargedIso_DR0p1To0p2_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("ChargedIso_DR0p1To0p2"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon Charged Isolation DR 0.2 To 0.3
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histChargedIso_DR0p2To0p3_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histChargedIso_DR0p2To0p3_2, legendLabel2.c_str(), "L");
  histChargedIso_DR0p2To0p3_1->SetLineColor(kBlue);
  histChargedIso_DR0p2To0p3_2->SetLineColor(kRed);
  histChargedIso_DR0p2To0p3_1->Draw("hist");
  histChargedIso_DR0p2To0p3_2->Draw("hist, same");
  histChargedIso_DR0p2To0p3_1->SetMaximum( 1.2 * fmax(histChargedIso_DR0p2To0p3_1->GetMaximum(),histChargedIso_DR0p2To0p3_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("ChargedIso_DR0p2To0p3"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon Charged Isolation DR 0.3 To 0.4
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histChargedIso_DR0p3To0p4_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histChargedIso_DR0p3To0p4_2, legendLabel2.c_str(), "L");
  histChargedIso_DR0p3To0p4_1->SetLineColor(kBlue);
  histChargedIso_DR0p3To0p4_2->SetLineColor(kRed);
  histChargedIso_DR0p3To0p4_1->Draw("hist");
  histChargedIso_DR0p3To0p4_2->Draw("hist, same");
  histChargedIso_DR0p3To0p4_1->SetMaximum( 1.2 * fmax(histChargedIso_DR0p3To0p4_1->GetMaximum(),histChargedIso_DR0p3To0p4_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("ChargedIso_DR0p3To0p4"+Label+".gif").c_str() );

 //*******************************************************************************************
  //Photon Charged Isolation 0.4 To 0.5
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histChargedIso_DR0p4To0p5_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histChargedIso_DR0p4To0p5_2, legendLabel2.c_str(), "L");
  histChargedIso_DR0p4To0p5_1->SetLineColor(kBlue);
  histChargedIso_DR0p4To0p5_2->SetLineColor(kRed);
  histChargedIso_DR0p4To0p5_1->Draw("hist");
  histChargedIso_DR0p4To0p5_2->Draw("hist, same");
  histChargedIso_DR0p4To0p5_1->SetMaximum( 1.2 * fmax(histChargedIso_DR0p4To0p5_1->GetMaximum(),histChargedIso_DR0p4To0p5_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("ChargedIso_DR0p4To0p5"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon Gamma Isolation DR 0.0 to 0.1
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histGammaIso_DR0p0To0p1_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histGammaIso_DR0p0To0p1_2, legendLabel2.c_str(), "L");
  histGammaIso_DR0p0To0p1_1->SetLineColor(kBlue);
  histGammaIso_DR0p0To0p1_2->SetLineColor(kRed);
  histGammaIso_DR0p0To0p1_1->Draw("hist");
  histGammaIso_DR0p0To0p1_2->Draw("hist, same");
  histGammaIso_DR0p0To0p1_1->SetMaximum( 1.2 * fmax(histGammaIso_DR0p0To0p1_1->GetMaximum(),histGammaIso_DR0p0To0p1_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("GammaIso_DR0p0To0p1"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon Gamma Isolation DR 0.1 to 0.2
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histGammaIso_DR0p1To0p2_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histGammaIso_DR0p1To0p2_2, legendLabel2.c_str(), "L");
  histGammaIso_DR0p1To0p2_1->SetLineColor(kBlue);
  histGammaIso_DR0p1To0p2_2->SetLineColor(kRed);
  histGammaIso_DR0p1To0p2_1->Draw("hist");
  histGammaIso_DR0p1To0p2_2->Draw("hist, same");
  histGammaIso_DR0p1To0p2_1->SetMaximum( 1.2 * fmax(histGammaIso_DR0p1To0p2_1->GetMaximum(),histGammaIso_DR0p1To0p2_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("GammaIso_DR0p1To0p2"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon Gamma Isolation DR 0.2 To 0.3
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histGammaIso_DR0p2To0p3_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histGammaIso_DR0p2To0p3_2, legendLabel2.c_str(), "L");
  histGammaIso_DR0p2To0p3_1->SetLineColor(kBlue);
  histGammaIso_DR0p2To0p3_2->SetLineColor(kRed);
  histGammaIso_DR0p2To0p3_1->Draw("hist");
  histGammaIso_DR0p2To0p3_2->Draw("hist, same");
  histGammaIso_DR0p2To0p3_1->SetMaximum( 1.2 * fmax(histGammaIso_DR0p2To0p3_1->GetMaximum(),histGammaIso_DR0p2To0p3_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("GammaIso_DR0p2To0p3"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon Gamma Isolation DR 0.3 To 0.4
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histGammaIso_DR0p3To0p4_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histGammaIso_DR0p3To0p4_2, legendLabel2.c_str(), "L");
  histGammaIso_DR0p3To0p4_1->SetLineColor(kBlue);
  histGammaIso_DR0p3To0p4_2->SetLineColor(kRed);
  histGammaIso_DR0p3To0p4_1->Draw("hist");
  histGammaIso_DR0p3To0p4_2->Draw("hist, same");
  histGammaIso_DR0p3To0p4_1->SetMaximum( 1.2 * fmax(histGammaIso_DR0p3To0p4_1->GetMaximum(),histGammaIso_DR0p3To0p4_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("GammaIso_DR0p3To0p4"+Label+".gif").c_str() );

 //*******************************************************************************************
  //Photon Gamma Isolation 0.4 To 0.5
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histGammaIso_DR0p4To0p5_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histGammaIso_DR0p4To0p5_2, legendLabel2.c_str(), "L");
  histGammaIso_DR0p4To0p5_1->SetLineColor(kBlue);
  histGammaIso_DR0p4To0p5_2->SetLineColor(kRed);
  histGammaIso_DR0p4To0p5_1->Draw("hist");
  histGammaIso_DR0p4To0p5_2->Draw("hist, same");
  histGammaIso_DR0p4To0p5_1->SetMaximum( 1.2 * fmax(histGammaIso_DR0p4To0p5_1->GetMaximum(),histGammaIso_DR0p4To0p5_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("GammaIso_DR0p4To0p5"+Label+".gif").c_str() );

 //*******************************************************************************************
  //Photon Neutral Hadron Isolation DR 0.0 to 0.1
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histNeutralHadronIso_DR0p0To0p1_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histNeutralHadronIso_DR0p0To0p1_2, legendLabel2.c_str(), "L");
  histNeutralHadronIso_DR0p0To0p1_1->SetLineColor(kBlue);
  histNeutralHadronIso_DR0p0To0p1_2->SetLineColor(kRed);
  histNeutralHadronIso_DR0p0To0p1_1->Draw("hist");
  histNeutralHadronIso_DR0p0To0p1_2->Draw("hist, same");
  histNeutralHadronIso_DR0p0To0p1_1->SetMaximum( 1.2 * fmax(histNeutralHadronIso_DR0p0To0p1_1->GetMaximum(),histNeutralHadronIso_DR0p0To0p1_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("NeutralHadronIso_DR0p0To0p1"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon Neutral Hadron Isolation DR 0.1 to 0.2
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histNeutralHadronIso_DR0p1To0p2_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histNeutralHadronIso_DR0p1To0p2_2, legendLabel2.c_str(), "L");
  histNeutralHadronIso_DR0p1To0p2_1->SetLineColor(kBlue);
  histNeutralHadronIso_DR0p1To0p2_2->SetLineColor(kRed);
  histNeutralHadronIso_DR0p1To0p2_1->Draw("hist");
  histNeutralHadronIso_DR0p1To0p2_2->Draw("hist, same");
  histNeutralHadronIso_DR0p1To0p2_1->SetMaximum( 1.2 * fmax(histNeutralHadronIso_DR0p1To0p2_1->GetMaximum(),histNeutralHadronIso_DR0p1To0p2_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("NeutralHadronIso_DR0p1To0p2"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon Neutral Hadron Isolation DR 0.2 To 0.3
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histNeutralHadronIso_DR0p2To0p3_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histNeutralHadronIso_DR0p2To0p3_2, legendLabel2.c_str(), "L");
  histNeutralHadronIso_DR0p2To0p3_1->SetLineColor(kBlue);
  histNeutralHadronIso_DR0p2To0p3_2->SetLineColor(kRed);
  histNeutralHadronIso_DR0p2To0p3_1->Draw("hist");
  histNeutralHadronIso_DR0p2To0p3_2->Draw("hist, same");
  histNeutralHadronIso_DR0p2To0p3_1->SetMaximum( 1.2 * fmax(histNeutralHadronIso_DR0p2To0p3_1->GetMaximum(),histNeutralHadronIso_DR0p2To0p3_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("NeutralHadronIso_DR0p2To0p3"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon Neutral Hadron Isolation DR 0.3 To 0.4
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histNeutralHadronIso_DR0p3To0p4_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histNeutralHadronIso_DR0p3To0p4_2, legendLabel2.c_str(), "L");
  histNeutralHadronIso_DR0p3To0p4_1->SetLineColor(kBlue);
  histNeutralHadronIso_DR0p3To0p4_2->SetLineColor(kRed);
  histNeutralHadronIso_DR0p3To0p4_1->Draw("hist");
  histNeutralHadronIso_DR0p3To0p4_2->Draw("hist, same");
  histNeutralHadronIso_DR0p3To0p4_1->SetMaximum( 1.2 * fmax(histNeutralHadronIso_DR0p3To0p4_1->GetMaximum(),histNeutralHadronIso_DR0p3To0p4_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("NeutralHadronIso_DR0p3To0p4"+Label+".gif").c_str() );

 //*******************************************************************************************
  //Photon Neutral Hadron Isolation 0.4 To 0.5
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histNeutralHadronIso_DR0p4To0p5_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histNeutralHadronIso_DR0p4To0p5_2, legendLabel2.c_str(), "L");
  histNeutralHadronIso_DR0p4To0p5_1->SetLineColor(kBlue);
  histNeutralHadronIso_DR0p4To0p5_2->SetLineColor(kRed);
  histNeutralHadronIso_DR0p4To0p5_1->Draw("hist");
  histNeutralHadronIso_DR0p4To0p5_2->Draw("hist, same");
  histNeutralHadronIso_DR0p4To0p5_1->SetMaximum( 1.2 * fmax(histNeutralHadronIso_DR0p4To0p5_1->GetMaximum(),histNeutralHadronIso_DR0p4To0p5_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("NeutralHadronIso_DR0p4To0p5"+Label+".gif").c_str() );

 //*******************************************************************************************
  //Phton SC Raw Energy
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histSCRawEnergy_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histSCRawEnergy_2, legendLabel2.c_str(), "L");
  histSCRawEnergy_1->SetLineColor(kBlue);
  histSCRawEnergy_2->SetLineColor(kRed);
  histSCRawEnergy_1->Draw("hist");
  histSCRawEnergy_2->Draw("hist, same");
  histSCRawEnergy_1->SetMaximum( 1.2 * fmax(histSCRawEnergy_1->GetMaximum(),histSCRawEnergy_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_SCRawEnergy"+Label+".gif").c_str() );

 //*******************************************************************************************
  //Photon Pre-shower Over Raw
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPreShowerOverRaw_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPreShowerOverRaw_2, legendLabel2.c_str(), "L");
  histPreShowerOverRaw_1->SetLineColor(kBlue);
  histPreShowerOverRaw_2->SetLineColor(kRed);
  histPreShowerOverRaw_1->Draw("hist");
  histPreShowerOverRaw_2->Draw("hist, same");
  histPreShowerOverRaw_1->SetMaximum( 1.2 * fmax(histPreShowerOverRaw_1->GetMaximum(),histPreShowerOverRaw_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_R9"+Label+".gif").c_str() );

 //*******************************************************************************************
  //Photon R9
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhoR9_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhoR9_2, legendLabel2.c_str(), "L");
  histPhoR9_1->SetLineColor(kBlue);
  histPhoR9_2->SetLineColor(kRed);
  histPhoR9_1->Draw("hist");
  histPhoR9_2->Draw("hist, same");
  histPhoR9_1->SetMaximum( 1.2 * fmax(histPhoR9_1->GetMaximum(),histPhoR9_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_R9"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon NClusters
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histNClusters_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histNClusters_2, legendLabel2.c_str(), "L");
  histNClusters_1->SetLineColor(kBlue);
  histNClusters_2->SetLineColor(kRed);
  histNClusters_1->Draw("hist");
  histNClusters_2->Draw("hist, same");
  histNClusters_1->SetMaximum( 1.2 * fmax(histNClusters_1->GetMaximum(),histNClusters_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_NClusters"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon Eta Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histEtaSeed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histEtaSeed_2, legendLabel2.c_str(), "L");
  histEtaSeed_1->SetLineColor(kBlue);
  histEtaSeed_2->SetLineColor(kRed);
  histEtaSeed_1->Draw("hist");
  histEtaSeed_2->Draw("hist, same");
  histEtaSeed_1->SetMaximum( 1.2 * fmax(histEtaSeed_1->GetMaximum(),histEtaSeed_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_EtaSeed"+Label+".gif").c_str() );

 //*******************************************************************************************
  //Photon Phi Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhiSeed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhiSeed_2, legendLabel2.c_str(), "L");
  histPhiSeed_1->SetLineColor(kBlue);
  histPhiSeed_2->SetLineColor(kRed);
  histPhiSeed_1->Draw("hist");
  histPhiSeed_2->Draw("hist, same");
  histPhiSeed_1->SetMaximum( 1.2 * fmax(histPhiSeed_1->GetMaximum(),histPhiSeed_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_PhiSeed"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon E Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histESeed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histESeed_2, legendLabel2.c_str(), "L");
  histESeed_1->SetLineColor(kBlue);
  histESeed_2->SetLineColor(kRed);
  histESeed_1->Draw("hist");
  histESeed_2->Draw("hist, same");
  histESeed_1->SetMaximum( 1.2 * fmax(histESeed_1->GetMaximum(),histESeed_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_ESeed"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon E3x3 Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histE3x3Seed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histE3x3Seed_2, legendLabel2.c_str(), "L");
  histE3x3Seed_1->SetLineColor(kBlue);
  histE3x3Seed_2->SetLineColor(kRed);
  histE3x3Seed_1->Draw("hist");
  histE3x3Seed_2->Draw("hist, same");
  histE3x3Seed_1->SetMaximum( 1.2 * fmax(histE3x3Seed_1->GetMaximum(),histE3x3Seed_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_E3x3Seed"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon E5x5 Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histE5x5Seed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histE5x5Seed_2, legendLabel2.c_str(), "L");
  histE5x5Seed_1->SetLineColor(kBlue);
  histE5x5Seed_2->SetLineColor(kRed);
  histE5x5Seed_1->Draw("hist");
  histE5x5Seed_2->Draw("hist, same");
  histE5x5Seed_1->SetMaximum( 1.2 * fmax(histE5x5Seed_1->GetMaximum(),histE5x5Seed_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_E5x5Seed"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon EMax Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histEMaxSeed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histEMaxSeed_2, legendLabel2.c_str(), "L");
  histEMaxSeed_1->SetLineColor(kBlue);
  histEMaxSeed_2->SetLineColor(kRed);
  histEMaxSeed_1->Draw("hist");
  histEMaxSeed_2->Draw("hist, same");
  histEMaxSeed_1->SetMaximum( 1.2 * fmax(histEMaxSeed_1->GetMaximum(),histEMaxSeed_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_EMaxSeed"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon E 2nd Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histE2ndSeed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histE2ndSeed_2, legendLabel2.c_str(), "L");
  histE2ndSeed_1->SetLineColor(kBlue);
  histE2ndSeed_2->SetLineColor(kRed);
  histE2ndSeed_1->Draw("hist");
  histE2ndSeed_2->Draw("hist, same");
  histE2ndSeed_1->SetMaximum( 1.2 * fmax(histE2ndSeed_1->GetMaximum(),histE2ndSeed_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_E2ndSeed"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon E Top Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histETopSeed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histETopSeed_2, legendLabel2.c_str(), "L");
  histETopSeed_1->SetLineColor(kBlue);
  histETopSeed_2->SetLineColor(kRed);
  histETopSeed_1->Draw("hist");
  histETopSeed_2->Draw("hist, same");
  histETopSeed_1->SetMaximum( 1.2 * fmax(histETopSeed_1->GetMaximum(),histETopSeed_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_ETopSeed"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon E Bottom Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histEBottomSeed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histEBottomSeed_2, legendLabel2.c_str(), "L");
  histEBottomSeed_1->SetLineColor(kBlue);
  histEBottomSeed_2->SetLineColor(kRed);
  histEBottomSeed_1->Draw("hist");
  histEBottomSeed_2->Draw("hist, same");
  histEBottomSeed_1->SetMaximum( 1.2 * fmax(histEBottomSeed_1->GetMaximum(),histEBottomSeed_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_EBottomSeed"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon E Left Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histELeftSeed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histELeftSeed_2, legendLabel2.c_str(), "L");
  histELeftSeed_1->SetLineColor(kBlue);
  histELeftSeed_2->SetLineColor(kRed);
  histELeftSeed_1->Draw("hist");
  histELeftSeed_2->Draw("hist, same");
  histELeftSeed_1->SetMaximum( 1.2 * fmax(histELeftSeed_1->GetMaximum(),histELeftSeed_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_ELeftSeed"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon E Right Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histERightSeed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histERightSeed_2, legendLabel2.c_str(), "L");
  histERightSeed_1->SetLineColor(kBlue);
  histERightSeed_2->SetLineColor(kRed);
  histERightSeed_1->Draw("hist");
  histERightSeed_2->Draw("hist, same");
  histERightSeed_1->SetMaximum( 1.2 * fmax(histERightSeed_1->GetMaximum(),histERightSeed_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_ERightSeed"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon E2x5Max Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histE2x5MaxSeed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histE2x5MaxSeed_2, legendLabel2.c_str(), "L");
  histE2x5MaxSeed_1->SetLineColor(kBlue);
  histE2x5MaxSeed_2->SetLineColor(kRed);
  histE2x5MaxSeed_1->Draw("hist");
  histE2x5MaxSeed_2->Draw("hist, same");
  histE2x5MaxSeed_1->SetMaximum( 1.2 * fmax(histE2x5MaxSeed_1->GetMaximum(),histE2x5MaxSeed_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_E2x5MaxSeed"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon E 2x5Top Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histE2x5TopSeed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histE2x5TopSeed_2, legendLabel2.c_str(), "L");
  histE2x5TopSeed_1->SetLineColor(kBlue);
  histE2x5TopSeed_2->SetLineColor(kRed);
  histE2x5TopSeed_1->Draw("hist");
  histE2x5TopSeed_2->Draw("hist, same");
  histE2x5TopSeed_1->SetMaximum( 1.2 * fmax(histE2x5TopSeed_1->GetMaximum(),histE2x5TopSeed_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_E2x5TopSeed"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon E 2x5Bottom Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histE2x5BottomSeed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histE2x5BottomSeed_2, legendLabel2.c_str(), "L");
  histE2x5BottomSeed_1->SetLineColor(kBlue);
  histE2x5BottomSeed_2->SetLineColor(kRed);
  histE2x5BottomSeed_1->Draw("hist");
  histE2x5BottomSeed_2->Draw("hist, same");
  histE2x5BottomSeed_1->SetMaximum( 1.2 * fmax(histE2x5BottomSeed_1->GetMaximum(),histE2x5BottomSeed_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_E2x5BottomSeed"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon E 2x5Left Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histE2x5LeftSeed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histE2x5LeftSeed_2, legendLabel2.c_str(), "L");
  histE2x5LeftSeed_1->SetLineColor(kBlue);
  histE2x5LeftSeed_2->SetLineColor(kRed);
  histE2x5LeftSeed_1->Draw("hist");
  histE2x5LeftSeed_2->Draw("hist, same");
  histE2x5LeftSeed_1->SetMaximum( 1.2 * fmax(histE2x5LeftSeed_1->GetMaximum(),histE2x5LeftSeed_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_E2x5LeftSeed"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon E 2x5Right Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histE2x5RightSeed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histE2x5RightSeed_2, legendLabel2.c_str(), "L");
  histE2x5RightSeed_1->SetLineColor(kBlue);
  histE2x5RightSeed_2->SetLineColor(kRed);
  histE2x5RightSeed_1->Draw("hist");
  histE2x5RightSeed_2->Draw("hist, same");
  histE2x5RightSeed_1->SetMaximum( 1.2 * fmax(histE2x5RightSeed_1->GetMaximum(),histE2x5RightSeed_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_E2x5RightSeed"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon IEta Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histIEtaSeed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histIEtaSeed_2, legendLabel2.c_str(), "L");
  histIEtaSeed_1->SetLineColor(kBlue);
  histIEtaSeed_2->SetLineColor(kRed);
  histIEtaSeed_1->Draw("hist");
  histIEtaSeed_2->Draw("hist, same");
  histIEtaSeed_1->SetMaximum( 1.2 * fmax(histIEtaSeed_1->GetMaximum(),histIEtaSeed_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_IEtaSeed"+Label+".gif").c_str() );

 //*******************************************************************************************
  //Photon IPhi Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histIPhiSeed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histIPhiSeed_2, legendLabel2.c_str(), "L");
  histIPhiSeed_1->SetLineColor(kBlue);
  histIPhiSeed_2->SetLineColor(kRed);
  histIPhiSeed_1->Draw("hist");
  histIPhiSeed_2->Draw("hist, same");
  histIPhiSeed_1->SetMaximum( 1.2 * fmax(histIPhiSeed_1->GetMaximum(),histIPhiSeed_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_IPhiSeed"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon EtaCry Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histEtaCrySeed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histEtaCrySeed_2, legendLabel2.c_str(), "L");
  histEtaCrySeed_1->SetLineColor(kBlue);
  histEtaCrySeed_2->SetLineColor(kRed);
  histEtaCrySeed_1->Draw("hist");
  histEtaCrySeed_2->Draw("hist, same");
  histEtaCrySeed_1->SetMaximum( 1.2 * fmax(histEtaCrySeed_1->GetMaximum(),histEtaCrySeed_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_EtaCrySeed"+Label+".gif").c_str() );

 //*******************************************************************************************
  //Photon PhiCry Seed
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhiCrySeed_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhiCrySeed_2, legendLabel2.c_str(), "L");
  histPhiCrySeed_1->SetLineColor(kBlue);
  histPhiCrySeed_2->SetLineColor(kRed);
  histPhiCrySeed_1->Draw("hist");
  histPhiCrySeed_2->Draw("hist, same");
  histPhiCrySeed_1->SetMaximum( 1.2 * fmax(histPhiCrySeed_1->GetMaximum(),histPhiCrySeed_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_PhiCrySeed"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon Generated Energy
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histGeneratedEnergy_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histGeneratedEnergy_2, legendLabel2.c_str(), "L");
  histGeneratedEnergy_1->SetLineColor(kBlue);
  histGeneratedEnergy_2->SetLineColor(kRed);
  histGeneratedEnergy_1->Draw("hist");
  histGeneratedEnergy_2->Draw("hist, same");
  histGeneratedEnergy_1->SetMaximum( 1.2 * fmax(histGeneratedEnergy_1->GetMaximum(),histGeneratedEnergy_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_GeneratedEnergy"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon Momentum Regression V0
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhoMomentum_Regression_V0_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhoMomentum_Regression_V0_2, legendLabel2.c_str(), "L");
  histPhoMomentum_Regression_V0_1->SetLineColor(kBlue);
  histPhoMomentum_Regression_V0_2->SetLineColor(kRed);
  histPhoMomentum_Regression_V0_1->Draw("hist");
  histPhoMomentum_Regression_V0_2->Draw("hist, same");
  histPhoMomentum_Regression_V0_1->SetMaximum( 1.2 * fmax(histPhoMomentum_Regression_V0_1->GetMaximum(),histPhoMomentum_Regression_V0_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_PhoMomentum_Regression_V0"+Label+".gif").c_str() );

  //*******************************************************************************************
  //Photon Momentum Error Regression V0
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhoMomentumError_Regression_V0_1, legendLabel1.c_str(), "L");
  legend->AddEntry(histPhoMomentumError_Regression_V0_2, legendLabel2.c_str(), "L");
  histPhoMomentumError_Regression_V0_1->SetLineColor(kBlue);
  histPhoMomentumError_Regression_V0_2->SetLineColor(kRed);
  histPhoMomentumError_Regression_V0_1->Draw("hist");
  histPhoMomentumError_Regression_V0_2->Draw("hist, same");
  histPhoMomentumError_Regression_V0_1->SetMaximum( 1.2 * fmax(histPhoMomentumError_Regression_V0_1->GetMaximum(),histPhoMomentumError_Regression_V0_2->GetMaximum()));
  legend->Draw();
  cv->SaveAs(("PhotonDistribution_PhoMomentumError_Regression_V0"+Label+".gif").c_str() );
 
}

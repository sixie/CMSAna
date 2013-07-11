#ifndef CMSANA_UTILS_LEPTONTOOLS_HH
#define CMSANA_UTILS_LEPTONTOOLS_HH

#include <TMath.h>
#include <TClonesArray.h>          
#include "CMSAna/DataTree/interface/TElectron.hh"
#include "CMSAna/DataTree/interface/TMuon.hh"
#include "CMSAna/DataTree/interface/TPFCandidate.hh"
#include "CMSAna/ObjectStudies/interface/ElectronTree.h"
#include "CMSAna/ObjectStudies/interface/MuonTree.h"
#include "CMSAna/Utils/LeptonID.hh"
// #include "CMSAna/DataTree/interface/DataTreeDefs.hh"



void FillElectronTree(cmsana::ElectronTree *eleTree,
                      const cmsana::TElectron *ele, 
                      Int_t eleIndex,
                      TClonesArray *pfCandidates, 
                      Double_t rho, UInt_t DataEra,
                      UInt_t NPV,
                      UInt_t runNum,
                      UInt_t lumiSec,
                      UInt_t evtNum,
                      Float_t weight = 1.0
                      
  ) {

  eleTree->fWeight = weight;
  eleTree->fRunNumber = runNum;
  eleTree->fLumiSectionNumber = lumiSec;
  eleTree->fEventNumber = evtNum;
  eleTree->fEleEventNumberParity = (evtNum % 2 == 0);
  eleTree->fElePt = ele->pt; 
  eleTree->fEleEta = ele->eta; 
  eleTree->fElePhi = ele->phi; 
  eleTree->fEleSCEt = ele->scEt; 
  eleTree->fEleSCEta = ele->scEta; 
  eleTree->fEleSCPhi = ele->scPhi; 
  eleTree->fEleEcalEnergy = ele->EcalEnergy; 
  eleTree->fEleIsEcalDriven = ele->isEcalDriven;
  eleTree->fEleTriggerBit = 0;
  eleTree->fRho = rho; 
  eleTree->fNVertices = NPV; 
  eleTree->fEleD0 = ele->d0; 
  eleTree->fEleDZ = ele->dz; 
  eleTree->fEleIP3d = ele->ip3d; 
  eleTree->fEleIP3dSig = ele->ip3dSig; 
  eleTree->fEleMatchedConversion = ele->isConv;
  eleTree->fEleConvDCot = ele->partnerDeltaCot;
  eleTree->fEleConvDist = ele->partnerDist;
  eleTree->fEleNMissHits = ele->nExpHitsInner;
  eleTree->fEleNBrem = ele->nBrem; 
  eleTree->fEleFBrem = ele->fBrem; 
  eleTree->fEleEOverP = ele->EOverP; 
  eleTree->fEleESeedClusterOverPIn = ele->ESeedClusterOverPIn; 
  eleTree->fEleESeedClusterOverPout = ele->ESeedClusterOverPout; 
  eleTree->fEleEEleClusterOverPout = ele->EEleClusterOverPout; 
  eleTree->fEleOneOverEMinusOneOverP = (1.0/(ele->scEt*cosh(ele->scEta)) - 1.0 / ele->pIn);
  eleTree->fEleDEtaIn = fabs(ele->deltaEtaIn);
  eleTree->fEleDPhiIn = ele->deltaPhiIn; 
  eleTree->fEledEtaCalo = ele->dEtaCalo;
  eleTree->fEledPhiCalo = ele->dPhiCalo;
  eleTree->fEleSigmaIEtaIEta = ele->sigiEtaiEta; 
  eleTree->fEleSigmaIPhiIPhi = ele->sigiPhiiPhi;
  eleTree->fEleSigmaIEtaIPhi = ele->sigiEtaiPhi;
  eleTree->fEleSCEtaWidth = ele->SCEtaWidth;
  eleTree->fEleSCPhiWidth = ele->SCPhiWidth;
  eleTree->fEleR9 = ele->R9;
  eleTree->fElePreShowerOverRaw = ele->PreShowerOverRaw;
  eleTree->fEleHoverE = ele->HoverE; 
  eleTree->fEleGsfTrackChi2OverNdof = ele->GsfTrackChi2OverNdof;
  eleTree->fEleKFTrackChi2OverNDoF = ele->KFTrackChi2OverNdof;
  eleTree->fEleKFTrackNHits = ele->KFTrackNHits;
  eleTree->fEleKFTrackNLayersWithMeasurement = ele->KFTrackNLayersWithMeasurement;
  eleTree->fEleOneMinusSeedE1x5OverE5x5 = 1.0 - ele->SeedE1x5OverE5x5;
//   eleTree->fElePFMVA = ele->mva;
  eleTree->fEleTrkIso03 = ele->trkIso03; 
  eleTree->fEleEMIso03 = ele->emIso03; 
  eleTree->fEleHadIso03 = ele->hadIso03; 
  eleTree->fEleTrkIso04 = ele->trkIso04; 
  eleTree->fEleEMIso04 = ele->emIso04; 
  eleTree->fEleHadIso04 = ele->hadIso04; 
  eleTree->fElePFIso04 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFIso, 0.0, 0.4);
  eleTree->fChargedIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFChargedIso, 0.0, 0.1)/ele->pt;
  eleTree->fChargedIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFChargedIso, 0.1, 0.2)/ele->pt;
  eleTree->fChargedIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFChargedIso, 0.2, 0.3)/ele->pt;
  eleTree->fChargedIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFChargedIso, 0.3, 0.4)/ele->pt;
  eleTree->fChargedIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFChargedIso, 0.4, 0.5)/ele->pt;
  eleTree->fGammaIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFGammaIso, 0.0, 0.1)/ele->pt;
  eleTree->fGammaIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFGammaIso, 0.1, 0.2)/ele->pt;
  eleTree->fGammaIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFGammaIso, 0.2, 0.3)/ele->pt;
  eleTree->fGammaIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFGammaIso, 0.3, 0.4)/ele->pt;
  eleTree->fGammaIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFGammaIso, 0.4, 0.5)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p0To0p1 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.0, 0.1)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p1To0p2 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.1, 0.2)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p2To0p3 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.2, 0.3)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p3To0p4 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.3, 0.4)/ele->pt;
  eleTree->fNeutralHadronIso_DR0p4To0p5 = ComputeElePFIsoRings(ele, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.4, 0.5)/ele->pt;
//   eleTree->fElePassTriggerDenominator = passElectronTriggerDenominator(ele);

  //**********************************************
  //Fill Variables for Regression
  //**********************************************
  eleTree->fIsEB = ele->isEB;
  eleTree->fIsEE = ele->isEE; 
  
  eleTree->fSCRawEnergy = ele->SCRawEnergy ;
  eleTree->fNClusters = ele->nBrem+1 ;
  eleTree->fEtaSeed = ele->EtaSeed ;
  eleTree->fPhiSeed = ele->PhiSeed ;
  eleTree->fESeed = ele->ESeed ;
  eleTree->fE3x3Seed = ele->E3x3Seed ;
  eleTree->fE5x5Seed = ele->E5x5Seed ;
  eleTree->fEMaxSeed = ele->EMaxSeed ;
  eleTree->fE2ndSeed = ele->E2ndSeed ;
  eleTree->fETopSeed = ele->ETopSeed ;
  eleTree->fEBottomSeed = ele->EBottomSeed ;
  eleTree->fELeftSeed = ele->ELeftSeed ;
  eleTree->fERightSeed = ele->ERightSeed ;
  eleTree->fE2x5MaxSeed = ele->E2x5MaxSeed ;
  eleTree->fE2x5TopSeed = ele->E2x5TopSeed ;
  eleTree->fE2x5BottomSeed = ele->E2x5BottomSeed ;
  eleTree->fE2x5LeftSeed = ele->E2x5LeftSeed ;
  eleTree->fE2x5RightSeed = ele->E2x5RightSeed ;
  eleTree->fIEtaSeed = ele->IEtaSeed ;
  eleTree->fIPhiSeed = ele->IPhiSeed ;
  eleTree->fEtaCrySeed = ele->EtaCrySeed ;
  eleTree->fPhiCrySeed = ele->PhiCrySeed ;
  eleTree->fEcalEnergyError = ele->EcalEnergyError;
  eleTree->fGsfTrackPIn = ele->pIn ;
  eleTree->fCharge = ele->q ;
  eleTree->fTrackMomentumError = fmin(ele->TrackMomentumError,500.0);
  eleTree->fEleClassification = ele->Classification;

  //**********************************************
  //Make final consistency checks before filling
  //**********************************************
  if (TMath::IsNaN(ele->sigiPhiiPhi)) {
    cout << "Problem with sigiPhiiPhi: NaN. exit.\n";
    assert(0);
  }

  //***********************
  //Fill Electron
  //***********************
  eleTree->tree_->Fill();
}








// void FillMuonTree(cmsana::MuonTree *muTree,
//                       const cmsana::TMuon *mu, 
//                       TClonesArray *pfCandidates, 
//                       Double_t rho, UInt_t DataEra,
//                       UInt_t NPV,
//                       UInt_t runNum,
//                       UInt_t lumiSec,
//                       UInt_t evtNum,
//                       Float_t weight = 1.0
                      
//   ) {


//   //Fill These Muons

//   muTree->fWeight = weight;
//   muTree->fRunNumber = runNum;
//   muTree->fLumiSectionNumber = lumiSec;
//   muTree->fEventNumber = evtNum;
//   muTree->fMuEventNumberParity = (evtNum % 2 == 0);
//   muTree->fRho = rho; 
//   muTree->fNVertices = NPV; 

//   muTree->fMuPt = mu->pt; 
//   muTree->fMuEta = mu->eta; 
//   muTree->fMuPhi = mu->phi; 

//   muTree->fMuTypeBits = mu->typeBits;
//   muTree->fIsAllArbitrated = ((mu->qualityBits & kAllArbitrated) == kAllArbitrated);
//   muTree->fMuTkNchi2 = mu->tkNchi2; 
//   muTree->fMuGlobalNchi2 = mu->muNchi2; 
//   muTree->fMuNValidHits = mu->nValidHits; 
//   muTree->fMuNTrackerHits = mu->nTkHits; 
//   muTree->fMuNPixelHits = mu->nPixHits; 
//   muTree->fMuNMatches = mu->nMatch ; 
//   muTree->fMuD0 = mu->d0; 

//   //Additional Vars 
//   muTree->fMuIP3d = mu->ip3d ; 
//   muTree->fMuIP3dSig = mu->ip3dSig ; 
//   muTree->fMuTrkKink = mu->TrkKink ; 
//   muTree->fMuGlobalKink = mu->GlobalKink ; 
//   muTree->fMuSegmentCompatibility = mu->SegmentCompatibility ; 
//   muTree->fMuCaloCompatibility = mu->CaloCompatilibity ; 
//   muTree->fMuHadEnergy = mu->HadEnergy; 
//   muTree->fMuHoEnergy = mu->HoEnergy; 
//   muTree->fMuEmEnergy = mu->EmEnergy; 
//   muTree->fMuHadS9Energy = mu->HadS9Energy; 
//   muTree->fMuHoS9Energy = mu->HoS9Energy; 
//   muTree->fMuEmS9Energy = mu->EmS9Energy; 

//   //Isolation Variables
//   muTree->fMuTrkIso03 = mu->trkIso03; 
//   muTree->fMuEMIso03 = mu->emIso03; 
//   muTree->fMuHadIso03 = mu->hadIso03; 
//   muTree->fMuTrkIso05 = mu->trkIso05; 
//   muTree->fMuEMIso05 = mu->emIso05; 
//   muTree->fMuHadIso05 = mu->hadIso05; 
//   muTree->fMuPFIso04 = ComputeMuonPFIso04( mu, pfCandidates, rho, DataEra); 

//   muTree->fChargedIso_DR0p0To0p1 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFChargedIso, 0.0, 0.1);
//   muTree->fChargedIso_DR0p1To0p2 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFChargedIso, 0.1, 0.2);
//   muTree->fChargedIso_DR0p2To0p3 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFChargedIso, 0.2, 0.3);
//   muTree->fChargedIso_DR0p3To0p4 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFChargedIso, 0.3, 0.4);
//   muTree->fChargedIso_DR0p4To0p5 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFChargedIso, 0.4, 0.5);
//   muTree->fGammaIso_DR0p0To0p1 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFGammaIso, 0.0, 0.1);
//   muTree->fGammaIso_DR0p1To0p2 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFGammaIso, 0.1, 0.2);
//   muTree->fGammaIso_DR0p2To0p3 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFGammaIso, 0.2, 0.3);
//   muTree->fGammaIso_DR0p3To0p4 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFGammaIso, 0.3, 0.4);
//   muTree->fGammaIso_DR0p4To0p5 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFGammaIso, 0.4, 0.5);
//   muTree->fNeutralHadronIso_DR0p0To0p1 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.0, 0.1);
//   muTree->fNeutralHadronIso_DR0p1To0p2 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.1, 0.2);
//   muTree->fNeutralHadronIso_DR0p2To0p3 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.2, 0.3);
//   muTree->fNeutralHadronIso_DR0p3To0p4 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.3, 0.4);
//   muTree->fNeutralHadronIso_DR0p4To0p5 = ComputeMuonPFIsoRings(mu, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.4, 0.5);
//   muTree->fMuPassTriggerDenominator = passMuonTriggerDenominator(mu, pfCandidates, rho, DataEra);

//   //***********************
//   //Fill Muon
//   //***********************
//   muTree->tree_->Fill();



// }



// void FillMuonTree(cmsana::MuonTree *muTree,
//                   const cmsana::TMuon *mu, 
//                   TClonesArray *pfCandidates, 
//                   vector<const cmsana::TPFCandidate*> particlesToVeto,
//                   Double_t rho, UInt_t DataEra,
//                   UInt_t NPV,
//                   UInt_t runNum,
//                   UInt_t lumiSec,
//                   UInt_t evtNum,
//                   Float_t weight = 1.0
                  
//   ) {


//   //Fill These Muons

//   muTree->fWeight = weight;
//   muTree->fRunNumber = runNum;
//   muTree->fLumiSectionNumber = lumiSec;
//   muTree->fEventNumber = evtNum;
//   muTree->fMuEventNumberParity = (evtNum % 2 == 0);
//   muTree->fRho = rho; 
//   muTree->fNVertices = NPV; 

//   muTree->fMuPt = mu->pt; 
//   muTree->fMuEta = mu->eta; 
//   muTree->fMuPhi = mu->phi; 

//   muTree->fMuTypeBits = mu->typeBits;
//   muTree->fIsAllArbitrated = ((mu->qualityBits & kAllArbitrated) == kAllArbitrated);
//   muTree->fMuTkNchi2 = mu->tkNchi2; 
//   muTree->fMuGlobalNchi2 = mu->muNchi2; 
//   muTree->fMuNValidHits = mu->nValidHits; 
//   muTree->fMuNTrackerHits = mu->nTkHits; 
//   muTree->fMuNPixelHits = mu->nPixHits; 
//   muTree->fMuNMatches = mu->nMatch ; 
//   muTree->fMuD0 = mu->d0; 

//   //Additional Vars 
//   muTree->fMuIP3d = mu->ip3d ; 
//   muTree->fMuIP3dSig = mu->ip3dSig ; 
//   muTree->fMuTrkKink = mu->TrkKink ; 
//   muTree->fMuGlobalKink = mu->GlobalKink ; 
//   muTree->fMuSegmentCompatibility = mu->SegmentCompatibility ; 
//   muTree->fMuCaloCompatibility = mu->CaloCompatilibity ; 
//   muTree->fMuHadEnergy = mu->HadEnergy; 
//   muTree->fMuHoEnergy = mu->HoEnergy; 
//   muTree->fMuEmEnergy = mu->EmEnergy; 
//   muTree->fMuHadS9Energy = mu->HadS9Energy; 
//   muTree->fMuHoS9Energy = mu->HoS9Energy; 
//   muTree->fMuEmS9Energy = mu->EmS9Energy; 

//   //Isolation Variables
//   muTree->fMuTrkIso03 = mu->trkIso03; 
//   muTree->fMuEMIso03 = mu->emIso03; 
//   muTree->fMuHadIso03 = mu->hadIso03; 
//   muTree->fMuTrkIso05 = mu->trkIso05; 
//   muTree->fMuEMIso05 = mu->emIso05; 
//   muTree->fMuHadIso05 = mu->hadIso05; 
//   muTree->fMuPFIso04 = ComputeMuonPFIso04( mu, pfCandidates, rho, DataEra); 

//   muTree->fChargedIso_DR0p0To0p1 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFChargedIso, 0.0, 0.1);
//   muTree->fChargedIso_DR0p1To0p2 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFChargedIso, 0.1, 0.2);
//   muTree->fChargedIso_DR0p2To0p3 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFChargedIso, 0.2, 0.3);
//   muTree->fChargedIso_DR0p3To0p4 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFChargedIso, 0.3, 0.4);
//   muTree->fChargedIso_DR0p4To0p5 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFChargedIso, 0.4, 0.5);
//   muTree->fGammaIso_DR0p0To0p1 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFGammaIso, 0.0, 0.1);
//   muTree->fGammaIso_DR0p1To0p2 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFGammaIso, 0.1, 0.2);
//   muTree->fGammaIso_DR0p2To0p3 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFGammaIso, 0.2, 0.3);
//   muTree->fGammaIso_DR0p3To0p4 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFGammaIso, 0.3, 0.4);
//   muTree->fGammaIso_DR0p4To0p5 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFGammaIso, 0.4, 0.5);
//   muTree->fNeutralHadronIso_DR0p0To0p1 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFNeutralHadronIso, 0.0, 0.1);
//   muTree->fNeutralHadronIso_DR0p1To0p2 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFNeutralHadronIso, 0.1, 0.2);
//   muTree->fNeutralHadronIso_DR0p2To0p3 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFNeutralHadronIso, 0.2, 0.3);
//   muTree->fNeutralHadronIso_DR0p3To0p4 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFNeutralHadronIso, 0.3, 0.4);
//   muTree->fNeutralHadronIso_DR0p4To0p5 = ComputeMuonPFIsoRings(mu, pfCandidates, particlesToVeto, rho, DataEra, kPFNeutralHadronIso, 0.4, 0.5);
//   muTree->fMuPassTriggerDenominator = passMuonTriggerDenominator(mu, pfCandidates, rho, DataEra);

//   //***********************
//   //Fill Muon
//   //***********************
//   muTree->tree_->Fill();



// }



#endif

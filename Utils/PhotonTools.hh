#ifndef CMSANA_UTILS_PHOTONTOOLS_HH
#define CMSANA_UTILS_PHOTONTOOLS_HH

#include <TMath.h>
#include <TClonesArray.h>          
#include "CMSAna/DataTree/interface/TPhoton.hh"
#include "CMSAna/DataTree/interface/TPFCandidate.hh"
#include "CMSAna/ObjectStudies/interface/PhotonTree.h"
#include "CMSAna/Utils/PhotonID.hh"
// #include "CMSAna/DataTree/interface/DataTreeDefs.hh"



void FillPhotonTree(cmsana::PhotonTree *phoTree,
                    const cmsana::TGenParticle *genpho, 
                    const cmsana::TPhoton *pho, 
                    TClonesArray *pfCandidates, 
                    TClonesArray *genParticles, 
                    Double_t rho, UInt_t DataEra,
                    cmsana::TEventInfo *info,
                    Float_t weight = 1.0
                    
  ) {

  if (genpho) {
    phoTree->fPhoGenPt = genpho->pt;
    phoTree->fPhoGenEta = genpho->eta;
    phoTree->fPhoGenPhi = genpho->phi;
    phoTree->fGeneratedEnergy = genpho->pt * cosh(genpho->eta);
  } else {
    phoTree->fPhoGenPt = -99.0;
    phoTree->fPhoGenEta = -99.0;
    phoTree->fPhoGenPhi = -99.0;
    phoTree->fGeneratedEnergy = -99.0;
  }

  if (pho) {
    phoTree->fWeight = weight;
    phoTree->fRunNumber = info->runNum;
    phoTree->fLumiSectionNumber = info->lumiSec;
    phoTree->fEventNumber = info->evtNum;
    phoTree->fPhoEventNumberParity = (info->evtNum % 2 == 0);
    phoTree->fPhoPt = pho->pt; 
    phoTree->fPhoEta = pho->eta; 
    phoTree->fPhoPhi = pho->phi; 
    phoTree->fPhoSCEt = pho->scEt; 
    phoTree->fPhoSCEta = pho->scEta; 
    phoTree->fPhoSCPhi = pho->scPhi; 
    phoTree->fPhoTriggerBit = 0;
    phoTree->fRho = rho; 
    phoTree->fNVertices = info->nGoodPV; 
    phoTree->fPhoHasPixelSeed = pho->hasPixelSeed; 
    phoTree->fPhoIsConversion = pho->isConversion;
    phoTree->fPhoPassEleVeto = pho->passEleVeto;
    phoTree->fPhoHoverE = pho->HoverE; 
    phoTree->fPhoHoverESingleTower = pho->HoverESingleTower; 
    phoTree->fPhoSigmaIEtaIEta = pho->sigiEtaiEta; 
    phoTree->fPhoSigmaIPhiIPhi = pho->sigiPhiiPhi;
    phoTree->fPhoSigmaIEtaIPhi = pho->sigiEtaiPhi;
    phoTree->fPhoSCEtaWidth = pho->SCEtaWidth;
    phoTree->fPhoSCPhiWidth = pho->SCPhiWidth;
    phoTree->fPhoR9 = pho->R9;
    phoTree->fPhoTrkIso03 = pho->trkIso03Hollow; 
    phoTree->fPhoEMIso03 = pho->emIso03; 
    phoTree->fPhoHadIso03 = pho->hadIso03; 
    phoTree->fPhoPFIso03 = ComputePhotonPFIsoRings(pho, pfCandidates, rho, DataEra, kPFIso, 0.0, 0.3);
    phoTree->fPhoPFIso04 = ComputePhotonPFIsoRings(pho, pfCandidates, rho, DataEra, kPFIso, 0.0, 0.4);
    phoTree->fChargedIso_DR0p0To0p1 = ComputePhotonPFIsoRings(pho, pfCandidates, rho, DataEra, kPFChargedIso, 0.0, 0.1);
    phoTree->fChargedIso_DR0p1To0p2 = ComputePhotonPFIsoRings(pho, pfCandidates, rho, DataEra, kPFChargedIso, 0.1, 0.2);
    phoTree->fChargedIso_DR0p2To0p3 = ComputePhotonPFIsoRings(pho, pfCandidates, rho, DataEra, kPFChargedIso, 0.2, 0.3);
    phoTree->fChargedIso_DR0p3To0p4 = ComputePhotonPFIsoRings(pho, pfCandidates, rho, DataEra, kPFChargedIso, 0.3, 0.4);
    phoTree->fChargedIso_DR0p4To0p5 = ComputePhotonPFIsoRings(pho, pfCandidates, rho, DataEra, kPFChargedIso, 0.4, 0.5);
    phoTree->fGammaIso_DR0p0To0p1 = ComputePhotonPFIsoRings(pho, pfCandidates, rho, DataEra, kPFGammaIso, 0.0, 0.1);
    phoTree->fGammaIso_DR0p1To0p2 = ComputePhotonPFIsoRings(pho, pfCandidates, rho, DataEra, kPFGammaIso, 0.1, 0.2);
    phoTree->fGammaIso_DR0p2To0p3 = ComputePhotonPFIsoRings(pho, pfCandidates, rho, DataEra, kPFGammaIso, 0.2, 0.3);
    phoTree->fGammaIso_DR0p3To0p4 = ComputePhotonPFIsoRings(pho, pfCandidates, rho, DataEra, kPFGammaIso, 0.3, 0.4);
    phoTree->fGammaIso_DR0p4To0p5 = ComputePhotonPFIsoRings(pho, pfCandidates, rho, DataEra, kPFGammaIso, 0.4, 0.5);
    phoTree->fNeutralHadronIso_DR0p0To0p1 = ComputePhotonPFIsoRings(pho, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.0, 0.1);
    phoTree->fNeutralHadronIso_DR0p1To0p2 = ComputePhotonPFIsoRings(pho, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.1, 0.2);
    phoTree->fNeutralHadronIso_DR0p2To0p3 = ComputePhotonPFIsoRings(pho, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.2, 0.3);
    phoTree->fNeutralHadronIso_DR0p3To0p4 = ComputePhotonPFIsoRings(pho, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.3, 0.4);
    phoTree->fNeutralHadronIso_DR0p4To0p5 = ComputePhotonPFIsoRings(pho, pfCandidates, rho, DataEra, kPFNeutralHadronIso, 0.4, 0.5);


    //**********************************************
    //Fill Variables for Regression
    //**********************************************
    phoTree->fSCRawEnergy = pho->SCRawEnergy ;
    phoTree->fPreShowerOverRaw = pho->PreShowerOverRaw;
    phoTree->fNClusters = pho->NClusters ;
    phoTree->fEtaSeed = pho->EtaSeed ;
    phoTree->fPhiSeed = pho->PhiSeed ;
    phoTree->fESeed = pho->ESeed ;
    phoTree->fE3x3Seed = pho->E3x3Seed ;
    phoTree->fE5x5Seed = pho->E5x5Seed ;
    phoTree->fEMaxSeed = pho->EMaxSeed ;
    phoTree->fE2ndSeed = pho->E2ndSeed ;
    phoTree->fETopSeed = pho->ETopSeed ;
    phoTree->fEBottomSeed = pho->EBottomSeed ;
    phoTree->fELeftSeed = pho->ELeftSeed ;
    phoTree->fERightSeed = pho->ERightSeed ;
    phoTree->fE2x5MaxSeed = pho->E2x5MaxSeed ;
    phoTree->fE2x5TopSeed = pho->E2x5TopSeed ;
    phoTree->fE2x5BottomSeed = pho->E2x5BottomSeed ;
    phoTree->fE2x5LeftSeed = pho->E2x5LeftSeed ;
    phoTree->fE2x5RightSeed = pho->E2x5RightSeed ;
    phoTree->fIEtaSeed = pho->IEtaSeed ;
    phoTree->fIPhiSeed = pho->IPhiSeed ;
    phoTree->fEtaCrySeed = pho->EtaCrySeed ;
    phoTree->fPhiCrySeed = pho->PhiCrySeed ;

    phoTree->fPhoMomentum_Regression_V0 = 0;
    phoTree->fPhoMomentumError_Regression_V0 = 0;
  } else {
    phoTree->fWeight = 0;
    phoTree->fRunNumber = 0;
    phoTree->fLumiSectionNumber = 0;
    phoTree->fEventNumber = 0;
    phoTree->fPhoEventNumberParity = 0;
    phoTree->fPhoPt = -99.0;
    phoTree->fPhoEta =  -99.0;
    phoTree->fPhoPhi =  -99.0;
    phoTree->fPhoSCEt =  -99.0;
    phoTree->fPhoSCEta =  -99.0;
    phoTree->fPhoSCPhi =  -99.0;
    phoTree->fPhoTriggerBit = 0;
    phoTree->fRho = 0;
    phoTree->fNVertices = 0;
    phoTree->fPhoHasPixelSeed = false;
    phoTree->fPhoIsConversion = false;
    phoTree->fPhoPassEleVeto = false;
    phoTree->fPhoHoverE =  -99.0;
    phoTree->fPhoHoverESingleTower =  -99.0;
    phoTree->fPhoSigmaIEtaIEta = -99.0;
    phoTree->fPhoSigmaIPhiIPhi =  -99.0;
    phoTree->fPhoSigmaIEtaIPhi =  -99.0;
    phoTree->fPhoSCEtaWidth =  -99.0;
    phoTree->fPhoSCPhiWidth =  -99.0;
    phoTree->fPhoR9 =  -99.0;
    phoTree->fPhoTrkIso03 =  -99.0;
    phoTree->fPhoEMIso03 =  -99.0;
    phoTree->fPhoHadIso03 =  -99.0;
    phoTree->fPhoPFIso03 =  -99.0;
    phoTree->fPhoPFIso04 =  -99.0;
    phoTree->fChargedIso_DR0p0To0p1 =  -99.0;
    phoTree->fChargedIso_DR0p1To0p2 =  -99.0;
    phoTree->fChargedIso_DR0p2To0p3 =  -99.0;
    phoTree->fChargedIso_DR0p3To0p4 =  -99.0;
    phoTree->fChargedIso_DR0p4To0p5 =  -99.0;
    phoTree->fGammaIso_DR0p0To0p1 =  -99.0;
    phoTree->fGammaIso_DR0p1To0p2 =  -99.0;
    phoTree->fGammaIso_DR0p2To0p3 =  -99.0;
    phoTree->fGammaIso_DR0p3To0p4 =  -99.0;
    phoTree->fGammaIso_DR0p4To0p5 =  -99.0;
    phoTree->fNeutralHadronIso_DR0p0To0p1 =  -99.0;
    phoTree->fNeutralHadronIso_DR0p1To0p2 =  -99.0;
    phoTree->fNeutralHadronIso_DR0p2To0p3 =  -99.0;
    phoTree->fNeutralHadronIso_DR0p3To0p4 =  -99.0;
    phoTree->fNeutralHadronIso_DR0p4To0p5 =  -99.0;
    phoTree->fSCRawEnergy =  -99.0;
    phoTree->fPreShowerOverRaw =  -99.0;
    phoTree->fNClusters =  0;
    phoTree->fEtaSeed = -99.0;
    phoTree->fPhiSeed = -99.0;
    phoTree->fESeed = -99.0;
    phoTree->fE3x3Seed = -99.0;
    phoTree->fE5x5Seed = -99.0;
    phoTree->fEMaxSeed = -99.0;
    phoTree->fE2ndSeed = -99.0;
    phoTree->fETopSeed = -99.0;
    phoTree->fEBottomSeed = -99.0;
    phoTree->fELeftSeed = -99.0;
    phoTree->fERightSeed = -99.0;
    phoTree->fE2x5MaxSeed = -99.0;
    phoTree->fE2x5TopSeed = -99.0;
    phoTree->fE2x5BottomSeed = -99.0;
    phoTree->fE2x5LeftSeed = -99.0;
    phoTree->fE2x5RightSeed = -99.0;
    phoTree->fIEtaSeed = -99.0;
    phoTree->fIPhiSeed = -99.0;
    phoTree->fEtaCrySeed = -99.0;
    phoTree->fPhiCrySeed = -99.0;
    phoTree->fPhoMomentum_Regression_V0 = -99.0;
    phoTree->fPhoMomentumError_Regression_V0 = -99.0;
  }

  //***********************
  //Fill Photon
  //***********************
  phoTree->tree_->Fill();
}



#endif

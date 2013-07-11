#ifndef PHOTONID_HH
#define PHOTONID_HH

#include "CMSAna/Utils/CommonDefs.hh"
#include "CMSAna/Utils/CommonTools.hh"
#include "CMSAna/Utils/IsolationPileupCorrections.hh"
#include "CMSAna/DataTree/interface/TPhoton.hh"
#include "CMSAna/DataTree/interface/DataTreeDefs.hh"
#include <cassert>
#include <iostream>
#include "TMath.h"

Double_t ComputePhotonPFIsoRings( const cmsana::TPhoton *photon,
                                  TClonesArray *pfCandidates, 
                                  Double_t rho, UInt_t DataEra,
                                  UInt_t pfType,
                                  Double_t minDR,
                                  Double_t maxDR,
                                  Bool_t printDebug = kFALSE);

Bool_t passPhotonIDSimpleMedium ( const cmsana::TPhoton *pho,
                                  TClonesArray *pfCandidates,  
                                  double rho, UInt_t DataEra,
                                  Bool_t printDebug = false );



//=== FUNCTION DEFINITIONS ======================================================================================





//--------------------------------------------------------------------------------------------------
Double_t ComputePhotonPFIsoRings( const cmsana::TPhoton *photon,
                                TClonesArray *pfCandidates, 
                                Double_t rho, UInt_t DataEra,
                                UInt_t pfIsoType,
                                Double_t minDR,
                                Double_t maxDR,
                                Bool_t printDebug) {

  Double_t fIso = 0.0;

  //
  //Loop over PF Candidates
  //
  for(UInt_t k=0; int(k) < pfCandidates->GetEntries(); ++k) {
    
    const cmsana::TPFCandidate *pf = (cmsana::TPFCandidate*)((*pfCandidates)[k]);

    //*****************************
    // veto Pileup PF Candidates
    //*****************************
    if (!pf->IsPFNoPU) continue;
    
    //*****************************
    //don't include pf electrons or pf muons
    //*****************************
    if (pf->pfType == eElectron || pf->pfType == eMuon ) continue;

    //************************************************
    // veto PF Candidates that match to the photon SC
    //************************************************
    bool SCMatch = false;
    for (UInt_t p=0; p < photon->NMatchedPFCandidates; ++p) {
      if (photon->MatchedPFCandidateIndex[p] == k) SCMatch = true;
    }
    if (SCMatch) continue;


    Double_t dr = cmsana::deltaR(photon->eta, photon->phi, pf->eta, pf->phi);
    if (dr < minDR) continue;
    if (dr >= maxDR) continue;

    if (printDebug) {
      cout << "PFCandidate type " << pf->pfType << " : " << pf->pt << " | " 
           << dr << " : " << fabs(pf->eta-photon->eta) << " | " << fabs(pf->dz) << " " ;
    }
    
    //***********************************************************************
    // Charged Iso 
    //***********************************************************************
    if (pf->q != 0) {      

      //inner DR veto
      if (dr > 0.02) {      
        if (pfIsoType == kPFChargedIso || pfIsoType == kPFIso) {
          fIso += pf->pt;
          if (printDebug) cout << " ADD";
        }
      }
    }
  
    //***********************************************************************
    // Gamma Iso
    //***********************************************************************
    else if (pf->pfType == eGamma) {

      if (fabs(photon->scEta) < 1.4442) {        
        //eta strip veto for barrel
        if( fabs(pf->eta-photon->eta) > 0.015 ) {
          if (pfIsoType == kPFGammaIso || pfIsoType == kPFIso) fIso += pf->pt;
          if (printDebug) cout << " ADD";
        }
      } else {        
        //DR veto for endcap
        if (dr > 0.00864*fabs(sinh(photon->scEta))*4) {
          if (pfIsoType == kPFGammaIso || pfIsoType == kPFIso) fIso += pf->pt;
          if (printDebug) cout << " ADD";
        }
      }

    }
    
    //***********************************************************************
    // Other Neutrals
    //***********************************************************************
    else {
      if( pf->pt > 0.5 ) {
	if (pfIsoType == kPFNeutralHadronIso || pfIsoType == kPFIso) fIso += pf->pt;
        if (printDebug) cout << " ADD";
      }
    }

    if (printDebug) cout << "\n";
  }



  if(printDebug) std::cout << "rho: " << rho << std::endl;

  UInt_t EAType = 0;
  if (pfIsoType == kPFChargedIso) {
    if (minDR == 0.0 && maxDR == 0.3) EAType = kPhotonChargedIso03;
    else EAType = kPhotonNoCorrection;
  } else if (pfIsoType == kPFGammaIso) {
    if (minDR == 0.0 && maxDR == 0.3) EAType = kPhotonGammaIso03;
    else EAType = kPhotonNoCorrection;
  } else if (pfIsoType == kPFNeutralHadronIso) {
    if (minDR == 0.0 && maxDR == 0.3) EAType = kPhotonNeutralHadronIso03;    
    else EAType = kPhotonNoCorrection;
  } else {
    EAType = kPhotonNoCorrection;
  }

  double iso = fmax(0.0, fIso - rho*PhotonEffectiveArea(EAType, photon->eta,DataEra));
  
  if( printDebug ) { 
    std::cout << "PFIso : type = " << pfIsoType << " minDR = " << minDR << " maxDR = " << maxDR << " | " 
              << "uncorrected iso: " << fIso
              << "\tcorrected iso: " << iso
              << std::endl;
  }

  return iso;
 
} 
  




Bool_t passPhotonIDSimpleLoose ( const cmsana::TPhoton *pho,
                                  TClonesArray *pfCandidates,  
                                  double rho, UInt_t DataEra,
                                  bool printDebug ) {
  
  Bool_t pass = true;
  double PFChargedHadronIso = ComputePhotonPFIsoRings(pho, pfCandidates, rho, DataEra, 
                                                      kPFChargedIso,
                                                      0.0,0.3,
                                                      printDebug);
  double PFPhotonIso = ComputePhotonPFIsoRings(pho, pfCandidates, rho, DataEra, 
                                               kPFGammaIso,
                                               0.0,0.3,
                                               printDebug);
  double PFNeutralHadronIso = ComputePhotonPFIsoRings(pho, pfCandidates, rho, DataEra, 
                                                      kPFNeutralHadronIso,
                                                      0.0,0.3,
                                                      printDebug);

  if (abs(pho->scEta) < 1.4442) {
      if (!( (0==0) 
             && pho->sigiEtaiEta < 0.012
             && pho->HoverESingleTower < 0.05
             && pho->passEleVeto
             && PFChargedHadronIso < 2.6
             && PFPhotonIso < 1.3 + 0.005 * pho->pt
             && PFNeutralHadronIso < 3.5 + 0.04 * pho->pt
            )) 
        pass = false;     
  } else {
    if (!( (0==0) 
           && pho->sigiEtaiEta < 0.034 
           && pho->HoverESingleTower < 0.05
           && pho->passEleVeto
           && PFChargedHadronIso < 2.3
           && PFNeutralHadronIso < 2.9 + 0.04 * pho->pt

          )) 
      pass = false;         
  }

  return pass;
}

#endif

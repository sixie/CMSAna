#ifndef LEPTONIDCUTS_HH
#define LEPTONIDCUTS_HH

#include "CMSAna/Utils/CommonDefs.hh"
#include "CMSAna/Utils/CommonTools.hh"
#include "CMSAna/Utils/IsolationPileupCorrections.hh"
#include "CMSAna/DataTree/interface/TElectron.hh"
#include "CMSAna/DataTree/interface/TMuon.hh"
#include "CMSAna/DataTree/interface/DataTreeDefs.hh"
#include <cassert>
#include "TMath.h"


Double_t ComputeMuonPFIsoRings( const cmsana::TMuon *muon,
                                TClonesArray *pfCandidates, 
                                Double_t rho, UInt_t DataEra,
                                UInt_t pfType,
                                Double_t minDR,
                                Double_t maxDR,
                                Bool_t printDebug = kFALSE);

Double_t ComputeElePFIsoRings( const cmsana::TElectron *ele, 
                               TClonesArray *pfCandidates, 
                               Double_t rho, UInt_t DataEra,
                               UInt_t pfIsoType,
                               Double_t minDR,
                               Double_t maxDR,
                               Bool_t printDebug = kFALSE);

Bool_t PassMuonIDVeto(const cmsana::TMuon *muon);



Bool_t PassEleSimpleCutsLooseID( const cmsana::TElectron *ele,   
                                 Bool_t printDebug = kFALSE);  

Bool_t PassEleSimpleCutsTightID( const cmsana::TElectron *ele,   
                                 Bool_t printDebug = kFALSE);  

Bool_t PassEleSimpleCutsVeto( const cmsana::TElectron *ele,   
                               TClonesArray *pfCandidates, 
                               Double_t rho, UInt_t DataEra, 
                               Bool_t printDebug = kFALSE);  



//=== FUNCTION DEFINITIONS ======================================================================================




//--------------------------------------------------------------------------------------------------
Double_t ComputeMuonPFIsoRings( const cmsana::TMuon *muon,
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
    
    //************************************************
    // veto PF Candidates that match to the muon track
    //************************************************
    bool match = false;
    for (UInt_t p=0; p < muon->NMatchedPFCandidates; ++p) {
      if (muon->MatchedPFCandidateIndex[p] == k) match = true;
    }
    if (match) continue;


    Double_t dr = cmsana::deltaR(muon->eta, muon->phi, pf->eta, pf->phi);
    if (dr < minDR) continue;
    if (dr >= maxDR) continue;

    //***********************************************************************
    // Charged Iso 
    //***********************************************************************
    if (pf->q != 0) {      
      //don't include pf electrons or pf muons
      if (pf->pfType == eElectron || pf->pfType == eMuon) continue;

      if (pfIsoType == kPFChargedIso || pfIsoType == kPFIso) fIso += pf->pt;
    }
    
    //***********************************************************************
    // Gamma Iso 
    //***********************************************************************
    else if (pf->pfType == eGamma) {
      if( pf->pt > 0.5 ) {
        if (pfIsoType == kPFGammaIso || pfIsoType == kPFIso) fIso += pf->pt;
      }
    }
    
    //***********************************************************************
    // Other Neutrals
    //***********************************************************************
    else {
      if( pf->pt > 0.5 ) {
	if (pfIsoType == kPFNeutralHadronIso || pfIsoType == kPFIso) fIso += pf->pt;
      }
    }
  }

  if(printDebug) cout << "rho: " << rho << endl;

  UInt_t EAType = 0;
  if (pfIsoType == kPFGammaIso) {
    if (minDR == 0.0 && maxDR == 0.1) EAType = kMuGammaIsoDR0p0To0p1;
    else if (minDR == 0.1 && maxDR == 0.2) EAType = kMuGammaIsoDR0p1To0p2;
    else if (minDR == 0.2 && maxDR == 0.3) EAType = kMuGammaIsoDR0p2To0p3;
    else if (minDR == 0.3 && maxDR == 0.4) EAType = kMuGammaIsoDR0p3To0p4;
    else if (minDR == 0.4 && maxDR == 0.5) EAType = kMuGammaIsoDR0p4To0p5;
    else EAType = kMuNoCorrection;
  } else if (pfIsoType == kPFNeutralHadronIso) {
    if (minDR == 0.0 && maxDR == 0.1) EAType = kMuNeutralHadronIsoDR0p0To0p1;
    else if (minDR == 0.1 && maxDR == 0.2) EAType = kMuNeutralHadronIsoDR0p1To0p2;
    else if (minDR == 0.2 && maxDR == 0.3) EAType = kMuNeutralHadronIsoDR0p2To0p3;
    else if (minDR == 0.3 && maxDR == 0.4) EAType = kMuNeutralHadronIsoDR0p3To0p4;
    else if (minDR == 0.4 && maxDR == 0.5) EAType = kMuNeutralHadronIsoDR0p4To0p5;
    else EAType = kMuNoCorrection;
  } else if (pfIsoType == kPFIso) {
    if (minDR == 0.0 && maxDR == 0.3) EAType = kMuGammaAndNeutralHadronIso03;
    else if (minDR == 0.0 && maxDR == 0.4) EAType = kMuGammaAndNeutralHadronIso04;
  } else {
    EAType = kMuNoCorrection;
  }

  double iso = fmax(0.0, fIso - rho*MuonEffectiveArea(EAType, muon->eta,DataEra));
  
  if( printDebug ) { 
    cout << "uncorrected iso: " << fIso
	 << "\tcorrected iso: " << iso
	 << endl;
  }

  return iso;

}




//--------------------------------------------------------------------------------------------------
Double_t ComputeElePFIsoRings( const cmsana::TElectron *ele, 
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

    //*************************************************************
    // veto PF Candidates that match to the electron SC or GsfTrack
    //*************************************************************
    bool match = false;
    for (UInt_t p=0; p < ele->NSCMatchedPFCandidates; ++p) {
      if (ele->SCMatchedPFCandidateIndex[p] == k) match = true;
    }
    for (UInt_t p=0; p < ele->NGsfTrackMatchedPFCandidates; ++p) {
      if (ele->GsfTrackMatchedPFCandidateIndex[p] == k) match = true;
    }
    if (match) continue;


    Double_t dr = cmsana::deltaR(ele->eta, ele->phi, pf->eta, pf->phi);
    if (dr < minDR) continue;
    if (dr >= maxDR) continue;

    if(printDebug) { 
      cout << "pf :: type: " << pf->pfType << "\tpt: " << pf->pt << "\tdR: " << dr;
      cout << "\tdZ: " << pf->dz;
      cout << endl;
    }

    //***********************************************************************
    // Charged Iso 
    //***********************************************************************
    if (pf->q != 0) {

      // Veto any PFmuon, or PFEle
      if ( pf->pfType == eElectron || pf->pfType == eMuon) {
        continue;
      }

      // Footprint Veto
      if (fabs(ele->scEta) > 1.479 && dr < 0.015) continue;

      if( printDebug) cout << "charged:: pt: " << pf->pt 
                           << "\ttype: " << pf->pfType
                           << endl;

      if (pfIsoType == kPFChargedIso || pfIsoType == kPFIso) fIso += pf->pt;
    }

    //***********************************************************************
    // Gamma Iso 
    //***********************************************************************
    else if (pf->pfType == eGamma) {

      if (fabs(ele->scEta) < 1.479) {
        if (fabs(ele->eta - pf->eta) < 0.015) continue;
      }

      if (fabs(ele->scEta) > 1.479) {
        if (cmsana::deltaR(ele->eta,ele->phi, pf->eta, pf->phi) < 0.08) continue;
      }
      if( printDebug) cout << "gamma:: " << pf->pt << " " 
                           << dr << endl;
      if (pfIsoType == kPFGammaIso || pfIsoType == kPFIso) fIso += pf->pt;
    }

    //***********************************************************************
    // Other Neutrals
    //***********************************************************************
    else {
      if( printDebug ) cout << "neutral:: " << pf->pt << " " 
                           << dr << endl;
      if (pfIsoType == kPFNeutralHadronIso || pfIsoType == kPFIso) fIso += pf->pt;
    }
  }

  if( printDebug ) cout << "rho: " << rho << endl;

  UInt_t EAType = 0;
  if (pfIsoType == kPFGammaIso) {
    if (minDR == 0.0 && maxDR == 0.1) EAType = kEleGammaIsoDR0p0To0p1;
    else if (minDR == 0.1 && maxDR == 0.2) EAType = kEleGammaIsoDR0p1To0p2;
    else if (minDR == 0.2 && maxDR == 0.3) EAType = kEleGammaIsoDR0p2To0p3;
    else if (minDR == 0.3 && maxDR == 0.4) EAType = kEleGammaIsoDR0p3To0p4;
    else if (minDR == 0.4 && maxDR == 0.5) EAType = kEleGammaIsoDR0p4To0p5;
    else EAType = kEleNoCorrection;
  } else if (pfIsoType == kPFNeutralHadronIso) {
    if (minDR == 0.0 && maxDR == 0.1) EAType = kEleNeutralHadronIsoDR0p0To0p1;
    else if (minDR == 0.1 && maxDR == 0.2) EAType = kEleNeutralHadronIsoDR0p1To0p2;
    else if (minDR == 0.2 && maxDR == 0.3) EAType = kEleNeutralHadronIsoDR0p2To0p3;
    else if (minDR == 0.3 && maxDR == 0.4) EAType = kEleNeutralHadronIsoDR0p3To0p4;
    else if (minDR == 0.4 && maxDR == 0.5) EAType = kEleNeutralHadronIsoDR0p4To0p5;
    else EAType = kEleNoCorrection;
  } else if (pfIsoType == kPFIso) {
    if (minDR == 0.0 && maxDR == 0.3) EAType = kEleGammaAndNeutralHadronIso03;
    else if (minDR == 0.0 && maxDR == 0.4) EAType = kEleGammaAndNeutralHadronIso04;
  } else {
    EAType = kEleNoCorrection;
  }

  double iso = fmax(0.0, fIso - rho*ElectronEffectiveArea(EAType, ele->scEta,DataEra));


  if( printDebug ) { 
    cout << "uncorrected iso: " << fIso
	 << "\tcorrected iso: " << iso
	 << endl;
  }

  return iso;

}


//--------------------------------------------------------------------------------------------------
Bool_t PassMuonIDVeto(const cmsana::TMuon *muon )
{

  if(muon->nTkHits	  < 11)    return false;
  if(muon->nPixHits	  < 1)     return false;
  if(fabs(muon->dz)       > 0.5)   return false;

  Bool_t isGlobal  = ((muon->typeBits & kGlobal) == kGlobal) && (muon->muNchi2 < 10) && (muon->nMatch > 1) && (muon->nValidHits > 0);
  Bool_t isTracker = ((muon->typeBits & kTracker) == kTracker) && (muon->qualityBits & kTMLastStationTight);

  if(!isGlobal && !isTracker) return false;

  return true;
}


Bool_t PassEleSimpleCutsLooseID( const cmsana::TElectron *ele,   
                                 Bool_t printDebug) {
  
  Bool_t pass = kTRUE;
  if (ele->isEB) {
    if (!(
          fabs(ele->deltaEtaIn) < 0.007
          && fabs(ele->deltaPhiIn) < 0.15
          && ele->sigiEtaiEta < 0.01
          && ele->HoverE < 0.12
          && fabs(ele->d0) < 0.02
//          && fabs(ele->dz) < 0.2
          && fabs(1.0/(ele->EcalEnergy) - 1.0/ele->pIn) < 0.05
          && !(ele->isConv)
          && ele->nExpHitsInner <= 1
          )
      ) pass = kFALSE;
  } else {
    if (!(
          fabs(ele->deltaEtaIn) < 0.01
          && fabs(ele->deltaPhiIn) < 0.10
          && ele->sigiEtaiEta < 0.03
          && ele->HoverE < 0.10
          && fabs(ele->d0) < 0.02
//          && fabs(ele->dz) < 0.2
          && fabs(1.0/(ele->EcalEnergy) - 1.0/ele->pIn) < 0.05
          && !(ele->isConv)
          && ele->nExpHitsInner <= 1
          )
      ) pass = kFALSE;
  }

  return pass;

}


Bool_t PassEleSimpleCutsTightID( const cmsana::TElectron *ele,   
                                 Bool_t printDebug) {
  
  Bool_t pass = kTRUE;
  if (ele->isEB) {
    if (!(
          fabs(ele->deltaEtaIn) < 0.004
          && fabs(ele->deltaPhiIn) < 0.03
          && ele->sigiEtaiEta < 0.01
          && ele->HoverE < 0.12
          && fabs(ele->d0) < 0.02
//          && fabs(ele->dz) < 0.2
          && fabs(1.0/(ele->EcalEnergy) - 1.0/ele->pIn) < 0.05
          && !(ele->isConv)
          && ele->nExpHitsInner <= 0
          )
      ) pass = kFALSE;
  } else {
    if (!(
          fabs(ele->deltaEtaIn) < 0.005
          && fabs(ele->deltaPhiIn) < 0.02
          && ele->sigiEtaiEta < 0.03
          && ele->HoverE < 0.10
          && fabs(ele->d0) < 0.02
//          && fabs(ele->dz) < 0.2
          && fabs(1.0/(ele->EcalEnergy) - 1.0/ele->pIn) < 0.05
          && !(ele->isConv)
          && ele->nExpHitsInner <= 0
          )
      ) pass = kFALSE;
  }

  return pass;

}

Bool_t PassEleSimpleCutsVeto( const cmsana::TElectron *ele,   
                              TClonesArray *pfCandidates, 
                              Double_t rho, UInt_t DataEra, 
                              Bool_t printDebug) {
  
  Bool_t pass = kTRUE;
  if (ele->isEB) {
    if (!(
          fabs(ele->deltaEtaIn) < 0.007
          && fabs(ele->deltaPhiIn) < 0.8
          && ele->sigiEtaiEta < 0.01
          && ele->HoverE < 0.15
          && fabs(ele->d0) < 0.04
          && fabs(ele->dz) < 0.2
          && ComputeElePFIsoRings( ele, pfCandidates, rho, DataEra, kPFIso, 
                                   0.0, 0.3, printDebug)/ele->pt < 0.15
          )
      ) pass = kFALSE;
  } else {
    if (!(
          fabs(ele->deltaEtaIn) < 0.01
          && fabs(ele->deltaPhiIn) < 0.70
          && ele->sigiEtaiEta < 0.03
          && fabs(ele->d0) < 0.04
          && fabs(ele->dz) < 0.2
          && ComputeElePFIsoRings( ele, pfCandidates, rho, DataEra, kPFIso,
                                   0.0, 0.3, printDebug)/ele->pt < 0.15
          )
      ) pass = kFALSE;
  }

  return pass;

}

#endif

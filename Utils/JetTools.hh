#ifndef CMSANA_UTILS_PHOTONTOOLS_HH
#define CMSANA_UTILS_PHOTONTOOLS_HH

#include <TMath.h>
#include <TClonesArray.h>          
#include "CMSAna/DataTree/interface/TJet.hh"
#include "CMSAna/DataTree/interface/TGenJet.hh"
#include "CMSAna/ObjectStudies/interface/JetTree.h"

void FillJetTree(cmsana::JetTree *jetTree,
                 const cmsana::TJet *jet, 
                 const cmsana::TGenJet *genjet, 
                 Double_t rho, UInt_t DataEra,
                 cmsana::TEventInfo *info,
                 Float_t weight = 1.0
                 
  ) {

  if (genjet) {
    jetTree->fJetGenPt = genjet->pt;
    jetTree->fJetGenEta = genjet->eta;
    jetTree->fJetGenPhi = genjet->phi;
  } else {
    jetTree->fJetGenPt = -99.0;
    jetTree->fJetGenEta = -99.0;
    jetTree->fJetGenPhi = -99.0;
  }

  if (jet) {
    jetTree->fWeight = weight;
    jetTree->fRunNumber = info->runNum;
    jetTree->fLumiSectionNumber = info->lumiSec;
    jetTree->fEventNumber = info->evtNum;
    jetTree->fJetEventNumberParity = (info->evtNum % 2 == 0);
    jetTree->fJetRawPt = jet->rawPt; 
    jetTree->fJetPt = jet->pt; 
    jetTree->fJetEta = jet->eta; 
    jetTree->fJetPhi = jet->phi; 
    jetTree->fJetTriggerBit = 0;
    jetTree->fRho = rho; 
    jetTree->fNVertices = info->nGoodPV; 
    jetTree->fJetTrackCountingHighEffBJetTagsDisc = jet->TrackCountingHighEffBJetTagsDisc; 
    jetTree->fJetTrackCountingHighPurBJetTagsDisc = jet->TrackCountingHighPurBJetTagsDisc; 
    jetTree->fJetSimpleSecondaryVertexHighPurBJetTagsDisc = jet->SimpleSecondaryVertexHighPurBJetTagsDisc; 
    jetTree->fJetSimpleSecondaryVertexHighEffBJetTagsDisc = jet->SimpleSecondaryVertexHighEffBJetTagsDisc; 
    jetTree->fJetCombinedSecondaryVertexBJetTagsDisc = jet->CombinedSecondaryVertexBJetTagsDisc; 
    jetTree->fJetCombinedSecondaryVertexMVABJetTagsDisc = jet->CombinedSecondaryVertexMVABJetTagsDisc; 
    jetTree->fJetProbabilityBJetTagsDisc = jet->JetProbabilityBJetTagsDisc; 
    jetTree->fJetBProbabilityBJetTagsDisc = jet->JetBProbabilityBJetTagsDisc; 
    jetTree->fJetArea = jet->JetArea; 
    jetTree->fJetMatchedPdgId = jet->matchedPdgId;
  }

  //***********************
  //Fill Jet
  //***********************
  jetTree->tree_->Fill();

}



#endif

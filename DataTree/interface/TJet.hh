#ifndef CMSANA_NTUPLER_TJET_HH
#define CMSANA_NTUPLER_TJET_HH

#include <TObject.h>

namespace cmsana
{
  class TJet : public TObject
  {
    public:
      TJet(){
        pt = 0; 
        eta = -999;
        phi = -999;
        mass = -999;
        rawPt = 0;
        L1JECScale = 1;
        L2JECScale = 1;
        L3JECScale = 1;        
        hltMatchBits = 0;
        NConstituents = 0;
        NeutralHadronFraction = -999;
        NeutralEMFraction = -999;
        ChargedHadronFraction = -999;
        ChargedEMFraction = -999;
        TrackCountingHighEffBJetTagsDisc = -999;
        TrackCountingHighPurBJetTagsDisc = -999;
        SoftElectronByPtBJetTagsDisc = -999;
        SoftElectronByIP3dBJetTagsDisc = -999;
        SoftMuonByPtBJetTagsDisc = -999;
        SoftMuonByIP3dBJetTagsDisc = -999;
        SoftMuonBJetTagsDisc = -999;
        SimpleSecondaryVertexHighPurBJetTagsDisc = -999;
        SimpleSecondaryVertexHighEffBJetTagsDisc = -999;
        CombinedSecondaryVertexBJetTagsDisc = -999;
        CombinedSecondaryVertexMVABJetTagsDisc = -999;
        JetProbabilityBJetTagsDisc = -999;
        JetBProbabilityBJetTagsDisc = -999;
        JetArea = -999;        
        matchedPdgId = -999;
      }
      ~TJet(){}

      Float_t pt, eta, phi, mass;  // kinematics
      Float_t rawPt;
      Float_t L1JECScale;
      Float_t L2JECScale;
      Float_t L3JECScale;

      ULong_t hltMatchBits;   // bits from matching with HLT primitives
      UInt_t NConstituents;
      Float_t NeutralHadronFraction;
      Float_t NeutralEMFraction;
      Float_t ChargedHadronFraction;
      Float_t ChargedEMFraction;

      Float_t TrackCountingHighEffBJetTagsDisc; //btag discriminator
      Float_t TrackCountingHighPurBJetTagsDisc;
      Float_t SoftElectronByPtBJetTagsDisc;
      Float_t SoftElectronByIP3dBJetTagsDisc;
      Float_t SoftMuonByPtBJetTagsDisc;
      Float_t SoftMuonByIP3dBJetTagsDisc;
      Float_t SoftMuonBJetTagsDisc;
      Float_t SimpleSecondaryVertexHighPurBJetTagsDisc;
      Float_t SimpleSecondaryVertexHighEffBJetTagsDisc;
      Float_t CombinedSecondaryVertexBJetTagsDisc;
      Float_t CombinedSecondaryVertexMVABJetTagsDisc;
      Float_t JetProbabilityBJetTagsDisc;
      Float_t JetBProbabilityBJetTagsDisc;
      Float_t JetArea;

      Int_t   matchedPdgId;

    ClassDef(TJet,2)
  };
}
#endif

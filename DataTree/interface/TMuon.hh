#ifndef CMSANA_NTUPLER_TMUON_HH
#define CMSANA_NTUPLER_TMUON_HH

#include <TObject.h>

namespace cmsana 
{
  class TMuon : public TObject
  {
    public:
      TMuon() {
        pt = -999;
        eta = -999;
        phi = -999;
        pterr = -999;                
        q = 0;		             
        typeBits = 0;	             
        qualityBits = 0;
        nValidHits = 0;
        nTkHits = 0;	
        nPixHits = 0;	
        nMatch = 0;      
        trkLayers = 0;    
        tkNchi2 = -999;	
        muNchi2 = -999;	
        hltMatchBits = 0;
        isMCReal = 0;	
        trkIso03 = -999;	
        emIso03 = -999;	
        hadIso03 = -999;	
        hoIso03 = -999;	
        pfIso04ChargedHadron = -999;   
        pfIso04ChargedParticle = -999; 
        pfIso04NeutralHadron = -999;   
        pfIso04Photon = -999;          
        pfIso04NeutralHadronHighThreshold = -999; 
        pfIso04PhotonHighThreshold = -999;        
        pfIso04PU = -999;      
        d0 = -999;
        dz = -999;            
        ip3d = -999;
        ip3dSig = -999;        
        ip3dBS = -999;
        ip3dSigBS = -999;   
        NMatchedPFCandidates = 0;
        for(int k=0; k<10;++k) MatchedPFCandidateIndex[k] = 0;
      }
      ~TMuon(){} 
  
      Float_t pt, eta, phi;           // kinematics
      Float_t pterr;                  // pt error
      Int_t q;		              // charge

      UInt_t typeBits;	              // global muon, tracker muon, or standalone muon
      UInt_t qualityBits;             // bits for various muon quality criteria

      Int_t nValidHits;	              // number of valid hits in muon system
      UInt_t nTkHits;	              // number of inner tracker hits
      UInt_t nPixHits;	              // number of pixel hits
      UInt_t nMatch;                  // number of muon chambers matched to segments
      UInt_t trkLayers;               // number of tracker layers with measured hits
      Float_t tkNchi2;	              // track chi^2/ndf 
      Float_t muNchi2;	              // global muon chi^2/ndf

      ULong_t hltMatchBits;           // bits for matching with HLT primitives 
      Int_t  isMCReal;	              // MC Truth Matching

      Float_t trkIso03;	              // track isolation
      Float_t emIso03;	              // ECAL-based isolation
      Float_t hadIso03;	              // HCAL-based isolation
      Float_t hoIso03;	              // HO-based isolation

      Float_t pfIso04ChargedHadron;   // PF Isolation
      Float_t pfIso04ChargedParticle; // PF Isolation
      Float_t pfIso04NeutralHadron;   // PF Isolation
      Float_t pfIso04Photon;          // PF Isolation
      Float_t pfIso04NeutralHadronHighThreshold; // PF Isolation
      Float_t pfIso04PhotonHighThreshold;        // PF Isolation
      Float_t pfIso04PU;
      
      Float_t d0, dz;                 // impact parameter
      Float_t ip3d, ip3dSig;          // IP3D PV
      Float_t ip3dBS, ip3dSigBS;      // IP3D BeamSpot constrained PV

      UInt_t NMatchedPFCandidates;
      UInt_t MatchedPFCandidateIndex[10];

      ClassDef(TMuon,2)
  };  
}
#endif

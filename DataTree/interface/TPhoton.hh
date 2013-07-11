#ifndef CMSANA_NTUPLER_TPHOTON_HH
#define CMSANA_NTUPLER_TPHOTON_HH

#include <TObject.h>
#include <vector>

namespace cmsana
{
  class TPhoton : public TObject
  {
    public:
      TPhoton(){
        pt = -999;
        eta = -999;
        phi = -999;
        scEt = -999;
        scEta = -999;
        scPhi = -999;
        trkIso03Hollow = -999;
        trkIso03Solid = -999;
        emIso03 = -999;
        hadIso03 = -999;
        HoverE = -999;
        HoverESingleTower = -999;
        R9 = -999;
        sigiEtaiEta = -999;       
        hltMatchBits = 0;
        scID = 0;        
        hasPixelSeed = false;
        isConversion = false;
        passEleVeto = false;

        SCRawEnergy = -999;
        PreShowerOverRaw = -999;
        sigiPhiiPhi = -999;         
        sigiEtaiPhi = -999;
        SCEtaWidth = -999;
        SCPhiWidth = -999;
        NClusters = 0;
        E5x5 = -999;
        EtaSeed = -999;
        PhiSeed = -999;
        ESeed = -999;
        E2x2Seed = -999;
        E3x3Seed = -999;
        E5x5Seed = -999;
        EMaxSeed = -999;
        E2ndSeed = -999;
        ETopSeed = -999;
        EBottomSeed = -999;
        ELeftSeed = -999;
        ERightSeed = -999;
        E2x5MaxSeed = -999;
        E2x5TopSeed = -999;
        E2x5BottomSeed = -999;
        E2x5LeftSeed = -999;
        E2x5RightSeed = -999;
        IEtaSeed = -999;
        IPhiSeed = -999;
        EtaCrySeed = -999;
        PhiCrySeed = -999;
        NMatchedPFCandidates = 0;
        for(int k=0; k<10;++k) MatchedPFCandidateIndex[k] = 0;
      }
      ~TPhoton(){}

      Float_t pt, eta, phi; 	              // kinematics
      Float_t scEt, scEta, scPhi;             // supercluster
      Float_t trkIso03Hollow, trkIso03Solid;  // track isolation
      Float_t emIso03;                        // ECAL-based isolation
      Float_t hadIso03;                       // HCAL-based isolation
      Float_t HoverE;		              // H/E
      Float_t HoverESingleTower;              // H/E Single Tower
      Float_t R9;		              // ratio of energies in 3x3 to SC
      Float_t sigiEtaiEta;                    // eta-width of shower in number of crystals

      UInt_t  hltMatchBits;  	              // bits from matching with HLT primitives
      UInt_t  scID;                           // supercluster ID (for matching to electron superclusters)

      Bool_t  hasPixelSeed;                   // supercluster has pixel seed?
      Bool_t  isConversion;
      Bool_t  passEleVeto;

      //regression related variables
      Float_t SCRawEnergy;
      Float_t PreShowerOverRaw;
      Float_t sigiPhiiPhi;         
      Float_t sigiEtaiPhi;
      Float_t SCEtaWidth;
      Float_t SCPhiWidth;
      UInt_t  NClusters;
      Float_t E5x5;
      Float_t EtaSeed;
      Float_t PhiSeed;
      Float_t ESeed;
      Float_t E2x2Seed;
      Float_t E3x3Seed;
      Float_t E5x5Seed;
      Float_t EMaxSeed;
      Float_t E2ndSeed;
      Float_t ETopSeed;
      Float_t EBottomSeed;
      Float_t ELeftSeed;
      Float_t ERightSeed;
      Float_t E2x5MaxSeed;
      Float_t E2x5TopSeed;
      Float_t E2x5BottomSeed;
      Float_t E2x5LeftSeed;
      Float_t E2x5RightSeed;
      Float_t IEtaSeed;
      Float_t IPhiSeed;
      Float_t EtaCrySeed;
      Float_t PhiCrySeed;

      //for match to pf candidates
      UInt_t NMatchedPFCandidates;
      UInt_t MatchedPFCandidateIndex[10];

    ClassDef(TPhoton,2)
  };  
}
#endif

#ifndef CMSANA_NTUPLER_TPFCANDIDATE_HH
#define CMSANA_NTUPLER_TPFCANDIDATE_HH

#include <TObject.h>

namespace cmsana
{
  class TPFCandidate : public TObject
  {
    public:

      TPFCandidate(){
        pt = -999;
        eta = -999;
        phi = -999;
        e = -999;
        q = -999;
        dz = -999;        
        IsPFNoPU = true;
        pfType = 0;
        matchedObjectType = 0;
        matchedObjectIndex = 0;
      }
      ~TPFCandidate(){}
    
      Float_t pt, eta, phi, e;     // kinematics      
      Int_t q;                     // charge    
      Float_t dz;                  // dz to the primary vertex
      Bool_t   IsPFNoPU;
      UInt_t   pfType;              
      UInt_t   matchedObjectType;   // matched to reco ele,mu,photons
      UInt_t   matchedObjectIndex;  // index of the matched object

    ClassDef(TPFCandidate,1)
  };
}
#endif

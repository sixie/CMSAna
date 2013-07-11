#ifndef CMSANA_NTUPLER_TGENJET_HH
#define CMSANA_NTUPLER_TGENJET_HH

#include <TObject.h>

namespace cmsana
{
  class TGenJet : public TObject
  {
    public:
      TGenJet(){
        pt = 0; 
        eta = -999;
        phi = -999;
        mass = -999;
        matchedPdgId = 0;
      }
      ~TGenJet(){}

      Float_t pt, eta, phi, mass;  // kinematics
      Int_t matchedPdgId;

    ClassDef(TGenJet,1)
  };
}
#endif

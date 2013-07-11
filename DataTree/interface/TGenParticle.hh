#ifndef CMSANA_NTUPLER_TGenParticle_HH
#define CMSANA_NTUPLER_TGenParticle_HH

#include <TObject.h>

namespace cmsana 
{
  class TGenParticle : public TObject
  {
    public:
      TGenParticle(){
        pt = -999;
        eta= -999;
        phi= -999;
        mass= -999;
        pdgid = -999;
        status = -999;
        motherPdgID = -999;
      }
      ~TGenParticle(){} 
  
      Float_t pt, eta, phi, mass;     // kinematics
      Int_t pdgid;
      Int_t status;
      Int_t motherPdgID;

      ClassDef(TGenParticle,1)
        };  
}
#endif

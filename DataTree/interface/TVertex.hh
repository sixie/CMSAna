#ifndef CMSANA_NTUPLER_TVERTEX_HH
#define CMSANA_NTUPLER_TVERTEX_HH

#include <TObject.h>

namespace cmsana
{
  class TVertex : public TObject
  {
    public:

      TVertex(){
        x = -999;
        y = -999;
        z = -999;
        isGoodVertex = false;
      }
      ~TVertex(){}
    
      Float_t x, y, z;
      Bool_t  isGoodVertex;

    ClassDef(TVertex,1)
  };
}
#endif

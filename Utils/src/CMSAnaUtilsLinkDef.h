#ifndef CMSANA_UTILS_LINKDEF_H
#define CMSANA_UTILS_LINKDEF_H

#include "CMSAna/Utils/interface/SimpleTable.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace cmsana;

#pragma link C++ class cmsana::SimpleTable+;
#pragma link C++ class cmsana::SimpleTable::MyParameter+;
#endif

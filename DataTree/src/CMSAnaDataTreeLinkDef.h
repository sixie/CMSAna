#ifndef CMSANA_NTUPLER_LINKDEF_H
#define CMSANA_NTUPLER_LINKDEF_H
#include "CMSAna/DataTree/interface/TEventInfo.hh"
#include "CMSAna/DataTree/interface/TGenParticle.hh"
#include "CMSAna/DataTree/interface/TGenJet.hh"
#include "CMSAna/DataTree/interface/TMuon.hh"
#include "CMSAna/DataTree/interface/TElectron.hh"
#include "CMSAna/DataTree/interface/TJet.hh"
#include "CMSAna/DataTree/interface/TPhoton.hh"
#include "CMSAna/DataTree/interface/TPFCandidate.hh"
#include "CMSAna/DataTree/interface/TVertex.hh"
#include "CMSAna/DataTree/interface/Types.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace cmsana;

#pragma link C++ class cmsana::TEventInfo+;
#pragma link C++ class cmsana::TGenParticle+;
#pragma link C++ class cmsana::TGenJet+;
#pragma link C++ class cmsana::TMuon+;
#pragma link C++ class cmsana::TElectron+;
#pragma link C++ class cmsana::TJet+;
#pragma link C++ class cmsana::TPhoton+;
#pragma link C++ class cmsana::TPFCandidate+;
#pragma link C++ class cmsana::TVertex+;
#pragma link C++ typedef cmsana::FourVector;
#pragma link C++ typedef cmsana::FourVectorM;
#pragma link C++ typedef cmsana::FourVectorE;
#endif

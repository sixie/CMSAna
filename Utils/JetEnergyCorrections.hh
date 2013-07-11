#ifndef JETENERGYCORRECTIONS_HH
#define JETENERGYCORRECTIONS_HH

#include "CMSAna/DataTree/interface/Types.h"
// #include "CMSAna/Utils/CommonDefs.hh"
// #include "CMSAna/Utils/CommonTools.hh"
#include "CMSAna/DataTree/interface/TJet.hh"
// #include "CMSAna/DataTree/interface/DataTreeDefs.hh"
#include "CMSAna/JetEnergyCorrections/interface/FactorizedJetCorrector.h"

#include <cassert>
#include <iostream>
#include "TMath.h"

Double_t JetEnergyCorrectionFactor( const cmsana::TJet *jet,
                                    cmsana::FactorizedJetCorrector jetcorrector,
                                    Bool_t printDebug = kFALSE);



//=== FUNCTION DEFINITIONS ======================================================================================





//--------------------------------------------------------------------------------------------------
Double_t JetEnergyCorrectionFactor( const cmsana::TJet *jet,
                                    cmsana::FactorizedJetCorrector *jetcorrector,  
                                    double rho,
                                    Bool_t printDebug) {
  cmsana::FourVectorM temp;
  temp.SetPt(jet->rawPt);
  temp.SetEta(jet->eta);
  temp.SetPhi(jet->phi);
  temp.SetM(jet->mass);

  jetcorrector->setJetEta(temp.Eta());
  jetcorrector->setJetPt(temp.Pt());
  jetcorrector->setJetPhi(temp.Phi());
  jetcorrector->setJetE(temp.E());
  jetcorrector->setRho(rho);
  jetcorrector->setJetA(jet->JetArea);

  std::vector<float> corrections;
  corrections = jetcorrector->getSubCorrections();

  if (printDebug) cout << "Computing Jet Energy Corrections for jet with raw momentum: " << jet->rawPt << " " << jet->eta << " " << jet->phi << "\n";

  Double_t cumulativeCorrection = 1.0;
  for (UInt_t j=0; j<corrections.size(); ++j) {
    Double_t currentCorrection = corrections.at(j)/cumulativeCorrection;
    cumulativeCorrection = corrections.at(j);
    if (printDebug) cout << "Correction Level " << j << " : current correction = " << currentCorrection << " , cumulative correction = " << cumulativeCorrection << "\n";
  }
  if (printDebug) cout << "Final Cumulative Correction: " << cumulativeCorrection << "\n";
  
  return cumulativeCorrection;

}

#endif

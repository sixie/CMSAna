#include "CMSAna/JetEnergyCorrections/interface/JetCorrectorParameters.h"
#include "CMSAna/JetEnergyCorrections/interface/FactorizedJetCorrector.h"
#include "CMSAna/JetEnergyCorrections/interface/SimpleJetCorrector.h"
#include <vector>
 
cmsana::JetCorrectorParameters corr;
cmsana::JetCorrectorParameters::Definitions def;
cmsana::JetCorrectorParameters::Record record;
std::vector<cmsana::JetCorrectorParameters> corrv;
std::vector<cmsana::JetCorrectorParameters::Record> recordv;
cmsana::JetCorrectorParametersCollection coll;
cmsana::JetCorrectorParametersCollection::pair_type pair_type;
cmsana::JetCorrectorParametersCollection::collection_type colltype;
std::vector<cmsana::JetCorrectorParametersCollection> collv;

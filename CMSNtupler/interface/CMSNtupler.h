// system include files
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TString.h"

using namespace std;
using namespace edm;


//*****************************************************************************
//BACON Data Formats
//*****************************************************************************
#include "CMSAna/DataTree/interface/DataTreeDefs.hh"
#include "CMSAna/DataTree/interface/TEventInfo.hh"
#include <TClonesArray.h>
#include "CMSAna/DataTree/interface/TGenParticle.hh"
#include "CMSAna/DataTree/interface/TGenJet.hh"
#include "CMSAna/DataTree/interface/TElectron.hh"
#include "CMSAna/DataTree/interface/TJet.hh"
#include "CMSAna/DataTree/interface/TMuon.hh"
#include "CMSAna/DataTree/interface/TPhoton.hh"
#include "CMSAna/DataTree/interface/TPFCandidate.hh"
#include "CMSAna/DataTree/interface/TVertex.hh"


class CMSNtupler : public edm::EDAnalyzer {
   public:
      explicit CMSNtupler(const edm::ParameterSet&);
      ~CMSNtupler();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
                   
   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
      double DeltaPhi(double v1, double v2); 
      double GetDeltaR(double eta1, double eta2, double phi1, double phi2); 

  
    //*****************************************************************
    // Input Collection Tags
    //*****************************************************************
    string fGenParticlesSrcName;
    string fGenJetsSrcName;
    string fPrimaryVerticesSrcName;
    string fPrimaryVerticesBSSrcName;
    string fMuonsSrcName;
    string fElectronsSrcName;
    string fPhotonsSrcName;
    string fPFJetsSrcName;
    string fPFCandidatesSrcName;

    //***************************************************************************
    //OUTPUT
    //***************************************************************************
    cmsana::TEventInfo      fEventInfo;       // general event information
    TClonesArray           *fGenParticleArr;  // genparticle array
    TClonesArray           *fGenJetArr;       // genjet array
    TClonesArray           *fElectronArr;     // electron array
    TClonesArray           *fMuonArr;         // muon array
    TClonesArray           *fPFJetArr;        // particle flow jet array
    TClonesArray           *fPhotonArr;       // photon array
    TClonesArray           *fPFCandidateArr;  // PFCandidate array
    TClonesArray           *fVertexArr;       // Vertex array
    
    vector<TString>         fTriggerNamesv;       // names of triggers we're interested in 
    vector<ULong_t>         fTriggerIdsv;     // corresponding ETriggerBit value
    vector<TString>         fFirstTriggerObjectModuleNamesv; // 
    vector<TString>         fSecondTriggerObjectModuleNamesv; // 
    vector<ULong_t>         fFirstTriggerObjectIdsv;   // ETriggerObjectBit
    vector<ULong_t>         fSecondTriggerObjectIdsv;  // ETriggerObjectBit
    vector<TString>         fL1TriggerNamesv; // names of L1 triggers we're interested in 
    vector<ULong_t>         fL1TriggerIdsv;   // corresponding L1TriggerBit value
    vector<TString>         fL1SeedModuleNamesv; // names of L1 Seed Modules
    vector<ULong_t>         fL1SeedModuleIdsv;   // corresponding L1TriggerBit value
    

    TFile                  *fOutputFile;      // output file handle
    string                  fOutputName;      // output file name
    TTree*                  fEventTree;       // event tree
    int nEventsProcessed; 
    
    //***************************************************************************
    //OPTIONS
    //***************************************************************************
    bool                  fPrintDebug;
    bool                  fUseGen;          // flag whether to look at generator info
    bool                  fPrintTable;      // flag whether to print out HLT table
    bool                  fFillGenOnly;     // flag to skip reco objects       
    bool                  fFillEGRegressionVars; //flat to fill geometry dependent variables

    Double_t                fGenJetPtMin;     // minimum genjet PT
    Double_t                fElePtMin;        // minimum supercluster ET
    Double_t                fElePtMax;        // maximum supercluster ET
    Double_t                fEleEtaMin;       // minimum supercluster eta
    Double_t                fEleEtaMax;       // maximum supercluster eta
    Double_t                fMuonPtMin;       // minimum reco muon pT
    Double_t                fMuonPtMax;       // maximum reco muon pT
    Double_t                fMuonEtaMin;      // minimum reco muon eta
    Double_t                fMuonEtaMax;      // maximum reco muon eta
    Double_t                fJetPtMin;        // minimum jet PT
    Double_t                fPhotonPtMin;     // minimum photon PT
    


  protected:
    
    template <typename TYPE>
    void                     GetProduct(const std::string name, edm::Handle<TYPE> &prod,
                                        const edm::Event &event) const;    

};



//--------------------------------------------------------------------------------------------------
template <typename TYPE>
inline void CMSNtupler::GetProduct(const std::string edmname, edm::Handle<TYPE> &prod,
                                   const edm::Event &event) const
{
  // Try to access data collection from EDM file. We check if we really get just one
  // product with the given name. If not we print an error and exit.

  try {
    event.getByLabel(edm::InputTag(edmname),prod);
    if (!prod.isValid()) 
      throw edm::Exception(edm::errors::Configuration, "BaseFiller::GetProduct()\n")
        << "Cannot get collection with label " << edmname <<  std::endl;
  } catch (...) {
    edm::LogError("CMSNtupler") << "Cannot get collection with label "
                                << edmname << std::endl;
//     PrintErrorAndExit(Form("Cannot get collection with label %s", 
//                            edmname.c_str()));
  }
}


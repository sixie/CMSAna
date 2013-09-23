 
#include "CMSAna/PATNtupler/interface/PATNtupler.h"

// #include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
// #include "TrackingTools/TransientTrack/plugins/TransientTrackBuilderESProducer.h"
// #include "TrackingTools/IPTools/interface/IPTools.h"
// #include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
// #include "DataFormats/VertexReco/interface/VertexFwd.h"
// #include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"


#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"


// #include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
// #include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
// #include "RecoEgamma/EgammaTools/interface/EcalClusterLocal.h"
// #include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
// #include "DataFormats/EgammaCandidates/interface/Photon.h"
// #include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
// #include "DataFormats/JetReco/interface/PFJet.h"
// #include "DataFormats/BTauReco/interface/JetTag.h"
// #include "JetMETCorrections/Objects/interface/JetCorrector.h"
// #include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
// #include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
// #include "DataFormats/METReco/interface/PFMETCollection.h"
// #include "DataFormats/METReco/interface/PFMET.h"
// #include "DataFormats/HLTReco/interface/TriggerEvent.h"


using namespace std;
using namespace edm;
using namespace reco;


//
// constructors and destructor
//
PATNtupler::PATNtupler(const edm::ParameterSet& iConfig)

{
  
  cout << "PATNtupler initialized... " << endl;
  
  fGenParticlesSrcName = iConfig.getParameter<string> ("GenParticlesSrcName");
  fGenJetsSrcName = iConfig.getParameter<string> ("GenJetsSrcName");
  fPrimaryVerticesSrcName = iConfig.getParameter<string> ("PrimaryVerticesSrcName");
  fPrimaryVerticesBSSrcName = iConfig.getParameter<string> ("PrimaryVerticesBSSrcName");
  fMuonsSrcName = iConfig.getParameter<string> ("MuonsSrcName");
  fElectronsSrcName = iConfig.getParameter<string> ("ElectronsSrcName");
  fPhotonsSrcName = iConfig.getParameter<string> ("PhotonsSrcName");
  fPFJetsSrcName = iConfig.getParameter<string> ("PFJetsSrcName");
  fPFCandidatesSrcName = iConfig.getParameter<string> ("PFCandidatesSrcName");

//   debug_ = iConfig.getParameter<int> ("debugLevel");
  fUseGen = iConfig.getParameter<bool> ("UseGen");
  fPrintDebug = iConfig.getParameter<bool> ("PrintDebug");
  fFillEGRegressionVars = iConfig.getParameter<bool> ("FillEGRegressionVars");

  fGenJetPtMin = iConfig.getParameter<double> ("GenJetPtMin");
  fElePtMin = iConfig.getParameter<double> ("ElePtMin");
  fElePtMax = iConfig.getParameter<double> ("ElePtMax");
  fEleEtaMin = iConfig.getParameter<double> ("EleEtaMin");
  fEleEtaMax = iConfig.getParameter<double> ("EleEtaMax");
  fMuonPtMin = iConfig.getParameter<double> ("MuonPtMin");
  fMuonPtMax = iConfig.getParameter<double> ("MuonPtMax");
  fMuonEtaMin = iConfig.getParameter<double> ("MuonEtaMin");
  fMuonEtaMax = iConfig.getParameter<double> ("MuonEtaMax");
  fJetPtMin = iConfig.getParameter<double> ("JetPtMin");
  fPhotonPtMin = iConfig.getParameter<double> ("PhotonPtMin");

  //output filename
  fOutputName = iConfig.getParameter<std::string>("outputFile");

  // Don't write TObject part of the objects
  cmsana::TEventInfo::Class()->IgnoreTObjectStreamer();
  cmsana::TGenParticle::Class()->IgnoreTObjectStreamer();
  cmsana::TGenJet::Class()->IgnoreTObjectStreamer();
  cmsana::TElectron::Class()->IgnoreTObjectStreamer();
  cmsana::TMuon::Class()->IgnoreTObjectStreamer();
  cmsana::TJet::Class()->IgnoreTObjectStreamer();
  cmsana::TPhoton::Class()->IgnoreTObjectStreamer();
  cmsana::TPFCandidate::Class()->IgnoreTObjectStreamer();
  cmsana::TVertex::Class()->IgnoreTObjectStreamer();
  
  nEventsProcessed  =0;

}


PATNtupler::~PATNtupler()
{
 
  // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PATNtupler::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
    
  
  //*****************************
  //Clear Arrays
  //*****************************
  nEventsProcessed++; 
  fGenParticleArr->Clear();
  fGenJetArr->Clear();
  fElectronArr->Clear();
  fMuonArr->Clear();
  fPFJetArr->Clear();
  fPhotonArr->Clear();
  fPFCandidateArr->Clear();
  fVertexArr->Clear();

  //***************************************************************
  //Fill Gen Particles
  //***************************************************************
  double genVertexX = -999;
  double genVertexY = -999;
  double genVertexZ = -999;
  
  if (fUseGen) {

    Handle<reco::GenParticleCollection> hGenPProduct;
    GetProduct(fGenParticlesSrcName, hGenPProduct, event);  
    const reco::GenParticleCollection genParticles = *(hGenPProduct.product());  
  
    bool foundGenVertex = false;
    // loop over all genparticles
    for (reco::GenParticleCollection::const_iterator pgen = genParticles.begin();
         pgen != genParticles.end(); ++pgen) {

      //find genVertex
      if (!foundGenVertex) {
        for (unsigned int i=0; i<pgen->numberOfDaughters(); ++i) {
          const reco::Candidate *dau = pgen->daughter(i);
          if (dau) {
            genVertexX = dau->vx();
            genVertexY = dau->vy();
            genVertexZ = dau->vz();
            foundGenVertex = true;
            break;
          }
        }
      }

      //only fill leptons, neutrinos, photons, bosons
      bool fillThisParticle = false;
      if (abs(pgen->pdgId()) >= 11 && abs(pgen->pdgId()) <= 16) 
        fillThisParticle = true;
      if (abs(pgen->pdgId()) == 23 || abs(pgen->pdgId()) == 24 || abs(pgen->pdgId()) == 25
          || abs(pgen->pdgId()) == 32 || abs(pgen->pdgId()) == 33 
          || abs(pgen->pdgId()) == 35 || abs(pgen->pdgId()) == 36 
          || abs(pgen->pdgId()) == 37 || abs(pgen->pdgId()) == 443
          || abs(pgen->pdgId()) == 553) 
        fillThisParticle = true;
      
      //don't fill photons for now
//       if (abs(pgen->pdgId()) == 22 && pgen->pt() > 0.0 && fabs(pgen->eta()) < 5.0)
//         fillThisParticle = true;

      if (abs(pgen->pdgId()) == 1 || abs(pgen->pdgId()) == 2 || abs(pgen->pdgId()) == 3
          || abs(pgen->pdgId()) == 4 || abs(pgen->pdgId()) == 5 
          || abs(pgen->pdgId()) == 6)
        fillThisParticle = true;


      //debug Gen-Level information
      //cout << pgen->pdgId() << " " << pgen->status() << " " << pgen->mass() << " "
      //     << pgen->pt() << " " << pgen->eta() << " " << pgen->phi() << " ";
      //if (pgen->mother()) cout << pgen->mother()->pdgId() << " ";
      //cout << " | " << fillThisParticle << " \n";

      if (fillThisParticle) {
        TClonesArray &rGenParticleArr = *fGenParticleArr;
        assert(rGenParticleArr.GetEntries() < rGenParticleArr.GetSize());
        const Int_t index = rGenParticleArr.GetEntries();  
        new(rGenParticleArr[index]) cmsana::TGenParticle();
        cmsana::TGenParticle *pGenParticle = (cmsana::TGenParticle*)rGenParticleArr[index];
      
        pGenParticle->pt	          = pgen->pt();
        pGenParticle->eta  	  = pgen->eta();
        pGenParticle->phi  	  = pgen->phi();
        pGenParticle->mass 	  = pgen->mass();
        pGenParticle->status        = pgen->status();
        pGenParticle->pdgid	  = pgen->pdgId();
        pGenParticle->motherPdgID   = 0;    
        if (pgen->mother()) {
          pGenParticle->motherPdgID = pgen->mother()->pdgId();
          if (pgen->mother()->pdgId() == pgen->pdgId()) {
            if (pgen->mother()->mother()) {
              pGenParticle->motherPdgID = pgen->mother()->mother()->pdgId();
            }
          }
        }
      }
    } //loop over gen particles
  } //if useGen == true

  //***************************************************************
  //Fill Gen Jets
  //***************************************************************
  if (fUseGen) {
    // handle for the jet collection
    Handle<reco::GenJetCollection> hGenJetProduct;
    GetProduct(fGenJetsSrcName, hGenJetProduct, event);
    const reco::GenJetCollection inGenJets = *(hGenJetProduct.product());  

    //handle for jet flavour association
    edm::Handle<reco::JetFlavourMatchingCollection> genJetFlavMatch;
    event.getByLabel ("myAK5GenJetFlavourAssociation", genJetFlavMatch);

    // loop through all jets
    for (reco::GenJetCollection::const_iterator inGenJet = inGenJets.begin(); 
         inGenJet != inGenJets.end(); ++inGenJet) {
    
      reco::GenJetRef jetRef(hGenJetProduct, inGenJet - inGenJets.begin());    
    
      if (!(inGenJet->pt() > fGenJetPtMin)) continue;

      TClonesArray &rGenJetArr = *fGenJetArr;
      assert(rGenJetArr.GetEntries() < rGenJetArr.GetSize());
      const Int_t index = rGenJetArr.GetEntries();  
      new(rGenJetArr[index]) cmsana::TGenJet();
      cmsana::TGenJet *pGenJet = (cmsana::TGenJet*)rGenJetArr[index];
    
      pGenJet->pt	          = inGenJet->pt();
      pGenJet->eta  	  = inGenJet->eta();
      pGenJet->phi  	  = inGenJet->phi();
      pGenJet->mass 	  = inGenJet->mass();
      pGenJet->matchedPdgId = (*genJetFlavMatch)[edm::RefToBase<reco::Jet>(jetRef)].getFlavour() ;
    
    }
  }




  //****************************************************************************************************
  //Rho
  //****************************************************************************************************
  Handle<double> hRhoKt6PFJets;
  event.getByLabel(edm::InputTag("kt6PFJets","rho"),hRhoKt6PFJets);

  Handle<double> hRhoKt6PFJetsCentralChargedPileup;
  event.getByLabel(edm::InputTag("kt6PFJetsCentralChargedPileUp","rho"),hRhoKt6PFJetsCentralChargedPileup);

  Handle<double> hRhoKt6PFJetsCentralNeutral;
  event.getByLabel(edm::InputTag("kt6PFJetsCentralNeutral","rho"),hRhoKt6PFJetsCentralNeutral);

  Handle<double> hRhoKt6PFJetsCentralNeutralTight;
  event.getByLabel(edm::InputTag("kt6PFJetsCentralNeutralTight","rho"),hRhoKt6PFJetsCentralNeutralTight);



  //****************************************************************************************************
  //Handles to Muons, Electrons, Photons, Taus, ...
  //****************************************************************************************************
    Handle<vector<pat::Muon>> hMuonProduct;
    GetProduct("goodMu", hMuonProduct, event);
    const vector<pat::Muon> inMuons = *(hMuonProduct.product());

    Handle<vector<pat::Electron>> hElectronProduct;
    GetProduct("goodEl", hElectronProduct, event);
    const vector<pat::Electron> inElectrons = *(hElectronProduct.product());
 
//   Handle<reco::PhotonCollection> hPhotonProduct;
//   GetProduct(fPhotonsSrcName, hPhotonProduct, event);
//   const reco::PhotonCollection inPhotons = *(hPhotonProduct.product());  


//   //****************************************************************************************************
//   //Muons
//   //****************************************************************************************************

    for (vector<pat::Muon>::const_iterator iM = inMuons.begin(); iM != inMuons.end(); ++iM) {

      TClonesArray &rMuonArr = *fMuonArr;
      assert(rMuonArr.GetEntries() < rMuonArr.GetSize());
      const Int_t index = rMuonArr.GetEntries();  
      new(rMuonArr[index]) cmsana::TMuon();
      cmsana::TMuon *pMuon = (cmsana::TMuon*)rMuonArr[index];
      
      pMuon->pt       = iM->pt();
      pMuon->eta      = iM->eta();
      pMuon->phi      = iM->phi();
      pMuon->q        = iM->charge();

    }



   for (vector<pat::Electron>::const_iterator iE = inElectrons.begin(); 
        iE != inElectrons.end(); ++iE) {

     TClonesArray &rElectronArr = *fElectronArr;
     assert(rElectronArr.GetEntries() < rElectronArr.GetSize());
     const Int_t index = rElectronArr.GetEntries();  
     new(rElectronArr[index]) cmsana::TElectron();
     cmsana::TElectron *pElectron = (cmsana::TElectron*)rElectronArr[index];

     pElectron->pt              = iE->pt();
     pElectron->eta             = iE->eta();
     pElectron->phi             = iE->phi();
     pElectron->q               = iE->charge();

   }


  //
  // Fill event info tree
  //
  fEventInfo.runNum       = event.id().run();
  fEventInfo.evtNum       = event.id().event();
  fEventInfo.lumiSec      = event.luminosityBlock();
  fEventInfo.eventweight  = 1.0;

  fEventInfo.RhoKt6PFJets = *hRhoKt6PFJets;
  fEventInfo.RhoKt6PFJetsCentralChargedPileup = *hRhoKt6PFJetsCentralChargedPileup;
  fEventInfo.RhoKt6PFJetsCentralNeutral = *hRhoKt6PFJetsCentralNeutral;
  fEventInfo.RhoKt6PFJetsCentralNeutralTight = *hRhoKt6PFJetsCentralNeutralTight;

  ///Fill the branches
  fEventTree->Fill();

  //cleanup
//   if (lazyTools) delete lazyTools;
//   if (local) delete local;  
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
PATNtupler::beginJob()
{


  cout<<" PATNtupler::beginJob( " <<endl; 
  
  //************************************************************************************************
  // Set up output arrays
  //************************************************************************************************
  fGenParticleArr = new TClonesArray("cmsana::TGenParticle"); assert(fGenParticleArr);
  fGenJetArr      = new TClonesArray("cmsana::TGenJet");      assert(fGenJetArr);
  fElectronArr    = new TClonesArray("cmsana::TElectron");    assert(fElectronArr);
  fMuonArr        = new TClonesArray("cmsana::TMuon");        assert(fMuonArr);
  fPFJetArr       = new TClonesArray("cmsana::TJet");         assert(fPFJetArr);
  fPhotonArr      = new TClonesArray("cmsana::TPhoton");      assert(fPhotonArr);
  fPFCandidateArr = new TClonesArray("cmsana::TPFCandidate",20000);     assert(fPFCandidateArr);
  fVertexArr      = new TClonesArray("cmsana::TVertex");      assert(fVertexArr);

  //************************************************************************************************
  // Create output file
  //************************************************************************************************
  fOutputFile = new TFile(fOutputName.c_str(),"RECREATE"); 

  //************************************************************************************************
  // Initialize data trees and structs
  //************************************************************************************************
  fEventTree = new TTree("Events","Events");

  fEventTree->Branch("Info",&fEventInfo);
  if(fUseGen) {
    fEventTree->Branch("GenParticle",&fGenParticleArr);  
    fEventTree->Branch("GenJet",&fGenJetArr);  
  }

  fEventTree->Branch("Electron",   &fElectronArr);
  fEventTree->Branch("Muon",       &fMuonArr);
  fEventTree->Branch("PFJet",      &fPFJetArr);
  fEventTree->Branch("Photon",     &fPhotonArr);
  fEventTree->Branch("PFCandidate",&fPFCandidateArr);
  fEventTree->Branch("Vertex",     &fVertexArr);

  cout<<" PATNtupler tree branches defined.. " <<endl; 
  
  

}







double PATNtupler::DeltaPhi(double v1, double v2)
{ // Computes the correctly normalized phi difference
  // v1, v2 = phi of object 1 and 2
  double diff = v1 - v2;
  
  while (diff >acos(-1)) diff -= 2*acos(-1);
  while (diff <= -acos(-1)) diff += 2*acos(-1);
  
  return diff; 
}


double PATNtupler::GetDeltaR(double eta1, double eta2, double phi1, double phi2){ 
  // Computes the DeltaR of two objects from their eta and phi values
  
  return sqrt( (eta1-eta2)*(eta1-eta2) + DeltaPhi(phi1, phi2)*DeltaPhi(phi1, phi2) );
  
}


// ------------ method called once each job just after ending the event loop  ------------
void 
PATNtupler::endJob() 
{


  //fEventTree->Print();
  //fEventTree->Write();

  fOutputFile->cd();
  fOutputFile->Write();
  fOutputFile->Close();
  
  delete fGenParticleArr;
  delete fGenJetArr;
  delete fElectronArr;
  delete fMuonArr;
  delete fPFJetArr;
  delete fPhotonArr;
  delete fPFCandidateArr;
  delete fVertexArr;

}

// ------------ method called when starting to processes a run  ------------
void 
PATNtupler::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
PATNtupler::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
PATNtupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
PATNtupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PATNtupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PATNtupler);

// $Id: CMSNtupler.cc,v 1.14 2013/06/30 13:27:51 sixie Exp $
//
 
#include "CMSAna/CMSNtupler/interface/CMSNtupler.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/plugins/TransientTrackBuilderESProducer.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/EgammaTools/interface/EcalClusterLocal.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"


using namespace std;
using namespace edm;
using namespace reco;


//
// constructors and destructor
//
CMSNtupler::CMSNtupler(const edm::ParameterSet& iConfig)

{
  
  cout << "CMSNtupler initialized... " << endl;
  
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


CMSNtupler::~CMSNtupler()
{
 
  // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
CMSNtupler::analyze(const edm::Event& event, const edm::EventSetup& setup)
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
      if (abs(pgen->pdgId()) == 22 && pgen->pt() > 0.0 && fabs(pgen->eta()) < 5.0)
        fillThisParticle = true;

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
  //Trigger information
  //****************************************************************************************************
//   string hltProcName_;
//   string hltResName_;
//   string hltEvtName_;

//  // get HLT trigger object information to be able to access the tag information
//   Handle<trigger::TriggerEvent> triggerEventHLT;
//   GetProduct("hltTriggerSummaryAOD", triggerEventHLT, event);

//   //set process name from retrieved product in case no process name was specified
//   //(such that we take the hlt from the last run process)
//   if (hltProcName_.empty()) {
//     hltProcName_ = triggerEventHLT.provenance()->processName();
//     cout << "process name: " << triggerEventHLT.provenance()->processName() << "\n";

//     //printf("Extracted most recent process name %s\n",hltProcName_.c_str());
    
//     //add extracted process name to input labels for future access
//     if (hltResName_.find(':')==string::npos)
//       hltResName_ += "::";
//     else 
//       hltResName_ += ":";
//     hltResName_ += hltProcName_;
//     if (hltEvtName_.find(':')==string::npos)
//       hltEvtName_ += "::";
//     else 
//       hltEvtName_ += ":";
//     hltEvtName_ += hltProcName_;
//   }

//   bool hltConfigChanged = false;
//   if (!fHLTConfig.init(event.getRun(), setup, hltProcName_, hltConfigChanged)) {
//     cout << "Error: Cannot access hlt config.\n";
//     return;
//   }
  
//   // loop over hlt paths
//   for(UInt_t i=0;i<fHLTConfig.size();++i) {

//     cout << "HLT " << i << " : " 
//          << fHLTConfig.triggerName(i) << " "
//          << "\n";

//     const vector<string> &mLabels(fHLTConfig.moduleLabels(i));
//     for (UInt_t j=0; j<mLabels.size(); ++j) {
//       //const string& label(mLabels[j]);
//       cout << j << " " << mLabels[j] << "\n";
//     }
//   }
//   cout << "HLT done\n";


//   // get HLT trigger information
//   Handle<TriggerResults> triggerResultsHLT;
//   GetProduct(hltResName_, triggerResultsHLT, event);

//   // get HLT trigger object information
//   Handle<trigger::TriggerEvent> triggerEventHLT;
//   GetProduct(hltEvtName_, triggerEventHLT, event);



  //****************************************************************************************************
  //Pileup information
  //****************************************************************************************************
  
  std::vector<PileupSummaryInfo> inInfos;
  Handle<std::vector< PileupSummaryInfo > >  hPileupInfoProduct;
  event.getByLabel("addPileupInfo", hPileupInfoProduct);
  inInfos = *hPileupInfoProduct.product();

  for (std::vector<PileupSummaryInfo>::const_iterator edmPUInfo = inInfos.begin(); edmPUInfo != inInfos.end(); ++edmPUInfo) {    
//     printf("filling puinfo for bx %i with %i interactions with mean %.2f\n",edmPUInfo->getBunchCrossing(), 
//            edmPUInfo->getPU_NumInteractions(), edmPUInfo->getTrueNumInteractions());
    if (edmPUInfo->getBunchCrossing() == 0) {
      fEventInfo.nPU     = edmPUInfo->getPU_NumInteractions();
      fEventInfo.nPUMean = edmPUInfo->getTrueNumInteractions();
    }
    else if (edmPUInfo->getBunchCrossing() == -1) {
      fEventInfo.nPUMinusOne = edmPUInfo->getPU_NumInteractions();
      fEventInfo.nPUMeanMinusOne = edmPUInfo->getTrueNumInteractions();
    }
    else if (edmPUInfo->getBunchCrossing() == 1) {
      fEventInfo.nPUPlusOne  = edmPUInfo->getPU_NumInteractions();
      fEventInfo.nPUMeanPlusOne = edmPUInfo->getTrueNumInteractions();
    }
  }




  //****************************************************************************************************
  //Transient Track Builder Service
  //****************************************************************************************************
  edm::ESHandle<TransientTrackBuilder> hTransientTrackBuilder;
  setup.get<TransientTrackRecord>().get("TransientTrackBuilder",hTransientTrackBuilder);
  const TransientTrackBuilder *transientTrackBuilder = hTransientTrackBuilder.product();

  //****************************************************************************************************
  //Vertices
  //****************************************************************************************************
  const reco::Vertex *thevtx     = 0;
  const reco::Vertex *thevtxbs   = 0;

  edm::Handle<reco::VertexCollection> hVertex;
  event.getByLabel(fPrimaryVerticesSrcName, hVertex);
  const reco::VertexCollection *pvCol = hVertex.product();
  
  edm::Handle<reco::VertexCollection> hVertexBS;
  event.getByLabel(fPrimaryVerticesBSSrcName, hVertexBS);
  const reco::VertexCollection *pvBSCol = hVertexBS.product();

  UInt_t NGoodPV = 0;
  for (reco::VertexCollection::const_iterator inV = pvCol->begin(); 
       inV != pvCol->end(); ++inV) {

    TClonesArray &rVertexArr = *fVertexArr;
    assert(rVertexArr.GetEntries() < rVertexArr.GetSize());
    const Int_t index = rVertexArr.GetEntries();  
    new(rVertexArr[index]) cmsana::TVertex();
    cmsana::TVertex *pVertex = (cmsana::TVertex*)rVertexArr[index];

    pVertex->x = inV->x();
    pVertex->y = inV->y();
    pVertex->z = inV->z();
    pVertex->isGoodVertex = false;

    if (inV->ndof() >= 4 &&
        fabs(inV->z()) <= 24.0 &&
        inV->position().Rho() <= 2.0
        ) {
      pVertex->isGoodVertex = true;
      if (!thevtx) thevtx = &(*inV);
      NGoodPV++;
    }

    //sumPt
    double sumpt = 0;
    for (reco::Vertex::trackRef_iterator iTrack = inV->tracks_begin(); iTrack!=inV->tracks_end(); ++iTrack) {
      const reco::Track* tmp = dynamic_cast<const reco::Track*>(iTrack->get());
      if (tmp) {
        sumpt += tmp->pt();
      } else {
        cout << "Error: track from vertex cannot be resolved.\n";
      }
    }
    //cout << "Vertex: " << inV->z() << " " << inV->tracksSize() << " " << sumpt << " | " << (inV->z() - genVertexX) << "\n";
  }
  //cout << "\n\n";

  for (reco::VertexCollection::const_iterator inV = pvBSCol->begin(); 
       inV != pvBSCol->end(); ++inV) {

    if (inV->ndof() >= 4 &&
        fabs(inV->z()) <= 24.0 &&
        inV->position().Rho() <= 2.0
        ) {
      if (!thevtxbs) thevtxbs = &(*inV);
      NGoodPV++;
    }
  }

  //****************************************************************************************************
  //Beamspot
  //****************************************************************************************************
  Handle<reco::BeamSpot> hBeamSpotProduct;
  GetProduct("offlineBeamSpot", hBeamSpotProduct, event);
  const reco::BeamSpot *inBeamSpot = hBeamSpotProduct.product();  

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
  Handle<reco::MuonCollection> hMuonProduct;
  GetProduct(fMuonsSrcName, hMuonProduct, event);
  const reco::MuonCollection inMuons = *(hMuonProduct.product());

  Handle<reco::GsfElectronCollection> hElectronProduct;
  GetProduct(fElectronsSrcName, hElectronProduct, event);
  const reco::GsfElectronCollection inElectrons = *(hElectronProduct.product());
 
  Handle<reco::PhotonCollection> hPhotonProduct;
  GetProduct(fPhotonsSrcName, hPhotonProduct, event);
  const reco::PhotonCollection inPhotons = *(hPhotonProduct.product());  

  //****************************************************************************************************
  //PF Candidates
  //****************************************************************************************************
  double customPFMetX = 0;
  double customPFMetY = 0;
  double trackPFMetX = 0;
  double trackPFMetY = 0;

  Handle<reco::PFCandidateCollection> hPfCandProduct;
  GetProduct(fPFCandidatesSrcName, hPfCandProduct, event);  
  const reco::PFCandidateCollection &inPfCands = *(hPfCandProduct.product());
  vector<const reco::PFCandidate*> PFCandidatesSaved;
  PFCandidatesSaved.clear();
  for (reco::PFCandidateCollection::const_iterator iP = inPfCands.begin(); 
       iP != inPfCands.end(); ++iP) {

    //***********************************************************************
    //Figure out if PFCandidate is PFNoPU
    //***********************************************************************
    bool isPFNoPU = true;
    if(iP->particleId() == reco::PFCandidate::h) {
      if(iP->trackRef().isNonnull() && thevtx &&
	 thevtx->trackWeight(iP->trackRef()) > 0) {
        isPFNoPU = true;
      } else { 

        bool vertexFound = false;
        const reco::Vertex *closestVtx = 0;
        double dzmin = 10000;

        // loop over vertices
        for (reco::VertexCollection::const_iterator inV = pvCol->begin(); 
             inV != pvCol->end(); ++inV) {
          if(iP->trackRef().isNonnull() && 
             inV->trackWeight(iP->trackRef()) > 0) {
            vertexFound = true;
            closestVtx = &(*inV);
            break;
          }
          double dz = fabs(iP->vertex().z() - inV->z());
          if(dz < dzmin) {
            closestVtx = &(*inV);
            dzmin = dz;
          }            
        }

        bool fCheckClosestZVertex = true; //we use option 1
        if(fCheckClosestZVertex) {
          // Fallback: if track is not associated with any vertex,
          // associate it with the vertex closest in z
          if(vertexFound || closestVtx != thevtx) {
            isPFNoPU = kFALSE;
          } else {
            isPFNoPU = kTRUE;
          }
        } else {
          if(vertexFound && closestVtx != thevtx) {
            isPFNoPU = kFALSE;
          } else {
            isPFNoPU = kTRUE;
          }
        }
      } //charged hadron & trk stuff
    } else { // neutrals 
      //
      isPFNoPU = kTRUE;
    }

    //******************************************************
    //Compute Custom MET
    //******************************************************
    customPFMetX -= iP->px();
    customPFMetY -= iP->py();
    if (iP->charge() != 0 && isPFNoPU) {
      trackPFMetX -= iP->px();
      trackPFMetY -= iP->py();
    }

    //******************************************************
    //Only Fill PF candidates that are within deltaR 0.5
    //to muons, electrons, photons, (add taus later)
    //******************************************************
    bool fillThisPFCandidate = false;
    for (reco::MuonCollection::const_iterator iM = inMuons.begin(); iM != inMuons.end(); ++iM) {
      if (deltaR(*iM,*iP) < 0.5) fillThisPFCandidate = true;
    }
    for (reco::GsfElectronCollection::const_iterator iE = inElectrons.begin(); 
         iE != inElectrons.end(); ++iE) {
      if (deltaR(*iE,*iP) < 0.5) fillThisPFCandidate = true;
    }
    for (reco::PhotonCollection::const_iterator iPh = inPhotons.begin(); 
         iPh != inPhotons.end(); ++iPh) {
      if (!(iPh->pt() > fPhotonPtMin)) continue;
      if (deltaR(*iPh,*iP) < 0.5) fillThisPFCandidate = true;
    }
    if (!fillThisPFCandidate) continue;
    //******************************************************

    PFCandidatesSaved.push_back(&(*iP));

    TClonesArray &rPFCandidateArr = *fPFCandidateArr;
    if (!(rPFCandidateArr.GetEntries() < rPFCandidateArr.GetSize())) {
      cout << "Error: PFCandidateArray has " << rPFCandidateArr.GetEntries() << " entries and is over the limit.\n";      
    }

    assert(rPFCandidateArr.GetEntries() < rPFCandidateArr.GetSize());
    const Int_t index = rPFCandidateArr.GetEntries();  
    new(rPFCandidateArr[index]) cmsana::TPFCandidate();
    cmsana::TPFCandidate *pPFCandidate = (cmsana::TPFCandidate*)rPFCandidateArr[index]; 

    pPFCandidate->pt	  = iP->pt();
    pPFCandidate->eta  	  = iP->eta();
    pPFCandidate->phi  	  = iP->phi();
    pPFCandidate->e  	  = iP->energy();  
    pPFCandidate->q        = iP->charge();

    if (thevtx && iP->gsfTrackRef().isNonnull()) {
      pPFCandidate->dz = iP->gsfTrackRef()->dz(thevtx->position());
    } else if (thevtx && iP->trackRef().isNonnull()) {
      pPFCandidate->dz = iP->trackRef()->dz(thevtx->position());
    } else { 
      pPFCandidate->dz       = 0.0;
    }

    pPFCandidate->pfType    = iP->particleId();
    pPFCandidate->IsPFNoPU = isPFNoPU;
  }


  //****************************************************************************************************
  //Muons
  //****************************************************************************************************

  for (reco::MuonCollection::const_iterator iM = inMuons.begin(); iM != inMuons.end(); ++iM) {
 
    TClonesArray &rMuonArr = *fMuonArr;
    assert(rMuonArr.GetEntries() < rMuonArr.GetSize());
    const Int_t index = rMuonArr.GetEntries();  
    new(rMuonArr[index]) cmsana::TMuon();
    cmsana::TMuon *pMuon = (cmsana::TMuon*)rMuonArr[index];


    pMuon->typeBits = 0;
    if(iM->isGlobalMuon())     { pMuon->typeBits |= kGlobal; }
    if(iM->isTrackerMuon())    { pMuon->typeBits |= kTracker; }
    if(iM->isStandAloneMuon()) { pMuon->typeBits |= kStandalone; }



    if (iM->track().isNonnull()) {
      pMuon->pt       = iM->track()->pt();
      pMuon->eta      = iM->track()->eta();
      pMuon->phi      = iM->track()->phi();
      pMuon->pterr    = iM->track()->ptError();
      pMuon->q        = iM->track()->charge();
      pMuon->nTkHits  = iM->track()->numberOfValidHits();
      pMuon->nPixHits = iM->track()->hitPattern().numberOfValidPixelHits();
      pMuon->tkNchi2  = iM->track()->normalizedChi2();
    } else if (iM->standAloneMuon().isNonnull()) {
      pMuon->pt       = iM->standAloneMuon()->pt();
      pMuon->eta      = iM->standAloneMuon()->eta();
      pMuon->phi      = iM->standAloneMuon()->phi();
      pMuon->pterr    = iM->standAloneMuon()->ptError();
      pMuon->q        = iM->standAloneMuon()->charge();
      pMuon->nTkHits  = iM->standAloneMuon()->numberOfValidHits();
      pMuon->nPixHits = iM->standAloneMuon()->hitPattern().numberOfValidPixelHits();
      pMuon->tkNchi2  = iM->standAloneMuon()->normalizedChi2();
    } else {
      cout << "Warning: muon has no tracker track and no standalone track.\n";
    }

    //Muon Normalized Chi^2
    if (iM->combinedMuon().isNonnull()) {
      pMuon->muNchi2 = iM->combinedMuon()->normalizedChi2();
    } else if (iM->track().isNonnull()) {
      pMuon->muNchi2 = iM->track()->normalizedChi2();
    } else if (iM->standAloneMuon().isNonnull()) {
      pMuon->muNchi2 = iM->standAloneMuon()->normalizedChi2();
    }

    pMuon->trkIso03  = iM->isolationR03().sumPt;
    pMuon->emIso03   = iM->isolationR03().emEt;
    pMuon->hadIso03  = iM->isolationR03().hadEt;
    pMuon->hoIso03   = iM->isolationR03().hoEt;

//     pMuon->trkIso05  = iM->isolationR05().sumPt;
//     pMuon->emIso05   = iM->isolationR05().emEt;
//     pMuon->hadIso05  = iM->isolationR05().hadEt;
//     pMuon->hoIso05   = iM->isolationR05().hoEt;

    pMuon->pfIso04ChargedHadron              = iM->pfIsolationR04().sumChargedHadronPt;
    pMuon->pfIso04ChargedParticle            = iM->pfIsolationR04().sumChargedParticlePt;
    pMuon->pfIso04NeutralHadron              = iM->pfIsolationR04().sumNeutralHadronEt;
    pMuon->pfIso04Photon                     = iM->pfIsolationR04().sumPhotonEt;
    pMuon->pfIso04NeutralHadronHighThreshold = iM->pfIsolationR04().sumNeutralHadronEtHighThreshold;
    pMuon->pfIso04PhotonHighThreshold        = iM->pfIsolationR04().sumPhotonEtHighThreshold;    
    pMuon->pfIso04PU                         = iM->pfIsolationR04().sumPUPt;

//     iM->pfIsolationR03().sumChargedHadronPt;
//     iM->pfIsolationR03().sumChargedParticlePt;
//     iM->pfIsolationR03().sumNeutralHadronEt;
//     iM->pfIsolationR03().sumPhotonEt;
//     iM->pfIsolationR03().sumNeutralHadronEtHighThreshold;
//     iM->pfIsolationR03().sumPhotonEtHighThreshold;    
//     iM->pfIsolationR03().sumPUPt;

    if (iM->combinedMuon().isNonnull()) {
      pMuon->nValidHits = iM->globalTrack()->hitPattern().numberOfValidMuonHits();
    }

    pMuon->nMatch    = iM->numberOfMatchedStations();

    if (thevtx && iM->muonBestTrack().isNonnull()) {      
      pMuon->d0      = iM->muonBestTrack()->dxy(thevtx->position());
      pMuon->dz      = iM->muonBestTrack()->dz(thevtx->position());
    }
    
    if (thevtx && iM->track().isNonnull()) {
      pMuon->trkLayers = iM->track()->hitPattern().trackerLayersWithMeasurement();

      const reco::TransientTrack &tt = transientTrackBuilder->build(iM->track());

      //preserve sign of transverse impact parameter (cross-product definition from track, not lifetime-signing)
      const double thesign   = ( (-iM->track()->dxy(thevtx->position()))   >=0 ) ? 1. : -1.;
      double thesignbs;
      if (thevtxbs) thesignbs = ( (-iM->track()->dxy(thevtxbs->position())) >=0 ) ? 1. : -1.;
      else thesignbs = thesign;

      const std::pair<bool,Measurement1D> &ip3dpv =  IPTools::absoluteImpactParameter3D(tt,*thevtx);
      if (ip3dpv.first) {
        pMuon->ip3d    =  thesign*ip3dpv.second.value();
        pMuon->ip3dSig =  thesign*ip3dpv.second.value() / ip3dpv.second.error();
      } else {
        pMuon->ip3d      = -99.0;
        pMuon->ip3dSig   = -99.0;
      }

      if (thevtxbs) {
	const std::pair<bool,Measurement1D> &ip3dpvbs =  IPTools::absoluteImpactParameter3D(tt,*thevtxbs);
	if (ip3dpvbs.first) {
	  pMuon->ip3dBS    =  thesignbs*ip3dpvbs.second.value();
	  pMuon->ip3dSigBS =  thesignbs*ip3dpvbs.second.value() / ip3dpvbs.second.error();
	} else {
	  pMuon->ip3dBS    = -99.0;
	  pMuon->ip3dSigBS = -99.0;
	}
      }

    }

    UInt_t muonQuality = 0;        
    muonQuality |= kAll;
    if (muon::isGoodMuon(*iM,muon::AllGlobalMuons))
      muonQuality |= kAllGlobalMuons;
    if (muon::isGoodMuon(*iM,muon::AllStandAloneMuons))
      muonQuality |= kAllStandAloneMuons;
    if (muon::isGoodMuon(*iM,muon::AllTrackerMuons))
      muonQuality |= kAllTrackerMuons;
    if (muon::isGoodMuon(*iM,muon::TrackerMuonArbitrated))
      muonQuality |= kTrackerMuonArbitrated;
    if (muon::isGoodMuon(*iM,muon::AllArbitrated))
      muonQuality |= kAllArbitrated;
    if (muon::isGoodMuon(*iM,muon::GlobalMuonPromptTight))
      muonQuality |= kGlobalMuonPromptTight;
    if (muon::isGoodMuon(*iM,muon::TMLastStationLoose))
      muonQuality |= kTMLastStationLoose;
    if (muon::isGoodMuon(*iM,muon::TMLastStationTight))
      muonQuality |= kTMLastStationTight;
    if (muon::isGoodMuon(*iM,muon::TM2DCompatibilityLoose))
      muonQuality |= kTM2DCompatibilityLoose;
    if (muon::isGoodMuon(*iM,muon::TM2DCompatibilityTight))
      muonQuality |= kTM2DCompatibilityTight;
    if (muon::isGoodMuon(*iM,muon::TMOneStationLoose))
      muonQuality |= kTMOneStationLoose;
    if (muon::isGoodMuon(*iM,muon::TMOneStationTight))
      muonQuality |= kTMOneStationTight;
    if (muon::isGoodMuon(*iM,muon::TMLastStationOptimizedLowPtLoose))
      muonQuality |= kTMLastStationOptimizedLowPtLoose;
    if (muon::isGoodMuon(*iM,muon::TMLastStationOptimizedLowPtTight))
      muonQuality |= kTMLastStationOptimizedLowPtTight;
    if (muon::isGoodMuon(*iM,muon::GMTkChiCompatibility))
      muonQuality |= kGMTkChiCompatibility;
    if (muon::isGoodMuon(*iM,muon::GMStaChiCompatibility))
      muonQuality |= kGMStaChiCompatibility;
    if (muon::isGoodMuon(*iM,muon::GMTkKinkTight))
      muonQuality |= kGMTkKinkTight;
    if (muon::isGoodMuon(*iM,muon::TMLastStationAngLoose))
      muonQuality |= kTMLastStationAngLoose;
    if (muon::isGoodMuon(*iM,muon::TMLastStationAngTight))
      muonQuality |= kTMLastStationAngTight;
    if (muon::isGoodMuon(*iM,muon::TMOneStationAngLoose))
      muonQuality |= kTMOneStationAngLoose;
    if (muon::isGoodMuon(*iM,muon::TMOneStationAngTight))
      muonQuality |= kTMOneStationAngTight;
    if (muon::isGoodMuon(*iM,muon::TMLastStationOptimizedBarrelLowPtLoose))
      muonQuality |= kTMLastStationOptimizedBarrelLowPtLoose;
    if (muon::isGoodMuon(*iM,muon::TMLastStationOptimizedBarrelLowPtTight))
      muonQuality |= kTMLastStationOptimizedBarrelLowPtTight;
    
//     pMuon->hltMatchBits = MatchHLTMuon(muTrk->Pt(),muTrk->Eta(),muTrk->Phi(),GetEventHeader()->RunNum(), GetEventHeader()->EvtNum() );    
//     pMuon->isMCReal = isMCMatched;

    //Save a list of indices of PFCandidates that match the supercluster used for the muon
    vector<UInt_t> matchedPFCandidateIndex;
    for (uint p = 0; p < PFCandidatesSaved.size() ; ++p) {
      if (iM->track().isNonnull() && iM->track() == PFCandidatesSaved[p]->trackRef() ) {
        matchedPFCandidateIndex.push_back(p);
      }
    }
    pMuon->NMatchedPFCandidates = matchedPFCandidateIndex.size();
    for (uint p=0; p<matchedPFCandidateIndex.size(); ++p) {
      if (p >= 10) break; //can store at most 10
      pMuon->MatchedPFCandidateIndex[p] = matchedPFCandidateIndex[p];
    }
  
  } //loop over muons
  


  EcalClusterLazyTools *lazyTools = 0;
  EcalClusterLocal *local = 0;  
  if (fFillEGRegressionVars) {
    lazyTools = new EcalClusterLazyTools(event, setup, edm::InputTag("reducedEcalRecHitsEB"), 
                                         edm::InputTag("reducedEcalRecHitsEE"));
    local = new EcalClusterLocal;
  }

  //****************************************************************************************************
  //Electrons
  //****************************************************************************************************
  edm::Handle<reco::ConversionCollection> hConversions;
  event.getByLabel("allConversions", hConversions);

  for (reco::GsfElectronCollection::const_iterator iE = inElectrons.begin(); 
       iE != inElectrons.end(); ++iE) {

    TClonesArray &rElectronArr = *fElectronArr;
    assert(rElectronArr.GetEntries() < rElectronArr.GetSize());
    const Int_t index = rElectronArr.GetEntries();  
    new(rElectronArr[index]) cmsana::TElectron();
    cmsana::TElectron *pElectron = (cmsana::TElectron*)rElectronArr[index];

    pElectron->pt              = iE->pt();
    pElectron->eta             = iE->eta();
    pElectron->phi             = iE->phi();

    if(iE->gsfTrack().isNonnull()) {
      pElectron->p = iE->gsfTrack()->p();
    } else if (iE->closestCtfTrackRef().isNonnull()) {
      pElectron->p = iE->closestCtfTrackRef()->p();
    }

    if (iE->superCluster().isNonnull()) {
      pElectron->scEt        = iE->superCluster()->energy() / cosh( iE->superCluster()->eta() ) ;
      pElectron->scEta           = iE->superCluster()->eta();
      pElectron->scPhi           = iE->superCluster()->phi();
    }

    pElectron->q               = iE->charge();


    pElectron->isEcalDriven    = iE->ecalDrivenSeed();
    pElectron->isTrackerDriven = iE->trackerDrivenSeed();
    pElectron->isEB            = iE->isEB();
    pElectron->isEE            = iE->isEE();
    pElectron->Classification  = iE->classification();

//     pElectron->isMCReal        = isMCMatched;
//     pElectron->hltMatchBits    = MatchHLT(ele->SCluster()->Eta(),ele->SCluster()->Phi(), GetEventHeader()->RunNum(), GetEventHeader()->EvtNum());  
//     pElectron->l1TriggerMatchBits = MatchL1(ele->SCluster()->Eta(),ele->SCluster()->Phi());
    
    pElectron->TrackMomentumError = iE->trackMomentumError();
    pElectron->nBrem           = iE->basicClustersSize() - 1;
    pElectron->fBrem           = iE->fbrem();
    pElectron->EOverP          = iE->eSuperClusterOverP();
    pElectron->pIn             = iE->trackMomentumAtVtx().R();
    pElectron->ESeedClusterOverPIn  = iE->superCluster()->seed()->energy() / iE->trackMomentumAtVtx().R();
    pElectron->ESeedClusterOverPout = iE->eSeedClusterOverPout();
    pElectron->EEleClusterOverPout  = iE->eEleClusterOverPout();
    pElectron->EcalEnergy      = iE->correctedEcalEnergy();
    pElectron->EcalEnergyError = iE->correctedEcalEnergyError();

    pElectron->deltaEtaIn      = iE->deltaEtaSuperClusterTrackAtVtx();
    pElectron->deltaPhiIn      = iE->deltaPhiSuperClusterTrackAtVtx();
    pElectron->dEtaCalo        = iE->deltaEtaSeedClusterTrackAtCalo();
    pElectron->dPhiCalo        = iE->deltaPhiSeedClusterTrackAtCalo();
    pElectron->sigiEtaiEta     = iE->sigmaIetaIeta();

    
    if (fFillEGRegressionVars) {
      pElectron->sigiPhiiPhi     = sqrt(lazyTools->covariances(*(iE->superCluster()->seed()))[2]);
      if (isnan(pElectron->sigiPhiiPhi))  pElectron->sigiPhiiPhi = 0.0;
      
      
      if (pElectron->sigiEtaiEta*pElectron->sigiPhiiPhi > 0) {
        pElectron->sigiEtaiPhi = lazyTools->covariances(*(iE->superCluster()->seed()))[1]/(pElectron->sigiEtaiEta*pElectron->sigiPhiiPhi);
      } else if (lazyTools->covariances(*(iE->superCluster()->seed()))[1]>0) {
        pElectron->sigiEtaiPhi = 1.0; 
      } else {
        pElectron->sigiEtaiPhi = -1.0; 
      }
      pElectron->R9 = lazyTools->e3x3(*(iE->superCluster()->seed())) / iE->superCluster()->rawEnergy();
      if (lazyTools->e5x5(*(iE->superCluster()->seed())) > 0) {
        pElectron->SeedE1x5OverE5x5 = lazyTools->e1x5(*(iE->superCluster()->seed())) / lazyTools->e5x5(*(iE->superCluster()->seed()));
      } 
    }
    
    pElectron->SCEtaWidth = iE->superCluster()->etaWidth();
    pElectron->SCPhiWidth = iE->superCluster()->phiWidth();
    pElectron->PreShowerOverRaw   = iE->superCluster()->preshowerEnergy() / iE->superCluster()->rawEnergy();
    pElectron->HoverE             = iE->hcalOverEcal();
    pElectron->HoverESingleTower  = iE->hcalOverEcalBc();

    if(iE->gsfTrack().isNonnull()) {
      pElectron->GsfTrackChi2OverNdof = iE->gsfTrack()->normalizedChi2();
    }
    if (iE->closestCtfTrackRef().isNonnull()) {
      pElectron->KFTrackChi2OverNdof           = iE->closestCtfTrackRef()->normalizedChi2();
      pElectron->KFTrackNHits                  = iE->closestCtfTrackRef()->numberOfValidHits();
      pElectron->KFTrackNLayersWithMeasurement = iE->closestCtfTrackRef()->hitPattern().trackerLayersWithMeasurement();
    }


    
    if(thevtx && iE->gsfTrack().isNonnull()) {

      const reco::TransientTrack &tt = transientTrackBuilder->build(iE->gsfTrack());

      pElectron->d0              = iE->gsfTrack()->dxy(thevtx->position());
      pElectron->dz              = iE->gsfTrack()->dz(thevtx->position());

      const double gsfsign   = ( (-iE->gsfTrack()->dxy(thevtx->position()))   >=0 ) ? 1. : -1.;
      double gsfsignbs;
      if (thevtxbs) gsfsignbs = ( (-iE->gsfTrack()->dxy(thevtxbs->position())) >=0 ) ? 1. : -1.;
      else gsfsignbs = gsfsign;
      
      const std::pair<bool,Measurement1D> &ip3dpv =  IPTools::absoluteImpactParameter3D(tt,*thevtx);
      if (ip3dpv.first) {
        pElectron->ip3d = gsfsign*ip3dpv.second.value();
        pElectron->ip3dSig = gsfsign*ip3dpv.second.value() / ip3dpv.second.error();
      } 
      
      if (thevtxbs) {
	const std::pair<bool,Measurement1D> &ip3dpvbs =  IPTools::absoluteImpactParameter3D(tt,*thevtxbs);
	if (ip3dpvbs.first) {
	  pElectron->ip3dBS = gsfsignbs*ip3dpvbs.second.value();
	  pElectron->ip3dSigBS = gsfsignbs*ip3dpvbs.second.value() / ip3dpvbs.second.error();
	}
      }
      pElectron->nExpHitsInner   = iE->gsfTrack()->trackerExpectedHitsInner().numberOfHits();     
    }

    //Do this later
    if (inBeamSpot) {
      pElectron->isConv = ConversionTools::hasMatchedConversion(*iE,hConversions,inBeamSpot->position());
    }


    pElectron->trkIso03        = iE->dr03TkSumPt();
    pElectron->emIso03         = iE->dr03EcalRecHitSumEt();
    pElectron->hadIso03        = iE->dr03HcalDepth1TowerSumEt() + iE->dr03HcalDepth2TowerSumEt();
    pElectron->trkIso04        = iE->dr04TkSumPt();
    pElectron->emIso04         = iE->dr04EcalRecHitSumEt();
    pElectron->hadIso04        = iE->dr04HcalDepth1TowerSumEt() + iE->dr04HcalDepth2TowerSumEt();
    //PFIsolation needs to be implemented later
  

    pElectron->SCRawEnergy = iE->superCluster()->rawEnergy();
    pElectron->E5x5 = iE->e5x5();
    pElectron->EtaSeed = iE->superCluster()->seed()->eta();
    pElectron->PhiSeed = iE->superCluster()->seed()->phi();
    pElectron->ESeed   = iE->superCluster()->seed()->energy();

    if (fFillEGRegressionVars) {
      pElectron->E3x3Seed = lazyTools->e3x3(*(iE->superCluster()->seed()));
      pElectron->E5x5Seed = lazyTools->e5x5(*(iE->superCluster()->seed()));
      pElectron->EMaxSeed = lazyTools->eMax(*(iE->superCluster()->seed()));
      pElectron->E2ndSeed = lazyTools->e2nd(*(iE->superCluster()->seed()));
      pElectron->ETopSeed= lazyTools->e2x5Top(*(iE->superCluster()->seed()));
      pElectron->EBottomSeed = lazyTools->eBottom(*(iE->superCluster()->seed()));
      pElectron->ELeftSeed = lazyTools->eLeft(*(iE->superCluster()->seed()));
      pElectron->ERightSeed = lazyTools->eRight(*(iE->superCluster()->seed()));
      pElectron->E2x5MaxSeed = lazyTools->e2x5Max(*(iE->superCluster()->seed()));
      pElectron->E2x5TopSeed = lazyTools->e2x5Top(*(iE->superCluster()->seed()));
      pElectron->E2x5BottomSeed = lazyTools->e2x5Bottom(*(iE->superCluster()->seed()));
      pElectron->E2x5LeftSeed = lazyTools->e2x5Left(*(iE->superCluster()->seed()));
      pElectron->E2x5RightSeed = lazyTools->e2x5Right(*(iE->superCluster()->seed()));


      //Turn this off for now because it depends on detector geometry
      //local coordinates
      if (iE->superCluster()->seed()->hitsAndFractions().at(0).first.subdetId()==EcalBarrel) {
        float etacry, phicry, thetatilt, phitilt;
        int ieta, iphi;
        local->localCoordsEB(*(iE->superCluster()->seed()),setup,etacry,phicry,ieta,iphi,thetatilt,phitilt);

        pElectron->IEtaSeed = ieta;
        pElectron->IPhiSeed = iphi;
        pElectron->EtaCrySeed = etacry;
        pElectron->PhiCrySeed = phicry;        
      }
      else {
        pElectron->IEtaSeed = -999;
        pElectron->IPhiSeed = -999;
        pElectron->EtaCrySeed = -999;
        pElectron->PhiCrySeed = -999;
      }
    }

    //Save a list of indices of PFCandidates that match the supercluster used for the ele
    vector<UInt_t> scMatchedPFCandidateIndex;
    vector<UInt_t> gsfTrackMatchedPFCandidateIndex;
    for (uint p = 0; p < PFCandidatesSaved.size() ; ++p) {
      if (iE->superCluster().isNonnull() && iE->superCluster() == PFCandidatesSaved[p]->superClusterRef()) {
        scMatchedPFCandidateIndex.push_back(p);
      }
      if (iE->gsfTrack().isNonnull() && iE->gsfTrack() == PFCandidatesSaved[p]->gsfTrackRef()) {
        gsfTrackMatchedPFCandidateIndex.push_back(p);
      }
    }
    pElectron->NSCMatchedPFCandidates = scMatchedPFCandidateIndex.size();
    for (uint p=0; p<scMatchedPFCandidateIndex.size(); ++p) {
      if (p >= 10) break; //can store at most 10
      pElectron->SCMatchedPFCandidateIndex[p] = scMatchedPFCandidateIndex[p];
    }
    pElectron->NGsfTrackMatchedPFCandidates = gsfTrackMatchedPFCandidateIndex.size();
    for (uint p=0; p<gsfTrackMatchedPFCandidateIndex.size(); ++p) {
      if (p >= 10) break; //can store at most 10
      pElectron->GsfTrackMatchedPFCandidateIndex[p] = gsfTrackMatchedPFCandidateIndex[p];
    }


  } //loop over electrons


  //****************************************************************************************************
  //Photons
  //****************************************************************************************************
  for (reco::PhotonCollection::const_iterator iP = inPhotons.begin(); 
       iP != inPhotons.end(); ++iP) {

    //min pt cuts for saving photons
    if (!(iP->pt() > fPhotonPtMin)) continue;
    
    TClonesArray &rPhotonArr = *fPhotonArr;
    assert(rPhotonArr.GetEntries() < rPhotonArr.GetSize());
    const Int_t index = rPhotonArr.GetEntries();  
    new(rPhotonArr[index]) cmsana::TPhoton();
    cmsana::TPhoton *pPhoton = (cmsana::TPhoton*)rPhotonArr[index];
    
    pPhoton->pt		       = iP->pt();
    pPhoton->eta  	       = iP->eta();
    pPhoton->phi  	       = iP->phi();
    if (iP->superCluster().isNonnull()) {
      pPhoton->scEt	       = iP->superCluster()->energy() / cosh( iP->superCluster()->eta() ) ;
      pPhoton->scEta  	       = iP->superCluster()->eta();
      pPhoton->scPhi  	       = iP->superCluster()->phi();
    }
    pPhoton->trkIso03Hollow    = iP->trkSumPtHollowConeDR03();
    pPhoton->trkIso03Solid     = iP->trkSumPtSolidConeDR03();
    pPhoton->emIso03           = iP->ecalRecHitSumEtConeDR03();
    pPhoton->hadIso03	       = iP->hcalTowerSumEtConeDR03();
    pPhoton->HoverE	       = iP->hadronicOverEm();
    pPhoton->HoverESingleTower = iP->hadTowOverEm();
    pPhoton->R9		       = iP->r9();
    pPhoton->sigiEtaiEta       = iP->sigmaEtaEta();
    pPhoton->hasPixelSeed      = iP->hasPixelSeed();
    pPhoton->isConversion      = iP->hasConversionTracks();
    if (inBeamSpot) {
      pPhoton->passEleVeto
	= !ConversionTools::hasMatchedPromptElectron(iP->superCluster(),hElectronProduct, 
						     hConversions, inBeamSpot->position());
    }
    
    pPhoton->SCRawEnergy       = iP->superCluster()->rawEnergy();
    pPhoton->PreShowerOverRaw  = iP->superCluster()->preshowerEnergy() / iP->superCluster()->rawEnergy();
    pPhoton->SCEtaWidth        = iP->superCluster()->etaWidth();
    pPhoton->SCPhiWidth        = iP->superCluster()->phiWidth();
    pPhoton->NClusters         = iP->superCluster()->clustersSize();
    pPhoton->E5x5              = iP->e5x5();
    pPhoton->ESeed             = iP->superCluster()->seed()->energy();
    pPhoton->EtaSeed           = iP->superCluster()->seed()->eta();
    pPhoton->PhiSeed           = iP->superCluster()->seed()->phi();

    pPhoton->sigiPhiiPhi       = sqrt(lazyTools->covariances(*(iP->superCluster()->seed()))[2]);
    if (isnan(pPhoton->sigiPhiiPhi))  pPhoton->sigiPhiiPhi = 0.0;

    if (pPhoton->sigiEtaiEta*pPhoton->sigiPhiiPhi > 0) {
      pPhoton->sigiEtaiPhi     = lazyTools->covariances(*(iP->superCluster()->seed()))[1]/(pPhoton->sigiEtaiEta*pPhoton->sigiPhiiPhi);
    } else if (lazyTools->covariances(*(iP->superCluster()->seed()))[1]>0) {
      pPhoton->sigiEtaiPhi     = 1.0; 
    } else {
      pPhoton->sigiEtaiPhi     = -1.0; 
    }

    if (fFillEGRegressionVars) {
      pPhoton->E2x2Seed          = lazyTools->e2x2(*(iP->superCluster()->seed()));
      pPhoton->E3x3Seed          = lazyTools->e3x3(*(iP->superCluster()->seed()));
      pPhoton->E5x5Seed          = lazyTools->e5x5(*(iP->superCluster()->seed()));
      pPhoton->EMaxSeed          = lazyTools->eMax(*(iP->superCluster()->seed()));
      pPhoton->E2ndSeed          = lazyTools->e2nd(*(iP->superCluster()->seed()));
      pPhoton->ETopSeed          = lazyTools->e2x5Top(*(iP->superCluster()->seed()));
      pPhoton->EBottomSeed       = lazyTools->eBottom(*(iP->superCluster()->seed()));
      pPhoton->ELeftSeed         = lazyTools->eLeft(*(iP->superCluster()->seed()));
      pPhoton->ERightSeed        = lazyTools->eRight(*(iP->superCluster()->seed()));
      pPhoton->E2x5MaxSeed       = lazyTools->e2x5Max(*(iP->superCluster()->seed()));
      pPhoton->E2x5TopSeed       = lazyTools->e2x5Top(*(iP->superCluster()->seed()));
      pPhoton->E2x5BottomSeed    = lazyTools->e2x5Bottom(*(iP->superCluster()->seed()));
      pPhoton->E2x5LeftSeed      = lazyTools->e2x5Left(*(iP->superCluster()->seed()));
      pPhoton->E2x5RightSeed     = lazyTools->e2x5Right(*(iP->superCluster()->seed()));

      if (iP->superCluster()->seed()->hitsAndFractions().at(0).first.subdetId()==EcalBarrel) {
        float etacry, phicry, thetatilt, phitilt;
        int ieta, iphi;
        local->localCoordsEB(*(iP->superCluster()->seed()),setup,etacry,phicry,ieta,iphi,thetatilt,phitilt);
      
        pPhoton->IEtaSeed        = ieta;
        pPhoton->IPhiSeed        = iphi;
        pPhoton->EtaCrySeed      = etacry;
        pPhoton->PhiCrySeed      = phicry;
      } else {
        pPhoton->IEtaSeed = -999;
        pPhoton->IPhiSeed = -999;
        pPhoton->EtaCrySeed = -999;
        pPhoton->PhiCrySeed = -999;
      }
    }

    //Save a list of indices of PFCandidates that match the supercluster used for the photon
    vector<UInt_t> matchedPFCandidateIndex;
    for (uint p = 0; p < PFCandidatesSaved.size() ; ++p) {
      if (iP->superCluster().isNonnull() && iP->superCluster() == PFCandidatesSaved[p]->superClusterRef()) {
        matchedPFCandidateIndex.push_back(p);
      }
    }
    pPhoton->NMatchedPFCandidates = matchedPFCandidateIndex.size();
    for (uint p=0; p<matchedPFCandidateIndex.size(); ++p) {
      if (p >= 10) break; //can store at most 10
      pPhoton->MatchedPFCandidateIndex[p] = matchedPFCandidateIndex[p];
    }

  } //loop over photons




  //****************************************************************************************************
  //PFJets
  //****************************************************************************************************
  Handle<reco::JetTagCollection> hJetProbabilityBJetTags;
  Handle<reco::JetTagCollection> hJetBProbabilityBJetTags;
  Handle<reco::JetTagCollection> hSimpleSecondaryVertexBJetTags;
  Handle<reco::JetTagCollection> hSimpleSecondaryVertexHighEffBJetTags;
  Handle<reco::JetTagCollection> hSimpleSecondaryVertexHighPurBJetTags;
  Handle<reco::JetTagCollection> hCombinedSecondaryVertexBJetTags;
  Handle<reco::JetTagCollection> hCombinedSecondaryVertexMVABJetTags;
  Handle<reco::JetTagCollection> hTrackCountingHighEffBJetTags;
  Handle<reco::JetTagCollection> hTrackCountingHighPurBJetTags;
  Handle<reco::JetTagCollection> hSoftMuonBJetTags;
  Handle<reco::JetTagCollection> hSoftMuonByIP3dBJetTags;
  Handle<reco::JetTagCollection> hSoftMuonByPtBJetTags;
  Handle<reco::JetTagCollection> hSoftElectronByIP3dBJetTags;
  Handle<reco::JetTagCollection> hSoftElectronByPtBJetTags;
  Handle<reco::JetTagCollection> hGhostTrackBJetTags;
  GetProduct("newJetProbabilityBJetTags", hJetProbabilityBJetTags, event);    
  GetProduct("newJetBProbabilityBJetTags", hJetBProbabilityBJetTags, event);        
  event.getByLabel("newSimpleSecondaryVertexHighEffBJetTags",hSimpleSecondaryVertexHighEffBJetTags);
  event.getByLabel("newSimpleSecondaryVertexHighPurBJetTags",hSimpleSecondaryVertexHighPurBJetTags);
  GetProduct("newCombinedSecondaryVertexBJetTags", hCombinedSecondaryVertexBJetTags, event);    
  GetProduct("newCombinedSecondaryVertexMVABJetTags", hCombinedSecondaryVertexMVABJetTags, event);
  GetProduct("newTrackCountingHighEffBJetTags", hTrackCountingHighEffBJetTags, event);    
  GetProduct("newTrackCountingHighPurBJetTags", hTrackCountingHighPurBJetTags, event);    
  GetProduct("newSoftMuonBJetTags", hSoftMuonBJetTags, event);    
  GetProduct("newSoftMuonByIP3dBJetTags", hSoftMuonByIP3dBJetTags, event);
  GetProduct("newSoftMuonByPtBJetTags", hSoftMuonByPtBJetTags, event);   
  GetProduct("newSoftElectronByIP3dBJetTags", hSoftElectronByIP3dBJetTags, event);
  GetProduct("newSoftElectronByPtBJetTags", hSoftElectronByPtBJetTags, event);    
  
  //handle for jet flavour association
  edm::Handle<reco::JetFlavourMatchingCollection> PFJetFlavMatch;
  if (fUseGen) {
    event.getByLabel ("myAK5PFJetFlavourAssociation", PFJetFlavMatch);
  }


  Handle<double> rhoForJEC;
  event.getByLabel(edm::InputTag("kt6PFJets","rho"),rhoForJEC);

  //Get Jet Collection
  Handle<reco::PFJetCollection> hPFJetProduct;
  GetProduct(fPFJetsSrcName, hPFJetProduct, event);
  const reco::PFJetCollection inJets = *(hPFJetProduct.product());  

  // loop through all jets
  for (reco::PFJetCollection::const_iterator inJet = inJets.begin(); 
       inJet != inJets.end(); ++inJet) {

    reco::PFJetRef jetRef(hPFJetProduct, inJet - inJets.begin());    
    reco::JetBaseRef jetBaseRef(jetRef);

    if (!(inJet->pt() > fJetPtMin)) continue;
    
    TClonesArray &rPFJetArr = *fPFJetArr;
    assert(rPFJetArr.GetEntries() < rPFJetArr.GetSize());
    const Int_t index = rPFJetArr.GetEntries();  
    new(rPFJetArr[index]) cmsana::TJet();
    cmsana::TJet *pPFJet = (cmsana::TJet*)rPFJetArr[index]; 
   
    //Jet Energy corrections

    pPFJet->rawPt        = inJet->pt();
    pPFJet->eta          = inJet->eta();
    pPFJet->phi          = inJet->phi();
    pPFJet->mass         = inJet->mass();

    double L1Scale = (inJet->pt() - (*rhoForJEC)*inJet->jetArea())/inJet->pt();
    L1Scale = (L1Scale>0) ? L1Scale : 0.0;
    double L2Scale = 1.0; //no corr for now
    double L3Scale = 1.0; //no corr for now

//     L2Scale = correctorL2->correction(inJet->p4());
//     L3Scale = correctorL3->correction(inJet->p4()*L2Scale);
    pPFJet->L1JECScale = L1Scale;
    pPFJet->L2JECScale = L2Scale;
    pPFJet->L3JECScale = L3Scale;
    pPFJet->pt         = inJet->pt()*L1Scale*L2Scale*L3Scale;    

    //pPFJet->hltMatchBits = MatchHLT(jet->Eta(),jet->Phi());  

    pPFJet->NConstituents = inJet->numberOfDaughters();
    pPFJet->ChargedEMFraction = inJet->chargedEmEnergy() / inJet->energy();
    pPFJet->ChargedHadronFraction = inJet->chargedHadronEnergy() / inJet->energy();
    pPFJet->NeutralEMFraction = inJet->neutralEmEnergy() / inJet->energy();
    pPFJet->NeutralHadronFraction = inJet->neutralHadronEnergy() / inJet->energy();

    pPFJet->TrackCountingHighEffBJetTagsDisc = (*(hTrackCountingHighEffBJetTags.product()))[jetBaseRef];
    pPFJet->TrackCountingHighPurBJetTagsDisc  = (*(hTrackCountingHighPurBJetTags.product()))[jetBaseRef];
    pPFJet->SoftElectronByPtBJetTagsDisc = (*(hSoftElectronByPtBJetTags.product()))[jetBaseRef];
    pPFJet->SoftElectronByIP3dBJetTagsDisc = (*(hSoftElectronByIP3dBJetTags.product()))[jetBaseRef];
    pPFJet->SoftMuonByPtBJetTagsDisc = (*(hSoftMuonByPtBJetTags.product()))[jetBaseRef];
    pPFJet->SoftMuonByIP3dBJetTagsDisc = (*(hSoftMuonByIP3dBJetTags.product()))[jetBaseRef];
    pPFJet->SoftMuonBJetTagsDisc = (*(hSoftMuonBJetTags.product()))[jetBaseRef];
    pPFJet->SimpleSecondaryVertexHighPurBJetTagsDisc = (*(hSimpleSecondaryVertexHighPurBJetTags.product()))[jetBaseRef];
    pPFJet->SimpleSecondaryVertexHighEffBJetTagsDisc = (*(hSimpleSecondaryVertexHighEffBJetTags.product()))[jetBaseRef];
    pPFJet->CombinedSecondaryVertexBJetTagsDisc = (*(hCombinedSecondaryVertexBJetTags.product()))[jetBaseRef];
    pPFJet->CombinedSecondaryVertexMVABJetTagsDisc = (*(hCombinedSecondaryVertexMVABJetTags.product()))[jetBaseRef];
    pPFJet->JetProbabilityBJetTagsDisc = (*(hJetProbabilityBJetTags.product()))[jetBaseRef];
    pPFJet->JetBProbabilityBJetTagsDisc = (*(hJetBProbabilityBJetTags.product()))[jetBaseRef];

    pPFJet->JetArea = inJet->jetArea();

    //jet flavor matching
    if (fUseGen) {
      pPFJet->matchedPdgId = (*PFJetFlavMatch)[edm::RefToBase<reco::Jet>(jetRef)].getFlavour();
    } else {
      pPFJet->matchedPdgId = -999;
    }

  }


  //****************************************************************************************************
  // MET
  //****************************************************************************************************
  Handle<reco::PFMETCollection> hMetProduct;
  GetProduct("pfMet", hMetProduct, event);
  const reco::PFMETCollection inMets = *(hMetProduct.product());  
  double pfMetX = 0;
  double pfMetY = 0;
  for (reco::PFMETCollection::const_iterator inMet = inMets.begin(); 
       inMet != inMets.end(); ++inMet) {
    pfMetX = inMet->px();
    pfMetY = inMet->py();
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

//   fEventInfo.triggerBits  = trigbits;
//   fEventInfo.l1triggerBits  = l1trigbits;

  fEventInfo.nGoodPV      = NGoodPV;
  if (thevtx) {
    fEventInfo.pvx          = thevtx->x();
    fEventInfo.pvy          = thevtx->y();
    fEventInfo.pvz          = thevtx->z();
  }
  if (inBeamSpot) {
    fEventInfo.bsx          = inBeamSpot->x0();
    fEventInfo.bsy          = inBeamSpot->y0();
    fEventInfo.bsz          = inBeamSpot->z0();
  }

  fEventInfo.genVertexX     = genVertexX;
  fEventInfo.genVertexY     = genVertexY;
  fEventInfo.genVertexZ     = genVertexZ;

  fEventInfo.pfMEx = pfMetX;
  fEventInfo.pfMEy = pfMetY;
  fEventInfo.pfTrackMEx = trackPFMetX;
  fEventInfo.pfTrackMEy = trackPFMetY;

  ///Fill the branches
  fEventTree->Fill();

  //cleanup
  if (lazyTools) delete lazyTools;
  if (local) delete local;  
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
CMSNtupler::beginJob()
{


  cout<<" CMSNtupler::beginJob( " <<endl; 
  
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

  cout<<" CMSNtupler tree branches defined.. " <<endl; 
  
  

}







double CMSNtupler::DeltaPhi(double v1, double v2)
{ // Computes the correctly normalized phi difference
  // v1, v2 = phi of object 1 and 2
  double diff = v1 - v2;
  
  while (diff >acos(-1)) diff -= 2*acos(-1);
  while (diff <= -acos(-1)) diff += 2*acos(-1);
  
  return diff; 
}


double CMSNtupler::GetDeltaR(double eta1, double eta2, double phi1, double phi2){ 
  // Computes the DeltaR of two objects from their eta and phi values
  
  return sqrt( (eta1-eta2)*(eta1-eta2) + DeltaPhi(phi1, phi2)*DeltaPhi(phi1, phi2) );
  
}


// ------------ method called once each job just after ending the event loop  ------------
void 
CMSNtupler::endJob() 
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
CMSNtupler::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
CMSNtupler::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CMSNtupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CMSNtupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CMSNtupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CMSNtupler);

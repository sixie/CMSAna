#ifndef JetTree_H
#define JetTree_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include <cmath>
#include "assert.h"

namespace cmsana
{
  class JetTree {

    public:

      /// bit map
      /// DON'T CHANGE ORDER

      //*******************************************
      //=== JetTriggerBits  ====
      //*******************************************
      enum JetTriggerBits { kJetTrigger_Jet                                    = 0x000001
      };

      /// variables
      Float_t                 fWeight;
      UInt_t                  fRunNumber;
      UInt_t                  fLumiSectionNumber;
      UInt_t                  fEventNumber;
      Bool_t                  fJetEventNumberParity;
      Float_t                 fJetGenPt;
      Float_t                 fJetGenEta; 
      Float_t                 fJetGenPhi;
      Float_t                 fJetRawPt;
      Float_t                 fJetPt; 
      Float_t                 fJetEta; 
      Float_t                 fJetPhi; 
      UInt_t                  fJetTriggerBit;
      Float_t                 fRho; 
      Float_t                 fNVertices; 
      Float_t                 fJetTrackCountingHighEffBJetTagsDisc;
      Float_t                 fJetTrackCountingHighPurBJetTagsDisc;
      Float_t                 fJetSimpleSecondaryVertexHighPurBJetTagsDisc;
      Float_t                 fJetSimpleSecondaryVertexHighEffBJetTagsDisc;
      Float_t                 fJetCombinedSecondaryVertexBJetTagsDisc;
      Float_t                 fJetCombinedSecondaryVertexMVABJetTagsDisc;
      Float_t                 fJetProbabilityBJetTagsDisc;
      Float_t                 fJetBProbabilityBJetTagsDisc;
      Float_t                 fJetArea;
      Int_t                   fJetMatchedPdgId;

    public:
      /// this is the main element
      TTree *tree_;
      TFile *f_;
  
      /// hold the names of variables to facilitate things (filled during Init)
      std::vector<std::string> variables_;

      /// default constructor  
      JetTree()  {};
      /// default destructor
      ~JetTree(){ 
        if (f_) f_->Close();  
      };
    
      /// initialize varibles and fill list of available variables
      void InitVariables() {
        fWeight			                     = 0.0;
        fRunNumber		                     = 0.0;
        fLumiSectionNumber	                     = 0.0;
        fEventNumber		                     = 0.0;
        fJetEventNumberParity 	                     = 0.0;
        fJetGenPt 		                     = 0.0;
        fJetGenEta 		                     = 0.0;
        fJetGenPhi 		                     = 0.0;
        fJetRawPt 		                     = 0.0;
        fJetPt 			                     = 0.0;
        fJetEta 		                     = 0.0;
        fJetPhi 		                     = 0.0;
        fJetTriggerBit		                     = 0.0;
        fRho  			                     = 0.0;
        fNVertices 		                     = 0.0;
        fJetTrackCountingHighEffBJetTagsDisc         = 0.0;
        fJetTrackCountingHighPurBJetTagsDisc         = 0.0;
        fJetSimpleSecondaryVertexHighPurBJetTagsDisc = 0.0;
        fJetSimpleSecondaryVertexHighEffBJetTagsDisc = 0.0;
        fJetCombinedSecondaryVertexBJetTagsDisc      = 0.0;
        fJetCombinedSecondaryVertexMVABJetTagsDisc   = 0.0;
        fJetProbabilityBJetTagsDisc                  = 0.0;
        fJetBProbabilityBJetTagsDisc                 = 0.0;
        fJetArea                                     = 0.0;
        fJetMatchedPdgId                             = 0;
      }
    
      /// load a JetTree
      void LoadTree(const char* file){
        f_ = TFile::Open(file);
        assert(f_);
        tree_ = dynamic_cast<TTree*>(f_->Get("Jets"));
        assert(tree_);
      }
    
      /// create a JetTree
      void CreateTree(){
        tree_ = new TTree("Jets","Jettons");
        f_ = 0;

        //book the branches
        tree_->Branch("weight",&fWeight,"weight/F");
        tree_->Branch("run",&fRunNumber,"run/i");
        tree_->Branch("lumi",&fLumiSectionNumber,"lumi/i");
        tree_->Branch("event",&fEventNumber,"event/i");
        tree_->Branch("EventNumberParity",&fJetEventNumberParity,"EventNumberParity/O"); 
        tree_->Branch("genpt",&fJetGenPt,"genpt/F"); 
        tree_->Branch("geneta",&fJetGenEta,"geneta/F"); 
        tree_->Branch("genphi",&fJetGenPhi,"genphi/F"); 
        tree_->Branch("rawpt",&fJetRawPt,"rawpt/F"); 
        tree_->Branch("pt",&fJetPt,"pt/F"); 
        tree_->Branch("eta",&fJetEta,"eta/F"); 
        tree_->Branch("phi",&fJetPhi,"phi/F"); 
        tree_->Branch("triggerBit",&fJetTriggerBit,"triggerBit/i"); 
        tree_->Branch("rho",&fRho,"rho/F"); 
        tree_->Branch("vertices",&fNVertices,"vertices/F"); 
        tree_->Branch("TrackCountingHighEffBJetTagsDisc",&fJetTrackCountingHighEffBJetTagsDisc,"TrackCountingHighEffBJetTagsDisc/F"); 
        tree_->Branch("TrackCountingHighPurBJetTagsDisc",&fJetTrackCountingHighPurBJetTagsDisc,"TrackCountingHighPurBJetTagsDisc/F"); 
        tree_->Branch("SimpleSecondaryVertexHighPurBJetTagsDisc",&fJetSimpleSecondaryVertexHighPurBJetTagsDisc,"SimpleSecondaryVertexHighPurBJetTagsDisc/F"); 
        tree_->Branch("SimpleSecondaryVertexHighEffBJetTagsDisc",&fJetSimpleSecondaryVertexHighEffBJetTagsDisc,"SimpleSecondaryVertexHighEffBJetTagsDisc/F"); 
        tree_->Branch("CombinedSecondaryVertexBJetTagsDisc",&fJetCombinedSecondaryVertexBJetTagsDisc,"CombinedSecondaryVertexBJetTagsDisc/F"); 
        tree_->Branch("CombinedSecondaryVertexMVABJetTagsDisc",&fJetCombinedSecondaryVertexMVABJetTagsDisc,"CombinedSecondaryVertexMVABJetTagsDisc/F"); 
        tree_->Branch("ProbabilityBJetTagsDisc",&fJetProbabilityBJetTagsDisc,"ProbabilityBJetTagsDisc/F"); 
        tree_->Branch("BProbabilityBJetTagsDisc",&fJetBProbabilityBJetTagsDisc,"BProbabilityBJetTagsDisc/F"); 
        tree_->Branch("Area",&fJetArea,"Area/F"); 
        tree_->Branch("MatchedPdgId",&fJetMatchedPdgId,"MatchedPdgId/I"); 

      } 

      // initialze a JettonTree
      void InitTree(){
        assert(tree_);
        // don't forget to set pointers to zero before you set address
        // or you will fully appreciate that "ROOT sucks" :)
        InitVariables();
        //Set branch address
        Int_t currentState = gErrorIgnoreLevel;

        tree_->SetBranchAddress("weight",&fWeight);
        tree_->SetBranchAddress("run",&fRunNumber);
        tree_->SetBranchAddress("lumi",&fLumiSectionNumber);
        tree_->SetBranchAddress("event",&fEventNumber);
        tree_->SetBranchAddress("EventNumberParity",&fJetEventNumberParity);
        tree_->SetBranchAddress("genpt",&fJetGenPt);
        tree_->SetBranchAddress("geneta",&fJetGenEta);
        tree_->SetBranchAddress("genphi",&fJetGenPhi);
        tree_->SetBranchAddress("rawpt",&fJetRawPt);
        tree_->SetBranchAddress("pt",&fJetPt);
        tree_->SetBranchAddress("eta",&fJetEta);
        tree_->SetBranchAddress("phi",&fJetPhi);
        tree_->SetBranchAddress("triggerBit",&fJetTriggerBit);
        tree_->SetBranchAddress("rho",&fRho);
        tree_->SetBranchAddress("vertices",&fNVertices);
        tree_->SetBranchAddress("TrackCountingHighEffBJetTagsDisc",&fJetTrackCountingHighEffBJetTagsDisc);
        tree_->SetBranchAddress("TrackCountingHighPurBJetTagsDisc",&fJetTrackCountingHighPurBJetTagsDisc);
        tree_->SetBranchAddress("SimpleSecondaryVertexHighPurBJetTagsDisc",&fJetSimpleSecondaryVertexHighPurBJetTagsDisc);
        tree_->SetBranchAddress("SimpleSecondaryVertexHighEffBJetTagsDisc",&fJetSimpleSecondaryVertexHighEffBJetTagsDisc);
        tree_->SetBranchAddress("CombinedSecondaryVertexBJetTagsDisc",&fJetCombinedSecondaryVertexBJetTagsDisc);
        tree_->SetBranchAddress("CombinedSecondaryVertexMVABJetTagsDisc",&fJetCombinedSecondaryVertexMVABJetTagsDisc);
        tree_->SetBranchAddress("ProbabilityBJetTagsDisc",&fJetProbabilityBJetTagsDisc);
        tree_->SetBranchAddress("BProbabilityBJetTagsDisc",&fJetBProbabilityBJetTagsDisc);
        tree_->SetBranchAddress("Area",&fJetArea);
        tree_->SetBranchAddress("MatchedPdgId",&fJetMatchedPdgId);

        gErrorIgnoreLevel = currentState;
      }

  }; 

}

#endif

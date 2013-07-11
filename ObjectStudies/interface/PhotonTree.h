#ifndef PhotonTree_H
#define PhotonTree_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include <cmath>
#include "assert.h"

namespace cmsana
{
  class PhotonTree {

    public:

      /// bit map
      /// DON'T CHANGE ORDER

      //*******************************************
      //=== PhotonTriggerBits  ====
      //*******************************************
      enum PhotonTriggerBits { kPhoTrigger_Pho                                    = 0x000001
      };

      /// variables
      Float_t                 fWeight;
      UInt_t                  fRunNumber;
      UInt_t                  fLumiSectionNumber;
      UInt_t                  fEventNumber;
      Bool_t                  fPhoEventNumberParity;
      Float_t                 fPhoGenPt;
      Float_t                 fPhoGenEta; 
      Float_t                 fPhoGenPhi;
      Float_t                 fPhoPt; 
      Float_t                 fPhoEta; 
      Float_t                 fPhoPhi; 
      Float_t                 fPhoSCEt; 
      Float_t                 fPhoSCEta; 
      Float_t                 fPhoSCPhi; 
      UInt_t                  fPhoTriggerBit;
      Float_t                 fRho; 
      Float_t                 fNVertices; 

      Bool_t                  fPhoHasPixelSeed;
      Bool_t                  fPhoIsConversion;
      Bool_t                  fPhoPassEleVeto;
      Float_t                 fPhoHoverE; 
      Float_t                 fPhoHoverESingleTower; 

      //shower shape
      Float_t                 fPhoSigmaIEtaIEta;
      Float_t                 fPhoSigmaIPhiIPhi;
      Float_t                 fPhoSigmaIEtaIPhi;
      Float_t                 fPhoSCEtaWidth;
      Float_t                 fPhoSCPhiWidth;
      Float_t                 fPhoR9;

      //Isolation Variables
      Float_t                 fPhoTrkIso03; 
      Float_t                 fPhoEMIso03; 
      Float_t                 fPhoHadIso03; 
      Float_t                 fPhoPFIso03;
      Float_t                 fPhoPFIso04;
      Float_t                 fChargedIso_DR0p0To0p1;
      Float_t                 fChargedIso_DR0p1To0p2;
      Float_t                 fChargedIso_DR0p2To0p3;
      Float_t                 fChargedIso_DR0p3To0p4;
      Float_t                 fChargedIso_DR0p4To0p5;
      Float_t                 fGammaIso_DR0p0To0p1;
      Float_t                 fGammaIso_DR0p1To0p2;
      Float_t                 fGammaIso_DR0p2To0p3;
      Float_t                 fGammaIso_DR0p3To0p4;
      Float_t                 fGammaIso_DR0p4To0p5;
      Float_t                 fNeutralHadronIso_DR0p0To0p1;
      Float_t                 fNeutralHadronIso_DR0p1To0p2;
      Float_t                 fNeutralHadronIso_DR0p2To0p3;
      Float_t                 fNeutralHadronIso_DR0p3To0p4;
      Float_t                 fNeutralHadronIso_DR0p4To0p5;

      //Regression Variables
      Float_t                 fSCRawEnergy;
      Float_t                 fPreShowerOverRaw; 
      Float_t                 fNClusters;
      Float_t                 fEtaSeed;
      Float_t                 fPhiSeed;
      Float_t                 fESeed;
      Float_t                 fE3x3Seed;
      Float_t                 fE5x5Seed;
      Float_t                 fEMaxSeed;
      Float_t                 fE2ndSeed;
      Float_t                 fETopSeed;
      Float_t                 fEBottomSeed;
      Float_t                 fELeftSeed;
      Float_t                 fERightSeed;
      Float_t                 fE2x5MaxSeed;
      Float_t                 fE2x5TopSeed;
      Float_t                 fE2x5BottomSeed;
      Float_t                 fE2x5LeftSeed;
      Float_t                 fE2x5RightSeed;
      Float_t                 fIEtaSeed;
      Float_t                 fIPhiSeed;
      Float_t                 fEtaCrySeed;
      Float_t                 fPhiCrySeed;
      Float_t                 fGeneratedEnergy;
      Float_t                 fPhoMomentum_Regression_V0;
      Float_t                 fPhoMomentumError_Regression_V0;

    public:
      /// this is the main element
      TTree *tree_;
      TFile *f_;
  
      /// hold the names of variables to facilitate things (filled during Init)
      std::vector<std::string> variables_;

      /// default constructor  
      PhotonTree()  {};
      /// default destructor
      ~PhotonTree(){ 
        if (f_) f_->Close();  
      };
    
      /// initialize varibles and fill list of available variables
      void InitVariables() {
        fWeight			       = 0.0;
        fRunNumber		       = 0.0;
        fLumiSectionNumber	       = 0.0;
        fEventNumber		       = 0.0;
        fPhoEventNumberParity 	       = 0.0;
        fPhoGenPt 		       = 0.0;
        fPhoGenEta 		       = 0.0;
        fPhoGenPhi 		       = 0.0;
        fPhoPt 			       = 0.0;
        fPhoEta 		       = 0.0;
        fPhoPhi 		       = 0.0;
        fPhoSCEt 		       = 0.0;
        fPhoSCEta 		       = 0.0;
        fPhoSCPhi 		       = 0.0;
        fPhoTriggerBit		       = 0.0;
        fRho  			       = 0.0;
        fNVertices 		       = 0.0;
        fPhoHasPixelSeed               = false;
        fPhoIsConversion	       = false;
        fPhoPassEleVeto		       = false;
        fPhoHoverE 		       = 0.0;
        fPhoHoverESingleTower 	       = 0.0;
        fPhoSigmaIEtaIEta 	       = 0.0;
        fPhoSigmaIPhiIPhi 	       = 0.0;
        fPhoSigmaIEtaIPhi	       = 0.0;
        fPhoSCEtaWidth		       = 0.0;
        fPhoSCPhiWidth		       = 0.0;
        fPhoR9			       = 0.0;
        fPhoTrkIso03  		       = 0.0;
        fPhoEMIso03 		       = 0.0;
        fPhoHadIso03  		       = 0.0;
        fPhoPFIso03 		       = 0.0;
        fPhoPFIso04 		       = 0.0;
        fChargedIso_DR0p0To0p1	       = 0.0;
        fChargedIso_DR0p1To0p2	       = 0.0;
        fChargedIso_DR0p2To0p3	       = 0.0;
        fChargedIso_DR0p3To0p4	       = 0.0;
        fChargedIso_DR0p4To0p5	       = 0.0;
        fGammaIso_DR0p0To0p1	       = 0.0;
        fGammaIso_DR0p1To0p2	       = 0.0;
        fGammaIso_DR0p2To0p3	       = 0.0;
        fGammaIso_DR0p3To0p4	       = 0.0;
        fGammaIso_DR0p4To0p5	       = 0.0;
        fNeutralHadronIso_DR0p0To0p1   = 0.0;
        fNeutralHadronIso_DR0p1To0p2   = 0.0;
        fNeutralHadronIso_DR0p2To0p3   = 0.0;
        fNeutralHadronIso_DR0p3To0p4   = 0.0;
        fNeutralHadronIso_DR0p4To0p5   = 0.0;
        fSCRawEnergy                   = 0.0;
        fPreShowerOverRaw	       = 0.0;
        fNClusters                     = 0.0;
        fEtaSeed                       = 0.0;
        fPhiSeed                       = 0.0;
        fESeed                         = 0.0;
        fE3x3Seed                      = 0.0;
        fE5x5Seed                      = 0.0;
        fEMaxSeed                      = 0.0;
        fE2ndSeed                      = 0.0;
        fETopSeed                      = 0.0;
        fEBottomSeed                   = 0.0;
        fELeftSeed                     = 0.0;
        fERightSeed                    = 0.0;
        fE2x5MaxSeed                   = 0.0;
        fE2x5TopSeed                   = 0.0;
        fE2x5BottomSeed                = 0.0;
        fE2x5LeftSeed                  = 0.0;
        fE2x5RightSeed                 = 0.0;
        fIEtaSeed                      = 0.0;
        fIPhiSeed                      = 0.0;
        fEtaCrySeed                    = 0.0;
        fPhiCrySeed                    = 0.0;
        fGeneratedEnergy               = 0.0;
        fPhoMomentum_Regression_V0     = 0.0;
        fPhoMomentumError_Regression_V0= 0.0;
      }
    
      /// load a PhotonTree
      void LoadTree(const char* file){
        f_ = TFile::Open(file);
        assert(f_);
        tree_ = dynamic_cast<TTree*>(f_->Get("Photons"));
        assert(tree_);
      }
    
      /// create a PhotonTree
      void CreateTree(){
        tree_ = new TTree("Photons","Photons");
        f_ = 0;

        //book the branches
        tree_->Branch("weight",&fWeight,"weight/F");
        tree_->Branch("run",&fRunNumber,"run/i");
        tree_->Branch("lumi",&fLumiSectionNumber,"lumi/i");
        tree_->Branch("event",&fEventNumber,"event/i");
        tree_->Branch("EventNumberParity",&fPhoEventNumberParity,"EventNumberParity/O"); 
        tree_->Branch("genpt",&fPhoGenPt,"genpt/F"); 
        tree_->Branch("geneta",&fPhoGenEta,"geneta/F"); 
        tree_->Branch("genphi",&fPhoGenPhi,"genphi/F"); 
        tree_->Branch("pt",&fPhoPt,"pt/F"); 
        tree_->Branch("eta",&fPhoEta,"eta/F"); 
        tree_->Branch("phi",&fPhoPhi,"phi/F"); 
        tree_->Branch("scEt",&fPhoSCEt,"scEt/F"); 
        tree_->Branch("scEta",&fPhoSCEta,"scEta/F"); 
        tree_->Branch("scPhi",&fPhoSCPhi,"scPhi/F"); 
        tree_->Branch("triggerBit",&fPhoTriggerBit,"triggerBit/i"); 
        tree_->Branch("rho",&fRho,"rho/F"); 
        tree_->Branch("vertices",&fNVertices,"vertices/F"); 
        tree_->Branch("hasPixelSeed",&fPhoHasPixelSeed,"hasPixelSeed/O"); 
        tree_->Branch("isConversion",&fPhoIsConversion,"isConversion/O"); 
        tree_->Branch("passEleVeto",&fPhoPassEleVeto,"passEleVeto/O"); 
        tree_->Branch("HoverE",&fPhoHoverE,"HoverE/F"); 
        tree_->Branch("HoverESingleTower",&fPhoHoverESingleTower,"HoverESingleTower/F"); 
        tree_->Branch("see",&fPhoSigmaIEtaIEta,"see/F"); 
        tree_->Branch("spp",&fPhoSigmaIPhiIPhi,"spp/F"); 
        tree_->Branch("sep",&fPhoSigmaIEtaIPhi,"sep/F"); 
        tree_->Branch("etawidth",&fPhoSCEtaWidth,"etawidth/F"); 
        tree_->Branch("phiwidth",&fPhoSCPhiWidth,"phiwidth/F"); 
        tree_->Branch("R9",&fPhoR9,"R9/F"); 
        tree_->Branch("trkIso03",&fPhoTrkIso03,"trkIso03/F"); 
        tree_->Branch("ecalIso03",&fPhoEMIso03,"ecalIso03/F"); 
        tree_->Branch("hcalIso03",&fPhoHadIso03,"hcalIso03/F"); 
        tree_->Branch("pfIso03",&fPhoPFIso03,"pfIso03/F"); 
        tree_->Branch("pfIso04",&fPhoPFIso04,"pfIso04/F"); 
        tree_->Branch("ChargedIso_DR0p0To0p1",&fChargedIso_DR0p0To0p1,"ChargedIso_DR0p0To0p1/F");
        tree_->Branch("ChargedIso_DR0p1To0p2",&fChargedIso_DR0p1To0p2,"ChargedIso_DR0p1To0p2/F");
        tree_->Branch("ChargedIso_DR0p2To0p3",&fChargedIso_DR0p2To0p3,"ChargedIso_DR0p2To0p3/F");
        tree_->Branch("ChargedIso_DR0p3To0p4",&fChargedIso_DR0p3To0p4,"ChargedIso_DR0p3To0p4/F");
        tree_->Branch("ChargedIso_DR0p4To0p5",&fChargedIso_DR0p4To0p5,"ChargedIso_DR0p4To0p5/F");
        tree_->Branch("GammaIso_DR0p0To0p1",&fGammaIso_DR0p0To0p1,"GammaIso_DR0p0To0p1/F");
        tree_->Branch("GammaIso_DR0p1To0p2",&fGammaIso_DR0p1To0p2,"GammaIso_DR0p1To0p2/F");
        tree_->Branch("GammaIso_DR0p2To0p3",&fGammaIso_DR0p2To0p3,"GammaIso_DR0p2To0p3/F");
        tree_->Branch("GammaIso_DR0p3To0p4",&fGammaIso_DR0p3To0p4,"GammaIso_DR0p3To0p4/F");
        tree_->Branch("GammaIso_DR0p4To0p5",&fGammaIso_DR0p4To0p5,"GammaIso_DR0p4To0p5/F");
        tree_->Branch("NeutralHadronIso_DR0p0To0p1",&fNeutralHadronIso_DR0p0To0p1,"NeutralHadronIso_DR0p0To0p1/F");
        tree_->Branch("NeutralHadronIso_DR0p1To0p2",&fNeutralHadronIso_DR0p1To0p2,"NeutralHadronIso_DR0p1To0p2/F");
        tree_->Branch("NeutralHadronIso_DR0p2To0p3",&fNeutralHadronIso_DR0p2To0p3,"NeutralHadronIso_DR0p2To0p3/F");
        tree_->Branch("NeutralHadronIso_DR0p3To0p4",&fNeutralHadronIso_DR0p3To0p4,"NeutralHadronIso_DR0p3To0p4/F");
        tree_->Branch("NeutralHadronIso_DR0p4To0p5",&fNeutralHadronIso_DR0p4To0p5,"NeutralHadronIso_DR0p4To0p5/F");
        tree_->Branch("SCRawEnergy",&fSCRawEnergy,"SCRawEnergy/F");
        tree_->Branch("PreShowerOverRaw",&fPreShowerOverRaw,"PreShowerOverRaw/F");
        tree_->Branch("NClusters",&fNClusters,"NClusters/F");
        tree_->Branch("EtaSeed",&fEtaSeed,"EtaSeed/F");
        tree_->Branch("PhiSeed",&fPhiSeed,"PhiSeed/F");
        tree_->Branch("ESeed",&fESeed,"ESeed/F");
        tree_->Branch("E3x3Seed",&fE3x3Seed,"E3x3Seed/F");
        tree_->Branch("E5x5Seed",&fE5x5Seed,"E5x5Seed/F");
        tree_->Branch("EMaxSeed",&fEMaxSeed,"EMaxSeed/F");
        tree_->Branch("E2ndSeed",&fE2ndSeed,"E2ndSeed/F");
        tree_->Branch("ETopSeed",&fETopSeed,"ETopSeed/F");
        tree_->Branch("EBottomSeed",&fEBottomSeed,"EBottomSeed/F");
        tree_->Branch("ELeftSeed",&fELeftSeed,"ELeftSeed/F");
        tree_->Branch("ERightSeed",&fERightSeed,"ERightSeed/F");
        tree_->Branch("E2x5MaxSeed",&fE2x5MaxSeed,"E2x5MaxSeed/F");
        tree_->Branch("E2x5TopSeed",&fE2x5TopSeed,"E2x5TopSeed/F");
        tree_->Branch("E2x5BottomSeed",&fE2x5BottomSeed,"E2x5BottomSeed/F");
        tree_->Branch("E2x5LeftSeed",&fE2x5LeftSeed,"E2x5LeftSeed/F");
        tree_->Branch("E2x5RightSeed",&fE2x5RightSeed,"E2x5RightSeed/F");
        tree_->Branch("IEtaSeed",&fIEtaSeed,"IEtaSeed/F");
        tree_->Branch("IPhiSeed",&fIPhiSeed,"IPhiSeed/F");
        tree_->Branch("EtaCrySeed",&fEtaCrySeed,"EtaCrySeed/F");
        tree_->Branch("PhiCrySeed",&fPhiCrySeed,"PhiCrySeed/F");
        tree_->Branch("GeneratedEnergy",&fGeneratedEnergy,"GeneratedEnergy/F");
        tree_->Branch("EleMomentum_Regression_V0",&fPhoMomentum_Regression_V0,"EleMomentum_Regression_V0/F");
        tree_->Branch("EleMomentumError_Regression_V0",&fPhoMomentumError_Regression_V0,"EleMomentumError_Regression_V0/F");

      } 

      // initialze a PhotonTree
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
        tree_->SetBranchAddress("EventNumberParity",&fPhoEventNumberParity);
        tree_->SetBranchAddress("genpt",&fPhoGenPt);
        tree_->SetBranchAddress("geneta",&fPhoGenEta);
        tree_->SetBranchAddress("genphi",&fPhoGenPhi);
        tree_->SetBranchAddress("pt",&fPhoPt);
        tree_->SetBranchAddress("eta",&fPhoEta);
        tree_->SetBranchAddress("phi",&fPhoPhi);
        tree_->SetBranchAddress("scEt",&fPhoSCEt);
        tree_->SetBranchAddress("scEta",&fPhoSCEta);
        tree_->SetBranchAddress("scPhi",&fPhoSCPhi);
        tree_->SetBranchAddress("triggerBit",&fPhoTriggerBit);
        tree_->SetBranchAddress("rho",&fRho);
        tree_->SetBranchAddress("vertices",&fNVertices);
        tree_->SetBranchAddress("hasPixelSeed",&fPhoHasPixelSeed);
        tree_->SetBranchAddress("isConversion",&fPhoIsConversion);
        tree_->SetBranchAddress("passEleVeto",&fPhoPassEleVeto);
        tree_->SetBranchAddress("HoverE",&fPhoHoverE);
        tree_->SetBranchAddress("HoverESingleTower",&fPhoHoverESingleTower);
        tree_->SetBranchAddress("see",&fPhoSigmaIEtaIEta);
        tree_->SetBranchAddress("spp",&fPhoSigmaIPhiIPhi);
        tree_->SetBranchAddress("sep",&fPhoSigmaIEtaIPhi);
        tree_->SetBranchAddress("etawidth",&fPhoSCEtaWidth);
        tree_->SetBranchAddress("phiwidth",&fPhoSCPhiWidth);
        tree_->SetBranchAddress("R9",&fPhoR9);
        tree_->SetBranchAddress("trkIso03",&fPhoTrkIso03);
        tree_->SetBranchAddress("ecalIso03",&fPhoEMIso03);
        tree_->SetBranchAddress("hcalIso03",&fPhoHadIso03);
        tree_->SetBranchAddress("pfIso03",&fPhoPFIso03);
        tree_->SetBranchAddress("pfIso04",&fPhoPFIso04);
        tree_->SetBranchAddress("ChargedIso_DR0p0To0p1",&fChargedIso_DR0p0To0p1);
        tree_->SetBranchAddress("ChargedIso_DR0p1To0p2",&fChargedIso_DR0p1To0p2);
        tree_->SetBranchAddress("ChargedIso_DR0p2To0p3",&fChargedIso_DR0p2To0p3);
        tree_->SetBranchAddress("ChargedIso_DR0p3To0p4",&fChargedIso_DR0p3To0p4);
        tree_->SetBranchAddress("ChargedIso_DR0p4To0p5",&fChargedIso_DR0p4To0p5);
        tree_->SetBranchAddress("GammaIso_DR0p0To0p1",&fGammaIso_DR0p0To0p1);
        tree_->SetBranchAddress("GammaIso_DR0p1To0p2",&fGammaIso_DR0p1To0p2);
        tree_->SetBranchAddress("GammaIso_DR0p2To0p3",&fGammaIso_DR0p2To0p3);
        tree_->SetBranchAddress("GammaIso_DR0p3To0p4",&fGammaIso_DR0p3To0p4);
        tree_->SetBranchAddress("GammaIso_DR0p4To0p5",&fGammaIso_DR0p4To0p5);
        tree_->SetBranchAddress("NeutralHadronIso_DR0p0To0p1",&fNeutralHadronIso_DR0p0To0p1);
        tree_->SetBranchAddress("NeutralHadronIso_DR0p1To0p2",&fNeutralHadronIso_DR0p1To0p2);
        tree_->SetBranchAddress("NeutralHadronIso_DR0p2To0p3",&fNeutralHadronIso_DR0p2To0p3);
        tree_->SetBranchAddress("NeutralHadronIso_DR0p3To0p4",&fNeutralHadronIso_DR0p3To0p4);
        tree_->SetBranchAddress("NeutralHadronIso_DR0p4To0p5",&fNeutralHadronIso_DR0p4To0p5);
        tree_->SetBranchAddress("SCRawEnergy",&fSCRawEnergy);
        tree_->SetBranchAddress("PreShowerOverRaw",&fPreShowerOverRaw);
        tree_->SetBranchAddress("NClusters",&fNClusters);
        tree_->SetBranchAddress("EtaSeed",&fEtaSeed);
        tree_->SetBranchAddress("PhiSeed",&fPhiSeed);
        tree_->SetBranchAddress("ESeed",&fESeed);
        tree_->SetBranchAddress("E3x3Seed",&fE3x3Seed);
        tree_->SetBranchAddress("E5x5Seed",&fE5x5Seed);
        tree_->SetBranchAddress("EMaxSeed",&fEMaxSeed);
        tree_->SetBranchAddress("E2ndSeed",&fE2ndSeed);
        tree_->SetBranchAddress("ETopSeed",&fETopSeed);
        tree_->SetBranchAddress("EBottomSeed",&fEBottomSeed);
        tree_->SetBranchAddress("ELeftSeed",&fELeftSeed);
        tree_->SetBranchAddress("ERightSeed",&fERightSeed);
        tree_->SetBranchAddress("E2x5MaxSeed",&fE2x5MaxSeed);
        tree_->SetBranchAddress("E2x5TopSeed",&fE2x5TopSeed);
        tree_->SetBranchAddress("E2x5BottomSeed",&fE2x5BottomSeed);
        tree_->SetBranchAddress("E2x5LeftSeed",&fE2x5LeftSeed);
        tree_->SetBranchAddress("E2x5RightSeed",&fE2x5RightSeed);
        tree_->SetBranchAddress("IEtaSeed",&fIEtaSeed);
        tree_->SetBranchAddress("IPhiSeed",&fIPhiSeed);
        tree_->SetBranchAddress("EtaCrySeed",&fEtaCrySeed);
        tree_->SetBranchAddress("PhiCrySeed",&fPhiCrySeed);
        tree_->SetBranchAddress("GeneratedEnergy",&fGeneratedEnergy);
        tree_->SetBranchAddress("EleMomentum_Regression_V0",&fPhoMomentum_Regression_V0);
        tree_->SetBranchAddress("EleMomentumError_Regression_V0",&fPhoMomentumError_Regression_V0);

        gErrorIgnoreLevel = currentState;
      }

  }; 

}

#endif

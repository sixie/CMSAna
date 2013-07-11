#ifndef HHToBBGGEventTree_H
#define HHToBBGGEventTree_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include "Math/LorentzVector.h"
#include <cmath>
#include "assert.h"
#include "CMSAna/DataTree/interface/Types.h"

namespace cmsana
{
  class HHToBBGGEventTree {

    public:
      enum SampleType {
        none = -1,
        data = 0,
        HHToBBGG = 1,
        ttHgg = 2,
        ZHgg = 3,
        ggHgg = 4,        
        ttbar = 5,
        diphotonjets = 6,
        qcd = 7
      };
   
      /// variables
      SampleType              sampletype;                  
      Float_t                 weight;
      UInt_t                  run;
      UInt_t                  lumi;
      UInt_t                  event;
      UInt_t                  npu;
      Float_t                 rho; 
      UInt_t                  nvtx;       
      Float_t                 pfmet;
      FourVector genpho1;
      FourVector genpho2;      
      FourVector genb1;
      FourVector genb2;
      FourVector genbjet1;
      FourVector genbjet2;

      FourVector pho1;
      FourVector pho2;
      FourVector diphoton;
      FourVector bjet1;
      FourVector bjet2;
      FourVector dibjet;
      FourVector bbgg;

      UInt_t                  nlep;
      UInt_t                  njets;
      UInt_t                  ncentraljets;
      Float_t                 DRgg;
      Float_t                 DRbb;
      Float_t                 minDRgb;
      Float_t                 HT;
      Float_t                 MET;
      Float_t                 pfTrackMET;

    public:
      /// this is the main element
      TTree *tree_;
      TFile *f_;
  
      /// hold the names of variables to facilitate things (filled during Init)
      std::vector<std::string> variables_;
      
      /// default constructor  
      HHToBBGGEventTree()  {
        genpho1Ptr     = &genpho1;
        genpho2Ptr     = &genpho2;
        genb1Ptr       = &genb1;
        genb2Ptr       = &genb2;
        genbjet1Ptr    = &genbjet1;
        genbjet2Ptr    = &genbjet2;
        pho1Ptr        = &pho1;
        pho2Ptr        = &pho2;
        diphotonPtr    = &diphoton;
        bjet1Ptr       = &bjet1;
        bjet2Ptr       = &bjet2;
        dibjetPtr      = &dibjet;
        bbggPtr        = &bbgg;

      };
      /// default destructor
      ~HHToBBGGEventTree(){ 
        if (f_) f_->Close();  
      };
    
      /// initialize varibles and fill list of available variables
      void InitVariables() {
        sampletype                     = HHToBBGGEventTree::none;
        weight			       = 0.0;
        run		               = 0.0;
        lumi	                       = 0.0;
        event     		       = 0.0;
        npu  			       = 0.0;
        rho  			       = 0.0;
        nvtx     		       = 0.0;
        pfmet     		       = 0.0;
        genpho1                        = FourVector();
        genpho2                        = FourVector();
        genb1                          = FourVector();
        genb2                          = FourVector();
        genbjet1                       = FourVector();
        genbjet2                       = FourVector();
        pho1                           = FourVector();
        pho2                           = FourVector();
        diphoton                       = FourVector();
        bjet1                          = FourVector();
        bjet2                          = FourVector();
        dibjet                         = FourVector();
        bbgg                           = FourVector();
        nlep                           = 0.0;
        njets                          = 0.0;
        ncentraljets                   = 0.0;
        DRgg                           = 0.0;
        DRbb                           = 0.0;
        minDRgb                        = 0.0;       
	HT                             = 0.0;        
	pfTrackMET                     = 0.0; 
      }
    
      /// load a HHToBBGGEventTree
      void LoadTree(const char* file){
        f_ = TFile::Open(file);
        assert(f_);
        tree_ = dynamic_cast<TTree*>(f_->Get("HHToBBGGEvent"));
        InitTree();
        assert(tree_);
      }
    
      /// create a HHToBBGGEventTree
      void CreateTree(){
        tree_ = new TTree("HHToBBGGEvent","HHToBBGGEvent");
        f_ = 0;

        //book the branches
        tree_->Branch("sampletype",&sampletype,"sampletype/I");
        tree_->Branch("weight",&weight,"weight/F");
        tree_->Branch("run",&run,"run/i");
        tree_->Branch("lumi",&lumi,"lumi/i");
        tree_->Branch("event",&event,"event/i");
        tree_->Branch("npu",&npu,"npu/i"); 
        tree_->Branch("rho",&rho,"rho/F"); 
        tree_->Branch("nvtx",&nvtx,"nvtx/i"); 
        tree_->Branch("pfmet",&pfmet,"pfmet/F"); 

        tree_->Branch("genpho1"     , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &genpho1Ptr);
        tree_->Branch("genpho2"     , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &genpho2Ptr);
        tree_->Branch("genb1"       , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &genb1Ptr);
        tree_->Branch("genb2"       , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &genb2Ptr);
        tree_->Branch("genbjet1"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &genbjet1Ptr);
        tree_->Branch("genbjet2"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &genbjet2Ptr);
        tree_->Branch("pho1"        , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &pho1Ptr);
        tree_->Branch("pho2"        , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &pho2Ptr);
        tree_->Branch("diphoton"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &diphotonPtr);
        tree_->Branch("bjet1"       , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &bjet1Ptr);
        tree_->Branch("bjet2"       , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &bjet2Ptr);
        tree_->Branch("dibjet"      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &dibjetPtr);
        tree_->Branch("bbgg"        , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &bbggPtr);
        tree_->Branch("nlep",&nlep,"nlep/i"); 
        tree_->Branch("njets",&njets,"njets/i"); 
        tree_->Branch("ncentraljets",&ncentraljets,"ncentraljets/i"); 
        tree_->Branch("DRgg",&DRgg,"DRgg/F"); 
        tree_->Branch("DRbb",&DRbb,"DRbb/F"); 
        tree_->Branch("minDRgb",&minDRgb,"minDRgb/F"); 
	tree_->Branch("HT",&HT,"HT/F"); 
	tree_->Branch("MET",&MET,"MET"); 
	tree_->Branch("pfTrackMET",&pfTrackMET,"pfTrackMET"); 

    
      }

      // initialze a HHToBBGGEventTree
      void InitTree(){
        assert(tree_);
        // don't forget to set pointers to zero before you set address
        // or you will fully appreciate that "ROOT sucks" :)
        InitVariables();
        //Set branch address
        Int_t currentState = gErrorIgnoreLevel;


        tree_->SetBranchAddress("sampletype",&sampletype);
        tree_->SetBranchAddress("weight",&weight);
        tree_->SetBranchAddress("run",&run);
        tree_->SetBranchAddress("lumi",&lumi);
        tree_->SetBranchAddress("event",&event);
        tree_->SetBranchAddress("npu",&npu);
        tree_->SetBranchAddress("rho",&rho);
        tree_->SetBranchAddress("nvtx",&nvtx);        
        tree_->SetBranchAddress("pfmet",&pfmet);        
        tree_->SetBranchAddress("genpho1",&genpho1Ptr);
        tree_->SetBranchAddress("genpho2",&genpho2Ptr);
        tree_->SetBranchAddress("genb1",&genb1Ptr);
        tree_->SetBranchAddress("genb2",&genb2Ptr);
        tree_->SetBranchAddress("genbjet1",&genbjet1Ptr);
        tree_->SetBranchAddress("genbjet2",&genbjet2Ptr);
        tree_->SetBranchAddress("pho1",&pho1Ptr);
        tree_->SetBranchAddress("pho2",&pho2Ptr);
        tree_->SetBranchAddress("diphoton",&diphotonPtr);
        tree_->SetBranchAddress("bjet1",&bjet1Ptr);
        tree_->SetBranchAddress("bjet2",&bjet2Ptr);
        tree_->SetBranchAddress("dibjet",&dibjetPtr);
        tree_->SetBranchAddress("bbgg",&bbggPtr);
        
        tree_->SetBranchAddress("nlep",&nlep);
        tree_->SetBranchAddress("njets",&njets);
        tree_->SetBranchAddress("ncentraljets",&ncentraljets);
        tree_->SetBranchAddress("DRgg",&DRgg);
        tree_->SetBranchAddress("DRbb",&DRbb);
        tree_->SetBranchAddress("minDRgb",&minDRgb);
	tree_->SetBranchAddress("HT",&HT);
	tree_->SetBranchAddress("MET",&MET);
	tree_->SetBranchAddress("pfTrackMET",&pfTrackMET);
        gErrorIgnoreLevel = currentState;
      }

    private:
      FourVector* genpho1Ptr;
      FourVector* genpho2Ptr;
      FourVector* genb1Ptr;
      FourVector* genb2Ptr;
      FourVector* genbjet1Ptr;
      FourVector* genbjet2Ptr;

      FourVector* pho1Ptr;
      FourVector* pho2Ptr;
      FourVector* diphotonPtr;
      FourVector* bjet1Ptr;
      FourVector* bjet2Ptr;
      FourVector* dibjetPtr;
      FourVector* bbggPtr;

  }; 
}


#endif

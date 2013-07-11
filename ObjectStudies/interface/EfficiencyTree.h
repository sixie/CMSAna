#ifndef EFF_DATA_HH
#define EFF_DATA_HH


#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include <cmath>
#include "assert.h"

namespace cmsana
{
  class EfficiencyTree {

    public:

      /// variables
      Float_t                 weight;
      Float_t                 mass;
      Float_t                 pt;
      Float_t                 eta;
      Float_t                 phi;
      Float_t                 rho;
      Int_t                   q;
      UInt_t                  npv;
      UInt_t                  npu;
      Int_t                   matchedPdgId;
      UInt_t                  run;
      UInt_t                  lumi;
      UInt_t                  event;
      Bool_t                  pass;

    public:
      /// this is the main element
      TTree *tree_;
      TFile *f_;
  
      /// hold the names of variables to facilitate things (filled during Init)
      std::vector<std::string> variables_;

      /// default constructor  
      EfficiencyTree()  {};
      /// default destructor
      ~EfficiencyTree(){ 
        if (f_) f_->Close();  
      };
    
      /// initialize varibles and fill list of available variables
      void InitVariables() {
        weight   		       = 0.0;
        mass     		       = -99;
        pt                             = -99;
        eta                            = -99;
        phi                            = -99;
        rho                            = -99;
        q                              = -99;
        npv                            = -99;
        npu                            = -99;
        matchedPdgId                   = 0;
        run                            = 0;
        lumi                           = 0;
        event                          = 0;
        pass                           = false;
      }
    
      /// load a EfficiencyTree
      void LoadTree(const char* file){
        f_ = TFile::Open(file);
        assert(f_);
        tree_ = dynamic_cast<TTree*>(f_->Get("eff"));
        assert(tree_);
      }
    
      /// create a EfficiencyTree
      void CreateTree(){
        tree_ = new TTree("eff","eff");
        f_ = 0;

        //book the branches
        tree_->Branch("weight",&weight,"weight/F");
        tree_->Branch("mass"  ,&mass,  "mass/F");
        tree_->Branch("pt"    ,&pt,    "pt/F");
        tree_->Branch("eta"   ,&eta,   "eta/F");
        tree_->Branch("phi"   ,&phi,   "phi/F");
        tree_->Branch("rho"   ,&rho,   "rho/F");
        tree_->Branch("q"     ,&q,     "q/I");
        tree_->Branch("npv"   ,&npv,   "npv/i");
        tree_->Branch("npu"   ,&npu,   "npu/i");
        tree_->Branch("matchedPdgId"     ,&matchedPdgId,     "matchedPdgId/I");
        tree_->Branch("run"   ,&run,   "run/i");
        tree_->Branch("lumi"  ,&lumi,  "lumi/i");
        tree_->Branch("event" ,&event, "event/i");
        tree_->Branch("pass"  ,&pass,  "pass/O");

      } 

      // initialze a EfficiencyTree
      void InitTree(){
        assert(tree_);
        // don't forget to set pointers to zero before you set address
        // or you will fully appreciate that "ROOT sucks" :)
        InitVariables();
        //Set branch address
        Int_t currentState = gErrorIgnoreLevel;

        tree_->SetBranchAddress("weight",&weight);
        tree_->SetBranchAddress("mass",&mass);
        tree_->SetBranchAddress("pt",&pt);
        tree_->SetBranchAddress("eta",&eta);
        tree_->SetBranchAddress("phi",&phi);
        tree_->SetBranchAddress("rho",&rho);
        tree_->SetBranchAddress("q",&q);
        tree_->SetBranchAddress("npv",&npv);
        tree_->SetBranchAddress("npu",&npu);
        tree_->SetBranchAddress("matchedPdgId",&matchedPdgId);
        tree_->SetBranchAddress("run",&run);
        tree_->SetBranchAddress("lumi",&lumi);
        tree_->SetBranchAddress("event",&event);
        tree_->SetBranchAddress("pass",&pass);
        gErrorIgnoreLevel = currentState;
      }

  }; 

}

#endif

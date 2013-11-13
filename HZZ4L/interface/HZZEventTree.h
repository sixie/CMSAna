#ifndef HZZEventTree_H
#define HZZEventTree_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include <cmath>
#include "assert.h"

class HZZEventTree {

 public:
    
    enum Channel { kFourE     = 0,
                   kFourMu    = 1,
                   kTwoETwoMu = 2,
                   kTwoMuTwoE = 3,                   
    };
    enum HZZEventTreeType { kStandardTree   = 0,
                            kCompactTree    = 1,
    };

  /// variables
  Float_t                 weight;
  UInt_t                  run;
  UInt_t                  lumi;
  UInt_t                  event;
  Float_t                 rho;
  UInt_t                  nvtx;
  Float_t                 met;

  Int_t                   genl1id;
  Float_t                 genl1pt; 
  Float_t                 genl1eta; 
  Float_t                 genl1phi; 
  Int_t                   genl2id;
  Float_t                 genl2pt; 
  Float_t                 genl2eta; 
  Float_t                 genl2phi; 
  Int_t                   genl3id;
  Float_t                 genl3pt; 
  Float_t                 genl3eta; 
  Float_t                 genl3phi; 
  Int_t                   genl4id;
  Float_t                 genl4pt; 
  Float_t                 genl4eta; 
  Float_t                 genl4phi; 

  Int_t                   genchannel;
  Float_t                 genz1mass; 
  Float_t                 genz1pt; 
  Float_t                 genz1eta; 
  Float_t                 genz1rapidity; 
  Float_t                 genz2mass; 
  Float_t                 genz2pt; 
  Float_t                 genz2eta; 
  Float_t                 genz2rapidity; 
  Float_t                 genzzmass; 
  Float_t                 genzzpt; 
  Float_t                 genzzeta; 
  Float_t                 genzzrapidity; 

  Int_t                   l1id;
  Float_t                 l1pt; 
  Float_t                 l1eta; 
  Float_t                 l1phi; 
  Int_t                   l2id;
  Float_t                 l2pt; 
  Float_t                 l2eta; 
  Float_t                 l2phi; 
  Int_t                   l3id;
  Float_t                 l3pt; 
  Float_t                 l3eta; 
  Float_t                 l3phi; 
  Int_t                   l4id;
  Float_t                 l4pt; 
  Float_t                 l4eta; 
  Float_t                 l4phi; 

  Int_t                   channel;
  Double_t                z1mass; 
  Double_t                z1pt; 
  Double_t                z1eta; 
  Double_t                z1rapidity; 
  Double_t                z2mass; 
  Double_t                z2pt; 
  Double_t                z2eta; 
  Double_t                z2rapidity; 
  Double_t                zzmass; 
  Double_t                zzpt; 
  Double_t                zzeta; 
  Double_t                zzrapidity; 

  Double_t                phi0;
  Double_t                theta0;
  Double_t                phi;
  Double_t                theta1;
  Double_t                theta2;


 public:
  /// this is the main element
  TTree *tree_;
  TFile *f_;
  TDirectory *dir_;
  
  /// hold the names of variables to facilitate things (filled during Init)
  std::vector<std::string> variables_;

  /// default constructor  
  HZZEventTree() {};
  /// default destructor
  ~HZZEventTree(){ 
    if (f_) f_->Close();  
  };
    
  /// initialize varibles and fill list of available variables
  void InitVariables() {

    weight                =  0.0;
    run                   =  0.0;
    lumi                  =  0.0;
    event                 =  0.0;
    rho                   =  0.0;
    nvtx                  =  0.0;
    met                   =  0.0;
    genl1id               =  0;
    genl1pt               =  0.0;
    genl1eta              =  0.0;
    genl1phi              =  0.0;
    genl2id               =  0;
    genl2pt               =  0.0;
    genl2eta              =  0.0;
    genl2phi              =  0.0;
    genl3id               =  0;
    genl3pt               =  0.0;
    genl3eta              =  0.0;
    genl3phi              =  0.0;
    genl4id               =  0;
    genl4pt               =  0.0;
    genl4eta              =  0.0;
    genl4phi              =  0.0;
    genchannel            = -1;
    genz1mass             =  0.0;
    genz1pt               =  0.0;
    genz1eta              =  0.0;
    genz1rapidity         =  0.0;
    genz2mass             =  0.0;
    genz2pt               =  0.0;
    genz2eta              =  0.0;
    genz2rapidity         =  0.0;
    genzzmass             =  0.0;
    genzzpt               =  0.0;
    genzzeta              =  0.0;
    genzzrapidity         =  0.0;
    l1id                  =  0;
    l1pt                  =  0.0;
    l1eta                 =  0.0;
    l1phi                 =  0.0;
    l2id                  =  0;
    l2pt                  =  0.0;
    l2eta                 =  0.0;
    l2phi                 =  0.0;
    l3id                  =  0;
    l3pt                  =  0.0;
    l3eta                 =  0.0;
    l3phi                 =  0.0;
    l4id                  =  0;
    l4pt                  =  0.0;
    l4eta                 =  0.0;
    l4phi                 =  0.0;
    channel               =  0.0;
    z1mass                =  0.0;
    z1pt                  =  0.0;
    z1eta                 =  0.0;
    z1rapidity            =  0.0;
    z2mass                =  0.0;
    z2pt                  =  0.0;
    z2eta                 =  0.0;
    z2rapidity            =  0.0;
    zzmass                =  0.0;
    zzpt                  =  0.0;
    zzeta                 =  0.0;
    zzrapidity            =  0.0;
    phi0                  =  0.0;
    theta0                =  0.0;
    phi                   =  0.0;
    theta1                =  0.0;
    theta2                =  0.0;
  }
    
  /// load a HZZEventTree
  void LoadTree(const char* file){
    f_ = TFile::Open(file);
    assert(f_);
    tree_ = dynamic_cast<TTree*>(f_->Get("zz4lTree"));
    assert(tree_);
  }
    
  /// create a HZZEventTree
  void CreateTree( int type = HZZEventTree::kStandardTree ){
    tree_ = new TTree("zz4lTree","zz4lTree");
    f_ = 0;

    //book the branches
    if (type == HZZEventTree::kCompactTree) {
      tree_->Branch("weight", &weight, "weight/F");
      tree_->Branch("phi0", &phi0, "phi0/F");
      tree_->Branch("theta0", &theta0, "theta0/F");
      tree_->Branch("phi", &phi, "phi/F");
      tree_->Branch("theta1", &theta1, "theta1/F");
      tree_->Branch("theta2", &theta2, "theta2/F");
      tree_->Branch("z1mass",&z1mass,"z1mass/F");
      tree_->Branch("z2mass",&z2mass,"z2mass/F");
      tree_->Branch("zzmass",&zzmass,"zzmass/F");
    }
    
    if (type == HZZEventTree::kStandardTree) {
      tree_->Branch("weight", &weight, "weight/F");
      tree_->Branch("run",&run,"run/i");     
      tree_->Branch("lumi",&lumi, "lumi/i");
      tree_->Branch("event",&event, "event/i");     
      tree_->Branch("rho",&rho,"rho/F");
      tree_->Branch("nvtx",&nvtx,"nvtx/i");
      tree_->Branch("met",&met,"met/F");    
      tree_->Branch("genl1id",&genl1id ,"genl1id/I");         
      tree_->Branch("genl1pt",&genl1pt ,"genl1pt/F");         
      tree_->Branch("genl1eta",&genl1eta ,"genl1eta/F");          
      tree_->Branch("genl1phi",&genl1phi ,"genl1phi/F");         
      tree_->Branch("genl2id",&genl2id ,"genl2id/I");         
      tree_->Branch("genl2pt",&genl2pt ,"genl2pt/F");         
      tree_->Branch("genl2eta",&genl2eta ,"genl2eta/F");          
      tree_->Branch("genl2phi",&genl2phi ,"genl2phi/F");         
      tree_->Branch("genl3id",&genl3id ,"genl3id/I");         
      tree_->Branch("genl3pt",&genl3pt ,"genl3pt/F");         
      tree_->Branch("genl3eta",&genl3eta ,"genl3eta/F");          
      tree_->Branch("genl3phi",&genl3phi ,"genl3phi/F");         
      tree_->Branch("genl4id",&genl4id ,"genl4id/I");         
      tree_->Branch("genl4pt",&genl4pt ,"genl4pt/F");         
      tree_->Branch("genl4eta",&genl4eta ,"genl4eta/F");          
      tree_->Branch("genl4phi",&genl4phi ,"genl4phi/F");         
      tree_->Branch("genchannel",&genchannel,"genchannel/I");         
      tree_->Branch("genz1mass",&genz1mass,"genz1mass/F");           
      tree_->Branch("genz1pt",&genz1pt,"genz1pt/F");           
      tree_->Branch("genz1eta",&genz1eta,"genz1eta/F");           
      tree_->Branch("genz1rapidity",&genz1rapidity,"genz1rapidity/F");           
      tree_->Branch("genz2mass",&genz2mass,"genz2mass/F");           
      tree_->Branch("genz2pt",&genz2pt,"genz2pt/F");           
      tree_->Branch("genz2eta",&genz2eta,"genz2eta/F");           
      tree_->Branch("genz2rapidity",&genz2rapidity,"genz2rapidity/F");           
      tree_->Branch("genzzmass",&genzzmass,"genzzmass/F");           
      tree_->Branch("genzzpt",&genzzpt,"genzzpt/F");           
      tree_->Branch("genzzeta",&genzzeta,"genzzeta/F");           
      tree_->Branch("genzzrapidity",&genzzrapidity,"genzzrapidity/F");           

      tree_->Branch("l1id",&l1id ,"l1id/I");         
      tree_->Branch("l1pt",&l1pt ,"l1pt/F");         
      tree_->Branch("l1eta",&l1eta ,"l1eta/F");          
      tree_->Branch("l1phi",&l1phi ,"l1phi/F");         
      tree_->Branch("l2id",&l2id ,"l2id/I");         
      tree_->Branch("l2pt",&l2pt ,"l2pt/F");         
      tree_->Branch("l2eta",&l2eta ,"l2eta/F");          
      tree_->Branch("l2phi",&l2phi ,"l2phi/F");         
      tree_->Branch("l3id",&l3id ,"l3id/I");         
      tree_->Branch("l3pt",&l3pt ,"l3pt/F");         
      tree_->Branch("l3eta",&l3eta ,"l3eta/F");          
      tree_->Branch("l3phi",&l3phi ,"l3phi/F");         
      tree_->Branch("l4id",&l4id ,"l4id/I");         
      tree_->Branch("l4pt",&l4pt ,"l4pt/F");         
      tree_->Branch("l4eta",&l4eta ,"l4eta/F");          
      tree_->Branch("l4phi",&l4phi ,"l4phi/F");         
      tree_->Branch("channel",&channel,"channel/I");         
      tree_->Branch("z1mass",&z1mass,"z1mass/D");           
      tree_->Branch("z1pt",&z1pt,"z1pt/D");           
      tree_->Branch("z1eta",&z1eta,"z1eta/D");           
      tree_->Branch("z1rapidity",&z1rapidity,"z1rapidity/D");           
      tree_->Branch("z2mass",&z2mass,"z2mass/D");           
      tree_->Branch("z2pt",&z2pt,"z2pt/D");           
      tree_->Branch("z2eta",&z2eta,"z2eta/D");           
      tree_->Branch("z2rapidity",&z2rapidity,"z2rapidity/D");           
      tree_->Branch("zzmass",&zzmass,"zzmass/D");           
      tree_->Branch("zzpt",&zzpt,"zzpt/D");           
      tree_->Branch("zzeta",&zzeta,"zzeta/D");           
      tree_->Branch("zzrapidity",&zzrapidity,"zzrapidity/D");           
    }
  } 

    // initialze a HZZEventTree
    void InitTree(int type = HZZEventTree::kStandardTree ){
      assert(tree_);
      // don't forget to set pointers to zero before you set address
      // or you will fully appreciate that "ROOT sucks" :)
      InitVariables();
      //Set branch address
      Int_t currentState = gErrorIgnoreLevel;
      // gErrorIgnoreLevel = kError;
      gErrorIgnoreLevel = kBreak;
    
      if (type == HZZEventTree::kCompactTree) {
        tree_->SetBranchAddress("weight", &weight);
        tree_->SetBranchAddress("phi0", &phi0);
        tree_->SetBranchAddress("theta0", &theta0);
        tree_->SetBranchAddress("phi", &phi);
        tree_->SetBranchAddress("theta1", &theta1);
        tree_->SetBranchAddress("theta2", &theta2);
        tree_->SetBranchAddress("z1mass",&z1mass);
        tree_->SetBranchAddress("z2mass",&z2mass);
        tree_->SetBranchAddress("zzmass",&zzmass);
      }
    
      if (type == HZZEventTree::kStandardTree) {

        tree_->SetBranchAddress("weight", &weight);
        tree_->SetBranchAddress("run",&run);
        tree_->SetBranchAddress("lumi",&lumi);
        tree_->SetBranchAddress("event",&event);     
        tree_->SetBranchAddress("rho",&rho);
        tree_->SetBranchAddress("nvtx",&nvtx);
        tree_->SetBranchAddress("met",&met);
        tree_->SetBranchAddress("genl1id",&genl1id);
        tree_->SetBranchAddress("genl1pt",&genl1pt);
        tree_->SetBranchAddress("genl1eta",&genl1eta);
        tree_->SetBranchAddress("genl1phi",&genl1phi);
        tree_->SetBranchAddress("genl2id",&genl2id);
        tree_->SetBranchAddress("genl2pt",&genl2pt);
        tree_->SetBranchAddress("genl2eta",&genl2eta);
        tree_->SetBranchAddress("genl2phi",&genl2phi);
        tree_->SetBranchAddress("genl3id",&genl3id);
        tree_->SetBranchAddress("genl3pt",&genl3pt);
        tree_->SetBranchAddress("genl3eta",&genl3eta);
        tree_->SetBranchAddress("genl3phi",&genl3phi);
        tree_->SetBranchAddress("genl4id",&genl4id);
        tree_->SetBranchAddress("genl4pt",&genl4pt);
        tree_->SetBranchAddress("genl4eta",&genl4eta);
        tree_->SetBranchAddress("genl4phi",&genl4phi);
        tree_->SetBranchAddress("genchannel",&genchannel);
        tree_->SetBranchAddress("genz1mass",&genz1mass);
        tree_->SetBranchAddress("genz1pt",&genz1pt);
        tree_->SetBranchAddress("genz1eta",&genz1eta);
        tree_->SetBranchAddress("genz1rapidity",&genz1rapidity);
        tree_->SetBranchAddress("genz2mass",&genz2mass);
        tree_->SetBranchAddress("genz2pt",&genz2pt);
        tree_->SetBranchAddress("genz2eta",&genz2eta);
        tree_->SetBranchAddress("genz2rapidity",&genz2rapidity);
        tree_->SetBranchAddress("genzzmass",&genzzmass);
        tree_->SetBranchAddress("genzzpt",&genzzpt);
        tree_->SetBranchAddress("genzzeta",&genzzeta);
        tree_->SetBranchAddress("genzzrapidity",&genzzrapidity);

        tree_->SetBranchAddress("l1id",&l1id);
        tree_->SetBranchAddress("l1pt",&l1pt);
        tree_->SetBranchAddress("l1eta",&l1eta);
        tree_->SetBranchAddress("l1phi",&l1phi);
        tree_->SetBranchAddress("l2id",&l2id);
        tree_->SetBranchAddress("l2pt",&l2pt);
        tree_->SetBranchAddress("l2eta",&l2eta);
        tree_->SetBranchAddress("l2phi",&l2phi);
        tree_->SetBranchAddress("l3id",&l3id);
        tree_->SetBranchAddress("l3pt",&l3pt);
        tree_->SetBranchAddress("l3eta",&l3eta);
        tree_->SetBranchAddress("l3phi",&l3phi);
        tree_->SetBranchAddress("l4id",&l4id);
        tree_->SetBranchAddress("l4pt",&l4pt);
        tree_->SetBranchAddress("l4eta",&l4eta);
        tree_->SetBranchAddress("l4phi",&l4phi);
        tree_->SetBranchAddress("channel",&channel);
        tree_->SetBranchAddress("z1mass",&z1mass);
        tree_->SetBranchAddress("z1pt",&z1pt);
        tree_->SetBranchAddress("z1eta",&z1eta);
        tree_->SetBranchAddress("z1rapidity",&z1rapidity);
        tree_->SetBranchAddress("z2mass",&z2mass);
        tree_->SetBranchAddress("z2pt",&z2pt);
        tree_->SetBranchAddress("z2eta",&z2eta);
        tree_->SetBranchAddress("z2rapidity",&z2rapidity);
        tree_->SetBranchAddress("zzmass",&zzmass);
        tree_->SetBranchAddress("zzpt",&zzpt);
        tree_->SetBranchAddress("zzeta",&zzeta);
        tree_->SetBranchAddress("zzrapidity",&zzrapidity);
      }

      gErrorIgnoreLevel = currentState;
    }

}; 



#endif


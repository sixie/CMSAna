
//================================================================================================
//
// HZZ4l selection macro
//
//  * plots distributions associated with selected events
//  * prints list of selected events from data
//  * outputs ROOT files of events passing selection for each sample, 
//    which can be processed by plotSelect.C
//
//________________________________________________________________________________________________


#include "CMSAna/HZZ4L/ZJets/YYDataFormat/rootheader.h"

TChain *fChain; 


#include "CMSAna/HZZ4L/ZJets/YYDataFormat/RecoAnalyzer.h"
#include "CMSAna/HZZ4L/ZJets/YYDataFormat/setbranchaddress.cc"
#include "CMSAna/HZZ4L/ZJets/YYDataFormat/utils.cc"


void MakeHZZEventNtupleForPythiaZJetsSampleExample(string label = "datasetname"){


  TH1F *mass4L = new TH1F ("mass4L", " ; mass4L [GeV/c^{2}]; Number of Events", 200, 0, 200);
  
  //***************************************************
  //Set up Output ntuple
  //***************************************************
  TFile *fOutputFile = new TFile(("HZZ4lEvent.PythiaZJets."+label+".root").c_str(), "RECREATE");
  TTree *fOutputTree = new TTree("HZZEventTree","HZZEventTree");

  Int_t                   fLep1Type;
  Float_t                 fLep1Pt; 
  Float_t                 fLep1Eta; 
  Float_t                 fLep1Phi; 
  Int_t                   fLep2Type;
  Float_t                 fLep2Pt; 
  Float_t                 fLep2Eta; 
  Float_t                 fLep2Phi; 
  Int_t                   fLep3Type;
  Float_t                 fLep3Pt; 
  Float_t                 fLep3Eta; 
  Float_t                 fLep3Phi; 
  Int_t                   fLep4Type;
  Float_t                 fLep4Pt; 
  Float_t                 fLep4Eta; 
  Float_t                 fLep4Phi; 

  fOutputTree->Branch("l1pdgId",&fLep1Type,"l1pdgId/I");         
  fOutputTree->Branch("l1pt",&fLep1Pt,"l1pt/F");           
  fOutputTree->Branch("l1eta",&fLep1Eta,"l1eta/F");           
  fOutputTree->Branch("l1phi",&fLep1Phi,"l1phi/F");           
  fOutputTree->Branch("l2pdgId",&fLep2Type,"l2pdgId/I");         
  fOutputTree->Branch("l2pt",&fLep2Pt,"l2pt/F");           
  fOutputTree->Branch("l2eta",&fLep2Eta,"l2eta/F");           
  fOutputTree->Branch("l2phi",&fLep2Phi,"l2phi/F");           
  fOutputTree->Branch("l3pdgId",&fLep3Type,"l3pdgId/I");         
  fOutputTree->Branch("l3pt",&fLep3Pt,"l3pt/F");           
  fOutputTree->Branch("l3eta",&fLep3Eta,"l3eta/F");           
  fOutputTree->Branch("l3phi",&fLep3Phi,"l3phi/F");           
  fOutputTree->Branch("l4pdgId",&fLep4Type,"l4pdgId/I");         
  fOutputTree->Branch("l4pt",&fLep4Pt,"l4pt/F");           
  fOutputTree->Branch("l4eta",&fLep4Eta,"l4eta/F");           
  fOutputTree->Branch("l4phi",&fLep4Phi,"l4phi/F");           
  


  int NEventsPass = 0;

  fChain = new TChain("Analysis");
  
  fChain->Add("root://eoscms//eos/cms//store/group/phys_higgs/cmshzz4l/sixie/pythiazjets/analysis_9988.root");
  setbranchaddress();

  cout<<"br set " <<endl; 
  
  int totalEntries = fChain->GetEntries();
  cout<<" totalEntries " << totalEntries <<endl; 


  //loop over objects
  for (int entry=0; entry< totalEntries ; entry ++) {
    if (entry % 100000 == 0) cout << "Entry : " << entry << "\n";
    fChain->GetEntry(entry);
    
    int nGenEleFromZ = 0;
    int nGenMuFromZ = 0;

    TLorentzVector *ZLepton1;
    TLorentzVector *ZLepton2;
    int ZLepton1PdgId = 0;
    int ZLepton2PdgId = 0;

    for(int j=0; j< nGenEle; j++){
      if( statusGenEle[j] !=3 ) continue; 
      if( pidmomGenEle[j] != 23) continue; 

      if (nGenEleFromZ == 0) {
        ZLepton1 = new TLorentzVector();
        ZLepton1->SetPtEtaPhiM( ptGenEle[j] , etaGenEle[j], phiGenEle[j], 0.51099892e-3);
        if (chaGenEle[j] == -1) ZLepton1PdgId = 11;
        else ZLepton1PdgId = -11;
      }
      if (nGenEleFromZ == 1) {
        ZLepton2 = new TLorentzVector();
        ZLepton2->SetPtEtaPhiM( ptGenEle[j] , etaGenEle[j], phiGenEle[j], 0.51099892e-3);
        if (chaGenEle[j] == -1) ZLepton2PdgId = 11;
        else ZLepton2PdgId = -11;
      }
      nGenEleFromZ++;
    }    
  
    for(int j=0; j< nGenMu; j++){
      if( statusGenMu[j] !=3 ) continue; 
      if( pidmomGenMu[j] != 23) continue; 

      if (nGenMuFromZ == 0) {
        ZLepton1 = new TLorentzVector();
        ZLepton1->SetPtEtaPhiM( ptGenMu[j] , etaGenMu[j], phiGenMu[j], 105.658369e-3);
        if (chaGenMu[j] == -1) ZLepton1PdgId = 13;
        else ZLepton1PdgId = -13;
      }
      if (nGenMuFromZ == 1) {
        ZLepton2 = new TLorentzVector();
        ZLepton2->SetPtEtaPhiM( ptGenMu[j] , etaGenMu[j], phiGenMu[j], 105.658369e-3);
        if (chaGenMu[j] == -1) ZLepton2PdgId = 13;
        else ZLepton2PdgId = -13;
      }
      nGenMuFromZ++;
    }
     
    if (
      !(
        (nGenEleFromZ == 2 && nGenMuFromZ == 0) ||
        (nGenEleFromZ == 0 && nGenMuFromZ == 2) 
        )
      ) {
      cout << "Warning: nGenEleFromZ == " << nGenEleFromZ << " and nGenMuFromZ == " << nGenMuFromZ << "\n";
    }
   
    assert(ZLepton1);
    assert(ZLepton2);

    bool passMassCut = false;

    if ( (*ZLepton1 + *ZLepton2).M() < 2) {

      cout << nGenEleFromZ << " " << nGenMuFromZ << "\n";
      cout << "ZLepton1 : " << ZLepton1->Pt() << " " << ZLepton1->Eta() << " " << ZLepton1->Phi() << "\n";
      cout << "ZLepton2 : " << ZLepton2->Pt() << " " << ZLepton2->Eta() << " " << ZLepton2->Phi() << "\n";

      for(int j=0; j< nGenEle; j++){
        if( statusGenEle[j] !=3 ) continue; 
        if( pidmomGenEle[j] != 23) continue; 
        cout << "Ele: " << ptGenEle[j] << " " << etaGenEle[j] << " " << phiGenEle[j] << "\n";
      }    
      
      for(int j=0; j< nGenMu; j++){
        if( statusGenMu[j] !=3 ) continue; 
        if( pidmomGenMu[j] != 23) continue; 
        cout << "Mu: " << ptGenMu[j] << " " << etaGenMu[j] << " " << phiGenMu[j] << "\n";
      }
      
      
    }


    for(int j=0; j< nGenJets; j++){
      if (ptGenJet[j] < 5) continue;
      for(int k=j+1; k< nGenJets; k++){
        if (ptGenJet[k] < 5) continue;

        TLorentzVector *fakeLepton1 = new TLorentzVector();
        TLorentzVector *fakeLepton2 = new TLorentzVector();
        fakeLepton1->SetPtEtaPhiM( ptGenJet[j] , etaGenJet[j], phiGenJet[j], 0);
        fakeLepton2->SetPtEtaPhiM( ptGenJet[k] , etaGenJet[k], phiGenJet[k], 0);

        float mass = (*ZLepton1 + *ZLepton2 + *fakeLepton1 + *fakeLepton2).M();
        mass4L->Fill(mass);
//         if (mass < 100) {
//           mass4L->Fill((*ZLepton1+*ZLepton2).M());
//         }

        if (mass > 120 && mass < 130) {
          passMassCut = true;
        }


        //Fill event ntuple
        if (passMassCut) {
          fLep1Type = ZLepton1PdgId;
          fLep1Pt = ZLepton1->Pt();
          fLep1Eta = ZLepton1->Eta();
          fLep1Phi = ZLepton1->Phi();
          fLep2Type = ZLepton2PdgId;
          fLep2Pt = ZLepton2->Pt();
          fLep2Eta = ZLepton2->Eta();
          fLep2Phi = ZLepton2->Phi();
          fLep3Type = 13;
          fLep3Pt = fakeLepton1->Pt();
          fLep3Eta = fakeLepton1->Eta();
          fLep3Phi = fakeLepton1->Phi();
          fLep4Type = -13;
          fLep4Pt = fakeLepton2->Pt();
          fLep4Eta = fakeLepton2->Eta();
          fLep4Phi = fakeLepton2->Phi();
          fOutputTree->Fill();
        }

        delete fakeLepton1;
        delete fakeLepton2;

      }
    }

    if (passMassCut) NEventsPass++;

    delete ZLepton1;
    delete ZLepton2;

  } //loop over events


  cout << "NEventsPass: " << NEventsPass << " / Total ( " << totalEntries << " )\n";
  mass4L->Draw();
  
  fOutputFile->Write();
  fOutputFile->Close();
  delete fOutputFile;

}

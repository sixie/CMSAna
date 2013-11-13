
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

#define PI 3.1415926
#include "CMSAna/HZZ4L/ZJets/YYDataFormat/rootheader.h"
#include <TRandom3.h>                

TChain *fChain; 


#include "CMSAna/HZZ4L/ZJets/YYDataFormat/RecoAnalyzer.h"
#include "CMSAna/HZZ4L/ZJets/YYDataFormat/setbranchaddress.cc"
#include "CMSAna/HZZ4L/ZJets/YYDataFormat/utils.cc"

#include "CMSAna/HZZ4L/HLL/LeptonResponseMap.hh"
#include "CMSAna/HZZ4L/HLL/LeptonResolutionData_LegacyPaper/ElectronEfficiencyMap.h"
#include "CMSAna/HZZ4L/HLL/LeptonResolutionData_LegacyPaper/MuonEfficiencyMap.h"
#include "CMSAna/HZZ4L/HLL/LeptonResolutionData_LegacyPaper/FakeElectronEfficiencyMap.h"
#include "CMSAna/HZZ4L/HLL/LeptonResolutionData_LegacyPaper/FakeMuonEfficiencyMap.h"
#include "CMSAna/HZZ4L/HLL/LeptonResolutionData_LegacyPaper/FakeElectronResponseMap.h"
#include "CMSAna/HZZ4L/HLL/LeptonResolutionData_LegacyPaper/FakeMuonResponseMap.h"

// output data structs
#include "CMSAna/HZZ4L/interface/HZZEventTree.h"
#include "CMSAna/Utils/CommonTools.hh"

void MakeHZZEventNtupleForPythiaZJetsSample(string label = "datasetname"){

  gBenchmark->Start("test");
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  TRandom3 *MyRandom = new TRandom3( time(NULL) );
  TFile *LeptonResponseFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/syncForRun1LegacyPaper/CMSSW_5_3_9_patch3/src/CMSAna/HZZ4L/HLL/LeptonResolutionData_LegacyPaper/PtResolutionModel_combined.LegacyPaper.root","READ");

  TH1F *mass4L = new TH1F ("mass4L", " ; mass4L [GeV/c^{2}]; Number of Events", 200, 0, 200);
  
  //***************************************************
  //Set up Output ntuple
  //***************************************************
  TFile *outputFile = new TFile(("HZZ4lEvent.PythiaZJets."+label+".root").c_str(), "RECREATE");
  HZZEventTree *hzzEventTree = new HZZEventTree;
  hzzEventTree->CreateTree(HZZEventTree::kStandardTree);
  hzzEventTree->tree_->SetAutoFlush(0);

  int NEventsPass = 0;
  int NEvents = 0;


  fChain = new TChain("Analysis");
  
  fChain->Add("root://eoscms//eos/cms//store/group/phys_higgs/cmshzz4l/sixie/pythiazjets/analysis_9988.root");
  setbranchaddress();

  cout<<"br set " <<endl; 
  
  int totalEntries = fChain->GetEntries();
  cout<<" totalEntries " << totalEntries <<endl; 


  //loop over objects
  for (int entry=0; entry< totalEntries ; entry ++) {
    if (entry % 10 == 0) cout << "Entry : " << entry << "\n";
    fChain->GetEntry(entry);

    int NRepeat = 20;
    for (int repeat = 0; repeat < NRepeat; ++repeat) {

    
      int nGenEleFromZ = 0;
      int nGenMuFromZ = 0;

      TLorentzVector *ZGenLepton1;
      TLorentzVector *ZGenLepton2;
      TLorentzVector *ZLepton1;
      TLorentzVector *ZLepton2;
      int ZLepton1PdgId = 0;
      int ZLepton2PdgId = 0;
      double ZLepton1Eff = 0;
      double ZLepton2Eff = 0;
      double ZLepton1Pt = 0;
      double ZLepton2Pt = 0;
      double ZLepton1GenPt = 0;
      double ZLepton2GenPt = 0;

      for(int j=0; j< nGenEle; j++){
        if( statusGenEle[j] !=3 ) continue; 
        if( pidmomGenEle[j] != 23) continue; 

        if (nGenEleFromZ == 0) {
          ZLepton1GenPt = ptGenEle[j];
          ZLepton1Eff = GetElectronEfficiencyPtEta(ptGenEle[j], fabs(etaGenEle[j]));
          ZLepton1Pt = GenerateLeptonPtFromGaussianModel(LeptonResponseFile, MyRandom, 11, ptGenEle[j], fabs(etaGenEle[j]));
          ZLepton1 = new TLorentzVector();
          ZLepton1->SetPtEtaPhiM( ZLepton1Pt , etaGenEle[j], phiGenEle[j], 0.51099892e-3);
          ZGenLepton1 = new TLorentzVector();
          ZGenLepton1->SetPtEtaPhiM( ptGenEle[j] , etaGenEle[j], phiGenEle[j], 0.51099892e-3);
          if (chaGenEle[j] == -1) ZLepton1PdgId = 11;
          else ZLepton1PdgId = -11;
        }
        if (nGenEleFromZ == 1) {
          ZLepton2GenPt = ptGenEle[j];
          ZLepton2Eff = GetElectronEfficiencyPtEta(ptGenEle[j], fabs(etaGenEle[j]));
          ZLepton2Pt = GenerateLeptonPtFromGaussianModel(LeptonResponseFile, MyRandom, 11, ptGenEle[j], fabs(etaGenEle[j]));
          ZLepton2 = new TLorentzVector();
          ZLepton2->SetPtEtaPhiM( ZLepton2Pt , etaGenEle[j], phiGenEle[j], 0.51099892e-3);
          ZGenLepton2 = new TLorentzVector();
          ZGenLepton2->SetPtEtaPhiM( ptGenEle[j] , etaGenEle[j], phiGenEle[j], 0.51099892e-3);
          if (chaGenEle[j] == -1) ZLepton2PdgId = 11;
          else ZLepton2PdgId = -11;
        }
        nGenEleFromZ++;
      }    

  
      for(int j=0; j< nGenMu; j++){
        if( statusGenMu[j] !=3 ) continue; 
        if( pidmomGenMu[j] != 23) continue; 

        if (nGenMuFromZ == 0) {
          ZLepton1GenPt = ptGenMu[j];
          ZLepton1Eff = GetElectronEfficiencyPtEta(ptGenMu[j], fabs(etaGenMu[j]));
          ZLepton1Pt = GenerateLeptonPtFromGaussianModel(LeptonResponseFile, MyRandom, 13, ptGenMu[j], fabs(etaGenMu[j]));
          ZLepton1 = new TLorentzVector();
          ZLepton1->SetPtEtaPhiM( ZLepton1Pt , etaGenMu[j], phiGenMu[j], 105.658369e-3);
          ZGenLepton1 = new TLorentzVector();
          ZGenLepton1->SetPtEtaPhiM( ptGenMu[j] , etaGenMu[j], phiGenMu[j], 105.658369e-3);
          if (chaGenMu[j] == -1) ZLepton1PdgId = 13;
          else ZLepton1PdgId = -13;
        }
        if (nGenMuFromZ == 1) {
          ZLepton2GenPt = ptGenMu[j];
          ZLepton2Eff = GetElectronEfficiencyPtEta(ptGenMu[j], fabs(etaGenMu[j]));
          ZLepton2Pt = GenerateLeptonPtFromGaussianModel(LeptonResponseFile, MyRandom, 13, ptGenMu[j], fabs(etaGenMu[j]));
          ZLepton2 = new TLorentzVector();
          ZLepton2->SetPtEtaPhiM( ZLepton2Pt , etaGenMu[j], phiGenMu[j], 105.658369e-3);
          ZGenLepton2 = new TLorentzVector();
          ZGenLepton2->SetPtEtaPhiM( ptGenMu[j] , etaGenMu[j], phiGenMu[j], 105.658369e-3);
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

        //if genjet is too close to leptons, skip them
        if (cmsana::deltaR(etaGenJet[j],phiGenJet[j],ZLepton1->Eta(),ZLepton1->Phi()) < 0.4) continue;
        if (cmsana::deltaR(etaGenJet[j],phiGenJet[j],ZLepton2->Eta(),ZLepton2->Phi()) < 0.4) continue;

        for(int k=j+1; k< nGenJets; k++){
          if (ptGenJet[k] < 5) continue;

          //if genjet is too close to leptons or other genjet, skip them
          if (cmsana::deltaR(etaGenJet[k],phiGenJet[k],ZLepton1->Eta(),ZLepton1->Phi()) < 0.4) continue;
          if (cmsana::deltaR(etaGenJet[k],phiGenJet[k],ZLepton2->Eta(),ZLepton2->Phi()) < 0.4) continue;
          if (cmsana::deltaR(etaGenJet[k],phiGenJet[k],etaGenJet[j],phiGenJet[j]) < 0.5) continue;

          double muonFakerate_j = GetFakeMuonEfficiencyPtEta(ptGenJet[j], fabs(etaGenJet[j]));
          double muonFakerate_k = GetFakeMuonEfficiencyPtEta(ptGenJet[k], fabs(etaGenJet[k]));
          double muPt_j = ptGenJet[j]*(1.0+MyRandom->Gaus(GetMuonResponseMeanPtEta(ptGenJet[j], fabs(etaGenJet[j])),
                                                          GetMuonResponseSigmaPtEta(ptGenJet[j], fabs(etaGenJet[j]))));
          double muEta_j = etaGenJet[j] + MyRandom->Gaus(0,0.035);
          double muPhi_j = phiGenJet[j] + MyRandom->Gaus(0,0.035);
          if (muPhi_j < -1*PI) muPhi_j = muPhi_j + 2*PI;
          if (muPhi_j > PI) muPhi_j = muPhi_j - 2*PI;

          double muPt_k = ptGenJet[k]*(1.0+MyRandom->Gaus(GetMuonResponseMeanPtEta(ptGenJet[k], fabs(etaGenJet[k])),
                                                          GetMuonResponseSigmaPtEta(ptGenJet[k], fabs(etaGenJet[k]))));
          double muEta_k = etaGenJet[k] + MyRandom->Gaus(0,0.035);
          double muPhi_k = phiGenJet[k] + MyRandom->Gaus(0,0.035);
          if (muPhi_k < -1*PI) muPhi_k = muPhi_k + 2*PI;
          if (muPhi_k > PI) muPhi_k = muPhi_k - 2*PI;


          double eleFakerate_j = GetFakeElectronEfficiencyPtEta(ptGenJet[j], fabs(etaGenJet[j]));
          double eleFakerate_k = GetFakeElectronEfficiencyPtEta(ptGenJet[k], fabs(etaGenJet[k]));
          double elePt_j = ptGenJet[j]*(1.0+MyRandom->Gaus(GetElectronResponseMeanPtEta(ptGenJet[j], fabs(etaGenJet[j])),
                                                           GetElectronResponseSigmaPtEta(ptGenJet[j], fabs(etaGenJet[j]))));
          double eleEta_j = etaGenJet[j] + MyRandom->Gaus(0,0.02);
          double elePhi_j = phiGenJet[j] + MyRandom->Gaus(0,0.025);
          if (elePhi_j < -1*PI) elePhi_j = elePhi_j + 2*PI;
          if (elePhi_j > PI) elePhi_j = elePhi_j - 2*PI;

          double elePt_k = ptGenJet[k]*(1.0+MyRandom->Gaus(GetElectronResponseMeanPtEta(ptGenJet[k], fabs(etaGenJet[k])),
                                                           GetElectronResponseSigmaPtEta(ptGenJet[k], fabs(etaGenJet[k]))));
          double eleEta_k = etaGenJet[k] + MyRandom->Gaus(0,0.02);
          double elePhi_k = phiGenJet[k] + MyRandom->Gaus(0,0.025);
          if (elePhi_k < -1*PI) elePhi_k = elePhi_k + 2*PI;
          if (elePhi_k > PI) elePhi_k = elePhi_k - 2*PI;


          double mass = 0;
          bool passMassCut = false;

          //*********************************************************************************************************
          //Fill Fake Muon + Fake Muon Case
          //*********************************************************************************************************
          passMassCut = false;
          TLorentzVector *fakeMu1 = new TLorentzVector();
          TLorentzVector *fakeMu2 = new TLorentzVector();
          fakeMu1->SetPtEtaPhiM( muPt_j , muEta_j, muPhi_j, 0);
          fakeMu2->SetPtEtaPhiM( muPt_k , muEta_k, muPhi_k, 0);
          TLorentzVector *GenFakeMu1 = new TLorentzVector();
          TLorentzVector *GenFakeMu2 = new TLorentzVector();
          GenFakeMu1->SetPtEtaPhiM( ptGenJet[j], etaGenJet[j], phiGenJet[j], 0);
          GenFakeMu2->SetPtEtaPhiM( ptGenJet[k] , etaGenJet[k], phiGenJet[k], 0);

          mass = (*ZLepton1 + *ZLepton2 + *fakeMu1 + *fakeMu2).M();
          mass4L->Fill(mass);

          if (mass > 100 && mass < 150) {
            passMassCut = true;
          }

          //Fill event ntuple
          if (passMassCut) {

            double leadingLepPt = ZLepton1Pt;
            double subleadingLepPt = ZLepton2Pt;
            if (ZLepton2Pt > ZLepton1Pt) {
              leadingLepPt = ZLepton2Pt;
              subleadingLepPt = ZLepton1Pt;
            }
            if (muPt_j > leadingLepPt) {
              subleadingLepPt = leadingLepPt;
              leadingLepPt = muPt_j;      
            } else if (muPt_j > subleadingLepPt) {
              subleadingLepPt = muPt_j;
            }
            if (muPt_k > leadingLepPt) {
              subleadingLepPt = leadingLepPt;
              leadingLepPt = muPt_k;      
            } else if (muPt_k > subleadingLepPt) {
              subleadingLepPt = muPt_k;
            }

            if (ZLepton1Pt > (abs(ZLepton1PdgId) == 11 ? 7 : 5)
                && ZLepton2Pt > (abs(ZLepton2PdgId) == 11 ? 7 : 5)
                && muPt_j >  5
                && muPt_k >  5
                && fabs(ZLepton1->Eta()) < (abs(ZLepton1PdgId) == 11 ? 2.5 : 2.4)
                && fabs(ZLepton2->Eta()) < (abs(ZLepton2PdgId) == 11 ? 2.5 : 2.4)
                && fabs(muEta_j) < 2.4
                && fabs(muEta_k) < 2.4
                && leadingLepPt > 20
                && subleadingLepPt > 10
              ) {

              //*******************************************
              // Fill HZZ Event
              //*******************************************
              hzzEventTree->weight = ZLepton1Eff*ZLepton2Eff*muonFakerate_j*muonFakerate_k;
              hzzEventTree->run = 0;
              hzzEventTree->lumi = 0;
              hzzEventTree->event = NEvents;
              hzzEventTree->rho = 0;
              hzzEventTree->nvtx = 0;
              hzzEventTree->met = 0;

              hzzEventTree->genl1id = ZLepton1PdgId;
              hzzEventTree->genl1pt = ZLepton1GenPt;
              hzzEventTree->genl1eta = ZLepton1->Eta();
              hzzEventTree->genl1phi = ZLepton1->Phi();
              hzzEventTree->genl2id = ZLepton2PdgId;
              hzzEventTree->genl2pt = ZLepton2GenPt;
              hzzEventTree->genl2eta = ZLepton2->Eta();
              hzzEventTree->genl2phi = ZLepton2->Phi();
              hzzEventTree->genl3id = 13;
              hzzEventTree->genl3pt = ptGenJet[j];
              hzzEventTree->genl3eta = etaGenJet[j];
              hzzEventTree->genl3phi = phiGenJet[j];
              hzzEventTree->genl4id = -13;
              hzzEventTree->genl4pt = ptGenJet[k];
              hzzEventTree->genl4eta = etaGenJet[k];
              hzzEventTree->genl4phi = phiGenJet[k];

              hzzEventTree->genchannel = (abs(ZLepton1PdgId) == 11) ? HZZEventTree::kTwoETwoMu : HZZEventTree::kFourMu;
              hzzEventTree->genz1mass = ((*ZGenLepton1)+(*ZGenLepton2)).M();
              hzzEventTree->genz1pt = ((*ZGenLepton1)+(*ZGenLepton2)).Pt();
              hzzEventTree->genz1eta = ((*ZGenLepton1)+(*ZGenLepton2)).Eta();
              hzzEventTree->genz1rapidity = ((*ZGenLepton1)+(*ZGenLepton2)).Rapidity();
              hzzEventTree->genz2mass = ((*GenFakeMu1)+(*GenFakeMu2)).M();
              hzzEventTree->genz2pt = ((*GenFakeMu1)+(*GenFakeMu2)).Pt();
              hzzEventTree->genz2eta = ((*GenFakeMu1)+(*GenFakeMu2)).Eta();
              hzzEventTree->genz2rapidity = ((*GenFakeMu1)+(*GenFakeMu2)).Rapidity();
              hzzEventTree->genzzmass = ((*ZGenLepton1)+(*ZGenLepton2)+(*GenFakeMu1)+(*GenFakeMu2)).M();
              hzzEventTree->genzzpt = ((*ZGenLepton1)+(*ZGenLepton2)+(*GenFakeMu1)+(*GenFakeMu2)).Pt();
              hzzEventTree->genzzeta = ((*ZGenLepton1)+(*ZGenLepton2)+(*GenFakeMu1)+(*GenFakeMu2)).Eta();
              hzzEventTree->genzzrapidity = ((*ZGenLepton1)+(*ZGenLepton2)+(*GenFakeMu1)+(*GenFakeMu2)).Rapidity();

              hzzEventTree->l1id = ZLepton1PdgId;
              hzzEventTree->l1pt = ZLepton1Pt;
              hzzEventTree->l1eta = ZLepton1->Eta();
              hzzEventTree->l1phi = ZLepton1->Phi();
              hzzEventTree->l2id = ZLepton2PdgId;
              hzzEventTree->l2pt = ZLepton2Pt;
              hzzEventTree->l2eta = ZLepton2->Eta();
              hzzEventTree->l2phi = ZLepton2->Phi();
              hzzEventTree->l3id = 13;
              hzzEventTree->l3pt = fakeMu1->Pt();
              hzzEventTree->l3eta = fakeMu1->Eta();
              hzzEventTree->l3phi = fakeMu1->Phi();
              hzzEventTree->l4id = -13;
              hzzEventTree->l4pt = fakeMu2->Pt();
              hzzEventTree->l4eta = fakeMu2->Eta();
              hzzEventTree->l4phi = fakeMu2->Phi();

              hzzEventTree->channel = (abs(ZLepton1PdgId) == 11) ? HZZEventTree::kTwoETwoMu : HZZEventTree::kFourMu;
              hzzEventTree->z1mass = ((*ZLepton1)+(*ZLepton2)).M();
              hzzEventTree->z1pt = ((*ZLepton1)+(*ZLepton2)).M();          
              hzzEventTree->z1eta = ((*ZLepton1)+(*ZLepton2)).Eta();   
              hzzEventTree->z1rapidity = ((*ZLepton1)+(*ZLepton2)).Rapidity();   
              hzzEventTree->z2mass = ((*fakeMu1)+(*fakeMu2)).M();
              hzzEventTree->z2pt = ((*fakeMu1)+(*fakeMu2)).Pt();
              hzzEventTree->z2eta = ((*fakeMu1)+(*fakeMu2)).Eta();
              hzzEventTree->z2rapidity = ((*fakeMu1)+(*fakeMu2)).Rapidity();
              hzzEventTree->zzmass = ((*ZLepton1)+(*ZLepton2)+(*fakeMu1)+(*fakeMu2)).M();
              hzzEventTree->zzpt = ((*ZLepton1)+(*ZLepton2)+(*fakeMu1)+(*fakeMu2)).Pt();
              hzzEventTree->zzeta = ((*ZLepton1)+(*ZLepton2)+(*fakeMu1)+(*fakeMu2)).Eta();
              hzzEventTree->zzrapidity = ((*ZLepton1)+(*ZLepton2)+(*fakeMu1)+(*fakeMu2)).Rapidity();

              hzzEventTree->tree_->Fill();
            }
          }
        
          delete fakeMu1;
          delete fakeMu2;
          delete GenFakeMu1;
          delete GenFakeMu2;

          NEvents++;
          if (passMassCut) NEventsPass++;

          //*********************************************************************************************************
          //Fill Fake Ele + Fake Ele Case
          //*********************************************************************************************************
          passMassCut = false;
          TLorentzVector *fakeEle1 = new TLorentzVector();
          TLorentzVector *fakeEle2 = new TLorentzVector();
          fakeEle1->SetPtEtaPhiM( elePt_j , eleEta_j, elePhi_j, 0);
          fakeEle2->SetPtEtaPhiM( elePt_k , eleEta_k, elePhi_k, 0);
          TLorentzVector *GenFakeEle1 = new TLorentzVector();
          TLorentzVector *GenFakeEle2 = new TLorentzVector();
          GenFakeEle1->SetPtEtaPhiM( ptGenJet[j], etaGenJet[j], phiGenJet[j], 0);
          GenFakeEle2->SetPtEtaPhiM( ptGenJet[k], etaGenJet[k], phiGenJet[k], 0);

          mass = (*ZLepton1 + *ZLepton2 + *fakeEle1 + *fakeEle2).M();
          mass4L->Fill(mass);

          if (mass > 100 && mass < 150) {
            passMassCut = true;
          }

          //Fill event ntuple
          if (passMassCut) {

            double leadingLepPt = ZLepton1Pt;
            double subleadingLepPt = ZLepton2Pt;
            if (ZLepton2Pt > ZLepton1Pt) {
              leadingLepPt = ZLepton2Pt;
              subleadingLepPt = ZLepton1Pt;
            }
            if (elePt_j > leadingLepPt) {
              subleadingLepPt = leadingLepPt;
              leadingLepPt = elePt_j;      
            } else if (elePt_j > subleadingLepPt) {
              subleadingLepPt = elePt_j;
            }
            if (elePt_k > leadingLepPt) {
              subleadingLepPt = leadingLepPt;
              leadingLepPt = elePt_k;      
            } else if (elePt_k > subleadingLepPt) {
              subleadingLepPt = elePt_k;
            }

            if (ZLepton1Pt > (abs(ZLepton1PdgId) == 11 ? 7 : 5)
                && ZLepton2Pt > (abs(ZLepton2PdgId) == 11 ? 7 : 5)
                && elePt_j >  7
                && elePt_k >  7
                && fabs(ZLepton1->Eta()) < (abs(ZLepton1PdgId) == 11 ? 2.5 : 2.4)
                && fabs(ZLepton2->Eta()) < (abs(ZLepton2PdgId) == 11 ? 2.5 : 2.4)
                && fabs(eleEta_j) < 2.5
                && fabs(eleEta_k) < 2.5
                && leadingLepPt > 20
                && subleadingLepPt > 10
              ) {

              //*******************************************
              // Fill HZZ Event
              //*******************************************
              hzzEventTree->weight = ZLepton1Eff*ZLepton2Eff*eleFakerate_j*eleFakerate_k;
              hzzEventTree->run = 0;
              hzzEventTree->lumi = 0;
              hzzEventTree->event = NEvents;
              hzzEventTree->rho = 0;
              hzzEventTree->nvtx = 0;
              hzzEventTree->met = 0;

              hzzEventTree->genl1id = ZLepton1PdgId;
              hzzEventTree->genl1pt = ZLepton1GenPt;
              hzzEventTree->genl1eta = ZLepton1->Eta();
              hzzEventTree->genl1phi = ZLepton1->Phi();
              hzzEventTree->genl2id = ZLepton2PdgId;
              hzzEventTree->genl2pt = ZLepton2GenPt;
              hzzEventTree->genl2eta = ZLepton2->Eta();
              hzzEventTree->genl2phi = ZLepton2->Phi();
              hzzEventTree->genl3id = 11;
              hzzEventTree->genl3pt = ptGenJet[j];
              hzzEventTree->genl3eta = etaGenJet[j];
              hzzEventTree->genl3phi = phiGenJet[j];
              hzzEventTree->genl4id = -11;
              hzzEventTree->genl4pt = ptGenJet[k];
              hzzEventTree->genl4eta = etaGenJet[k];
              hzzEventTree->genl4phi = phiGenJet[k];

              hzzEventTree->genchannel = (abs(ZLepton1PdgId) == 11) ? HZZEventTree::kFourE : HZZEventTree::kTwoMuTwoE;
              hzzEventTree->genz1mass = ((*ZGenLepton1)+(*ZGenLepton2)).M();
              hzzEventTree->genz1pt = ((*ZGenLepton1)+(*ZGenLepton2)).Pt();
              hzzEventTree->genz1eta = ((*ZGenLepton1)+(*ZGenLepton2)).Eta();
              hzzEventTree->genz1rapidity = ((*ZGenLepton1)+(*ZGenLepton2)).Rapidity();
              hzzEventTree->genz2mass = ((*GenFakeEle1)+(*GenFakeEle2)).M();
              hzzEventTree->genz2pt = ((*GenFakeEle1)+(*GenFakeEle2)).Pt();
              hzzEventTree->genz2eta = ((*GenFakeEle1)+(*GenFakeEle2)).Eta();
              hzzEventTree->genz2rapidity = ((*GenFakeEle1)+(*GenFakeEle2)).Rapidity();
              hzzEventTree->genzzmass = ((*ZGenLepton1)+(*ZGenLepton2)+(*GenFakeEle1)+(*GenFakeEle2)).M();
              hzzEventTree->genzzpt = ((*ZGenLepton1)+(*ZGenLepton2)+(*GenFakeEle1)+(*GenFakeEle2)).Pt();
              hzzEventTree->genzzeta = ((*ZGenLepton1)+(*ZGenLepton2)+(*GenFakeEle1)+(*GenFakeEle2)).Eta();
              hzzEventTree->genzzrapidity = ((*ZGenLepton1)+(*ZGenLepton2)+(*GenFakeEle1)+(*GenFakeEle2)).Rapidity();

              hzzEventTree->l1id = ZLepton1PdgId;
              hzzEventTree->l1pt = ZLepton1Pt;
              hzzEventTree->l1eta = ZLepton1->Eta();
              hzzEventTree->l1phi = ZLepton1->Phi();
              hzzEventTree->l2id = ZLepton2PdgId;
              hzzEventTree->l2pt = ZLepton2Pt;
              hzzEventTree->l2eta = ZLepton2->Eta();
              hzzEventTree->l2phi = ZLepton2->Phi();
              hzzEventTree->l3id = 11;
              hzzEventTree->l3pt = fakeEle1->Pt();
              hzzEventTree->l3eta = fakeEle1->Eta();
              hzzEventTree->l3phi = fakeEle1->Phi();
              hzzEventTree->l4id = -11;
              hzzEventTree->l4pt = fakeEle2->Pt();
              hzzEventTree->l4eta = fakeEle2->Eta();
              hzzEventTree->l4phi = fakeEle2->Phi();

              hzzEventTree->channel = (abs(ZLepton1PdgId) == 11) ? HZZEventTree::kFourE : HZZEventTree::kTwoMuTwoE;
              hzzEventTree->z1mass = ((*ZLepton1)+(*ZLepton2)).M();
              hzzEventTree->z1pt = ((*ZLepton1)+(*ZLepton2)).M();          
              hzzEventTree->z1eta = ((*ZLepton1)+(*ZLepton2)).Eta();   
              hzzEventTree->z1rapidity = ((*ZLepton1)+(*ZLepton2)).Rapidity();   
              hzzEventTree->z2mass = ((*fakeEle1)+(*fakeEle2)).M();
              hzzEventTree->z2pt = ((*fakeEle1)+(*fakeEle2)).Pt();
              hzzEventTree->z2eta = ((*fakeEle1)+(*fakeEle2)).Eta();
              hzzEventTree->z2rapidity = ((*fakeEle1)+(*fakeEle2)).Rapidity();
              hzzEventTree->zzmass = ((*ZLepton1)+(*ZLepton2)+(*fakeEle1)+(*fakeEle2)).M();
              hzzEventTree->zzpt = ((*ZLepton1)+(*ZLepton2)+(*fakeEle1)+(*fakeEle2)).Pt();
              hzzEventTree->zzeta = ((*ZLepton1)+(*ZLepton2)+(*fakeEle1)+(*fakeEle2)).Eta();
              hzzEventTree->zzrapidity = ((*ZLepton1)+(*ZLepton2)+(*fakeEle1)+(*fakeEle2)).Rapidity();

              hzzEventTree->tree_->Fill();

              if (hzzEventTree->weight > 1e-6) {
                //cout << eleFakerate_j << " , " << eleFakerate_k << " : " << ptGenJet[j] << " " << etaGenJet[j] << " " << phiGenJet[j] << " , " << ptGenJet[k] << " " << etaGenJet[k] << " " << phiGenJet[k] << " : " << hzzEventTree->genl1pt << " " << hzzEventTree->genl1eta << " " << hzzEventTree->genl1phi << " | " << hzzEventTree->genl2pt << " " << hzzEventTree->genl2eta << " " << hzzEventTree->genl2phi << "\n"; 
              }
            }
          }

          delete fakeEle1;
          delete fakeEle2;
          delete GenFakeEle1;
          delete GenFakeEle2;

          NEvents++;
          if (passMassCut) NEventsPass++;

        } //loop over genjet k
      } //loop over genjet j


      if (ZLepton1) delete ZLepton1;
      if (ZLepton2) delete ZLepton2;
      if (ZGenLepton1) delete ZGenLepton1;
      if (ZGenLepton2) delete ZGenLepton2;

    } // loop over repeats

  } //loop over events


  cout << "NEventsPass: " << NEventsPass << " / Total ( " << NEvents << " )\n";
  mass4L->Draw();
  
  outputFile->cd();
  outputFile->Write();
  outputFile->Close();
  delete outputFile;
  if (hzzEventTree) delete hzzEventTree;
  gBenchmark->Show("test");

}

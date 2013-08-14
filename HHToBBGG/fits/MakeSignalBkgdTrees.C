//
//
// root -l CMSAna/HHToBBGG/fits/MakeSignalBkgdTrees.C+'("/afs/cern.ch/work/d/daan/public/HHToBBGG/normalizedNtuples/HHToBBGGNtuple.All_Updated.root")'
// 
// /afs/cern.ch/work/v/vlambert/public/HHbbggBackground/WorkingNtuples/signal_bkgd.root
//

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TGraphAsymmErrors.h>  
#include <TLegend.h>
#include <TCanvas.h>
#include "CMSAna/HHToBBGG/interface/HHToBBGGEventTree.hh"
#include <vector>
#include <set>
#include <map>
#include <string>
#include <sstream>
//edit
#include "CMSAna/Utils/CommonTools.hh"

void MakeSignalBkgdTrees(const string filename = "/afs/cern.ch/work/d/daan/public/HHToBBGG/normalizedNtuples/HHToBBGGNtuple.All_Updated.root") {

  double intLumi = 3000; //in units of fb^-1
  //histograms for diphoton mass
  TH1F *sigPhoMass = new TH1F("sigPhoMass", ";M_{#gamma#gamma} [GeV/c^{2}];Number of Events", 30, 100, 150);
  sigPhoMass->SetFillColor(kRed); sigPhoMass->SetLineColor(kRed);
  TH1F *bkgPhoMass = new TH1F("bkgPhoMass", ";M_{#gamma#gamma} [GeV/c^{2}];Number of Events", 30, 100, 150);
  bkgPhoMass->SetFillColor(kBlue); bkgPhoMass->SetLineColor(kBlue);
  TH1F *resPhoMass = new TH1F("resPhoMass", ";M_{#gamma#gamma} [GeV/c^{2}];Number of Events", 30, 100, 150);
  resPhoMass->SetFillColor(kOrange); resPhoMass->SetLineColor(kOrange);
  TH1F *nonresPhoMass = new TH1F("nonresPhoMass", ";M_{#gamma#gamma} [GeV/c^{2}];Number of Events", 30, 100, 150);
  nonresPhoMass->SetFillColor(kGreen); nonresPhoMass->SetLineColor(kGreen);
  TH1F *allPhoMass = new TH1F("allPhoMass", ";M_{#gamma#gamma} [GeV/c^{2}];Number of Events", 30, 100, 150);
  allPhoMass->SetFillColor(kBlack); allPhoMass->SetLineColor(kBlack);
  
  //histograms for dibjet mass
  TH1F *sigBjetMass = new TH1F("sigBjetMass", ";M_{bb} [GeV/c^{2}];Number of Events", 30, 70, 200);
  sigBjetMass->SetFillColor(kRed); sigBjetMass->SetLineColor(kRed);
  TH1F *bkgBjetMass = new TH1F("bkgBjetMass", ";M_{bb} [GeV/c^{2}];Number of Events", 30, 70, 200);
  bkgBjetMass->SetFillColor(kBlue); bkgBjetMass->SetLineColor(kBlue);
  TH1F *resBjetMass = new TH1F("resBjetMass", ";M_{bb} [GeV/c^{2}];Number of Events", 30, 70, 200);
  resBjetMass->SetFillColor(kOrange); resBjetMass->SetLineColor(kOrange);
  TH1F *resBjetMassExt = new TH1F("resBjetMassExt", ";M_{bb} [GeV/c^{2}];Number of Events", 30, 40, 200);
  resBjetMassExt->SetFillColor(kOrange); resBjetMassExt->SetLineColor(kOrange);
  TH1F *nonresBjetMass = new TH1F("nonresBjetMass", ";M_{bb} [GeV/c^{2}];Number of Events", 30, 70, 200);
  nonresBjetMass->SetFillColor(kGreen); nonresBjetMass->SetLineColor(kGreen);
  TH1F *allBjetMass = new TH1F("allBjetMass", ";M_{bb} [GeV/c^{2}];Number of Events", 30, 70, 200);
  allBjetMass->SetFillColor(kBlack); allBjetMass->SetLineColor(kBlack);
  
  //histograms for dibjet mass with photon window
  TH1F *sigBjetMassPhoWin = new TH1F("sigBjetMass", ";M_{bb} [GeV/c^{2}];Number of Events", 30, 70, 200);
  sigBjetMassPhoWin->SetFillColor(kRed); sigBjetMassPhoWin->SetLineColor(kRed);
  TH1F *resBjetMassPhoWin = new TH1F("resBjetMass", ";M_{bb} [GeV/c^{2}];Number of Events", 30, 70, 200);
  resBjetMassPhoWin->SetFillColor(kOrange); resBjetMassPhoWin->SetLineColor(kOrange);
  TH1F *resBjetMassExtPhoWin = new TH1F("resBjetMassExtPhoWin", ";M_{bb} [GeV/c^{2}];Number of Events", 30, 40, 200);
  resBjetMassExtPhoWin->SetFillColor(kOrange); resBjetMassExtPhoWin->SetLineColor(kOrange);
  TH1F *nonresBjetMassPhoWin = new TH1F("nonresBjetMass", ";M_{bb} [GeV/c^{2}];Number of Events", 30, 70, 200);
  nonresBjetMassPhoWin->SetFillColor(kGreen); nonresBjetMassPhoWin->SetLineColor(kGreen);
  TH1F *allBjetMassPhoWin = new TH1F("allBjetMass", ";M_{bb} [GeV/c^{2}];Number of Events", 30, 70, 200);
  allBjetMassPhoWin->SetFillColor(kBlack); allBjetMassPhoWin->SetLineColor(kBlack);
  
  //2D histogram of Mgg and Mbb
  TH2F *sigMassPhoBjet = new TH2F ("sigMassPhoBjet",";M_{bb} [GeV/c^{2}]; M_{#gamma#gamma} [GeV/c^{2}]; Number of Events", 30, 70, 200, 30, 100, 150);
  sigMassPhoBjet->SetFillColor(kRed); sigMassPhoBjet->SetLineColor(kRed);
  TH2F *resMassPhoBjet = new TH2F ("resMassPhoBjet",";M_{bb} [GeV/c^{2}]; M_{#gamma#gamma} [GeV/c^{2}]; Number of Events", 30, 70, 200, 30, 100, 150);
  resMassPhoBjet->SetFillColor(kOrange); resMassPhoBjet->SetLineColor(kOrange);
  TH2F *nonresMassPhoBjet = new TH2F ("nonresMassPhoBjet",";M_{bb} [GeV/c^{2}]; M_{#gamma#gamma} [GeV/c^{2}]; Number of Events", 30, 70, 200, 30, 100, 150);
  nonresMassPhoBjet->SetFillColor(kGreen); nonresMassPhoBjet->SetLineColor(kGreen);
  
  int eventsPassPhoWindow = 0;
  int eventsNormally = 0;
  //Fill Trees
  cmsana::HHToBBGGEventTree event;
  event.LoadTree(filename.c_str());
  event.InitTree();
  double totalweight = 0;
  for (int n=0;n<event.tree_->GetEntries();n++) { 
    event.tree_->GetEntry(n);
    double weight = event.weight*intLumi;
    totalweight += weight;
    if (event.bjet1.Pt() > 25.0 && event.bjet2.Pt() > 25.0
          && event.pho1.Pt() > 25.0 && event.pho2.Pt() > 25.0
          && fmax(event.pho1.Pt(),event.pho2.Pt()) > 40.0
          && event.nlep <= 0 && event.ncentraljets <= 3
          && event.DRgg < 2.0 && event.DRbb < 2.0 && event.minDRgb > 1.5
          && event.dibjet.M() > 40 && event.dibjet.M() < 200
          && event.diphoton.M() > 100 && event.diphoton.M() < 150)
    {
      eventsNormally++;
      if (event.sampletype >= 2 && event.sampletype <= 4 || event.sampletype == 20) {
        resBjetMassExt->Fill( event.dibjet.M() , weight);
        if (event.diphoton.M() > 120 && event.diphoton.M() < 130)
          resBjetMassExtPhoWin->Fill( event.dibjet.M() , weight);
      }
      
      if (event.dibjet.M() > 70) {
        //fill the full mass ranges of diphotons and dibjets
        allPhoMass->Fill( event.diphoton.M() , weight);
        allBjetMass->Fill( event.dibjet.M() , weight);
        if (event.sampletype == 1) {
          sigPhoMass->Fill( event.diphoton.M() , weight);
          sigBjetMass->Fill( event.dibjet.M() , weight);
          sigMassPhoBjet->Fill( event.dibjet.M() , event.diphoton.M() , weight);
        }
        if (event.sampletype >= 2 && event.sampletype <= 10 || event.sampletype == 20) {
          bkgPhoMass->Fill( event.diphoton.M() , weight);
          bkgBjetMass->Fill( event.dibjet.M() , weight);
        }
        if (event.sampletype >= 2 && event.sampletype <= 4 || event.sampletype == 20) {
          resPhoMass->Fill( event.diphoton.M() , weight);
          resBjetMass->Fill( event.dibjet.M() , weight);
          resMassPhoBjet->Fill( event.dibjet.M() , event.diphoton.M() , weight);
        }
        if (event.sampletype >= 5 && event.sampletype <= 10) {
          nonresPhoMass->Fill( event.diphoton.M() , weight);
          nonresBjetMass->Fill( event.dibjet.M() , weight);
          nonresMassPhoBjet->Fill( event.dibjet.M() , event.diphoton.M() , weight);
        }
        
        //fill dibjet mass after applying diphoton mass window
        if (event.diphoton.M() > 120 && event.diphoton.M() < 130) {
          eventsPassPhoWindow++;
          allBjetMassPhoWin->Fill( event.dibjet.M() , weight);
            if (event.sampletype == 1) {
           sigBjetMassPhoWin->Fill( event.dibjet.M() , weight);
          }
            if (event.sampletype >= 2 && event.sampletype <= 4 || event.sampletype == 20) {
            resBjetMassPhoWin->Fill( event.dibjet.M() , weight);
          }
          if (event.sampletype >= 5 && event.sampletype <= 10) {
            nonresBjetMassPhoWin->Fill( event.dibjet.M() , weight);
          }
        }
      }
    }
  }
  
  cout << eventsNormally << " | " << eventsPassPhoWindow << endl;
  cout << "diPhoton {sig, bkg, res, nonres, all}: { " << sigPhoMass->Integral() << ", " <<  bkgPhoMass->Integral() << ", " <<  resPhoMass->Integral() << ", " <<  nonresPhoMass->Integral() << ", " <<  allPhoMass->Integral() << " }" << endl;
  cout << "diBjet {sig, bkg, res, nonres, all}: { " << sigBjetMass->Integral() << ", " <<  bkgBjetMass->Integral() << ", " << resBjetMass->Integral() << ", " <<  nonresBjetMass->Integral() << ", " <<  allBjetMass->Integral() << " }" << endl;
  cout << "diBjet w/ PhoWin {sig, res, nonres, all}: { " << sigBjetMassPhoWin->Integral() << ", " <<  resBjetMassPhoWin->Integral() << ", " <<  nonresBjetMassPhoWin->Integral() << ", " <<  allBjetMassPhoWin->Integral() << " }" << endl;
  
  
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open("CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_diphotonMass.root", "UPDATE");
  file->cd();
  file->WriteTObject(sigPhoMass, sigPhoMass->GetName(), "WriteDelete");
  file->WriteTObject(bkgPhoMass, bkgPhoMass->GetName(), "WriteDelete");
  file->WriteTObject(resPhoMass, resPhoMass->GetName(), "WriteDelete");
  file->WriteTObject(nonresPhoMass, nonresPhoMass->GetName(), "WriteDelete");
  file->WriteTObject(allPhoMass, allPhoMass->GetName(), "WriteDelete");
  file->Close();
  delete file;
  
  TFile *file1 = TFile::Open("CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_dibjetMass.root", "UPDATE");
  file1->cd();
  file1->WriteTObject(sigBjetMass, sigBjetMass->GetName(), "WriteDelete");
  file1->WriteTObject(bkgBjetMass, bkgBjetMass->GetName(), "WriteDelete");
  file1->WriteTObject(resBjetMass, resBjetMass->GetName(), "WriteDelete");
  file1->WriteTObject(resBjetMassExt, resBjetMassExt->GetName(), "WriteDelete");
  file1->WriteTObject(nonresBjetMass, nonresBjetMass->GetName(), "WriteDelete");
  file1->WriteTObject(allBjetMass, allBjetMass->GetName(), "WriteDelete");
  file1->Close();
  delete file1;
  
  TFile *file2 = TFile::Open("CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_dibjetMassPhoWin.root", "UPDATE");
  file2->cd();
  file2->WriteTObject(sigBjetMassPhoWin, sigBjetMassPhoWin->GetName(), "WriteDelete");
  file2->WriteTObject(resBjetMassPhoWin, resBjetMassPhoWin->GetName(), "WriteDelete");
  file2->WriteTObject(resBjetMassExtPhoWin, resBjetMassExtPhoWin->GetName(), "WriteDelete");
  file2->WriteTObject(nonresBjetMassPhoWin, nonresBjetMassPhoWin->GetName(), "WriteDelete");
  file2->WriteTObject(allBjetMassPhoWin, allBjetMassPhoWin->GetName(), "WriteDelete");
  file2->Close();
  delete file2;
  
  TFile *file3 = TFile::Open("CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_twoDMass.root", "UPDATE");
  file3->cd();
  file3->WriteTObject(sigMassPhoBjet, sigMassPhoBjet->GetName(), "WriteDelete");
  file3->WriteTObject(resMassPhoBjet, resMassPhoBjet->GetName(), "WriteDelete");
  file3->WriteTObject(nonresMassPhoBjet, nonresMassPhoBjet->GetName(), "WriteDelete");
  file3->Close();
  delete file3;
  
}

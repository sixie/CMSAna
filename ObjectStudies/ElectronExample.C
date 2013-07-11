#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include "CMSAna/ObjectStudies/interface/ElectronTree.h"
#include <vector>
#include <map>


//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
void NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return;
}



//------------------------------------------------------------------------------
// PlotHiggsRes_LP
//------------------------------------------------------------------------------
void ElectronExample ( string InputFilename    = "ElectronNtuple.DYToEEAge0START.1.test.root")
{

  const double intLumi = 3000;

  //*******************************************************************************************
  //Define Histograms
  //*******************************************************************************************
  TH1F *hist_electronSigmaIEtaIEta =  new TH1F( "hist_electronSigmaIEtaIEta", "; #sigma_{i#eta i#eta};Number of Events", 50, 0, 0.05);


  //*******************************************************************************************
  //Do analysis
  //*******************************************************************************************
                
  cmsana::ElectronTree eletree;
  eletree.LoadTree(InputFilename.c_str());
  eletree.InitTree();

  for (int n=0;n<eletree.tree_->GetEntries();n++) { 
    eletree.tree_->GetEntry(n);

    if (eletree.fElePt > 20) {
      hist_electronSigmaIEtaIEta->Fill(eletree.fEleSigmaIEtaIEta);      
    }
  }
  
  
  //Draw Plots
  TCanvas *cv = 0;
  TLegend *legend = 0;

  //*******************************************************************************************
  //Normalize histograms
  //*******************************************************************************************
  NormalizeHist(hist_electronSigmaIEtaIEta);


  //*******************************************************************************************
  //Photon Eta
  //*******************************************************************************************

  cv = new TCanvas("cv","cv", 800,600);
  hist_electronSigmaIEtaIEta->Draw("hist");
  cv->SaveAs("electronSigmaIEtaIEta.gif");


  
}

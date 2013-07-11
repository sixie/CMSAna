#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include "CMSAna/ObjectStudies/interface/PhotonTree.h"
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
void PhotonExample ( string InputFilename    = "PhotonNtuple.root")
{

  const double intLumi = 3000;

  //*******************************************************************************************
  //Define Histograms
  //*******************************************************************************************
  TH1F *hist_photonSigmaIEtaIEta =  new TH1F( "hist_photonSigmaIEtaIEta", "; #sigma_{i#eta i#eta};Number of Events", 100, 0, 0.05);


  //*******************************************************************************************
  //Do analysis
  //*******************************************************************************************
                
  cmsana::PhotonTree photree;
  photree.LoadTree(InputFilename.c_str());
  photree.InitTree();

  for (int n=0;n<photree.tree_->GetEntries();n++) { 
    photree.tree_->GetEntry(n);

    if (photree.fPhoPt > 20 && fabs(photree.fPhoEta) < 1.5 ) {
      hist_photonSigmaIEtaIEta->Fill(photree.fPhoSigmaIEtaIEta);
    }
  }
  
  
  //Draw Plots
  TCanvas *cv = 0;
  TLegend *legend = 0;

  //*******************************************************************************************
  //Normalize histograms
  //*******************************************************************************************
  NormalizeHist(hist_photonSigmaIEtaIEta);


  //*******************************************************************************************
  //Photon Eta
  //*******************************************************************************************

  cv = new TCanvas("cv","cv", 800,600);
  hist_photonSigmaIEtaIEta->Draw("hist");
  cv->SaveAs("photonSigmaIEtaIEta.gif");


  
}

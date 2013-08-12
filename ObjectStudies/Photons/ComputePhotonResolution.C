//root -l CMSAna/ObjectStudies/Photons/ComputePhotonResolution.C+'("/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonStudies/ntuples/PhotonNtuple.real.combinedSummer12.root")'

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
void ComputePhotonResolution ( string filename = "PhotonNtuple1.root"  )
{

//   string Label = "";
//   if (label != "") Label = "_" + label;

  //*******************************************************************************************
  //Define Histograms
  //*******************************************************************************************
 
       //*** variables _1****
  vector<vector<TH1F* > > photonResponseFunctions;
  vector<double> PtBins;
  PtBins.push_back(25);
  PtBins.push_back(35);
  PtBins.push_back(45);
  PtBins.push_back(55);
  PtBins.push_back(65);
  PtBins.push_back(75);
  vector<double> EtaBins;
  EtaBins.push_back(0.0);
  EtaBins.push_back(0.25);
  EtaBins.push_back(0.50);
  EtaBins.push_back(0.75);
  EtaBins.push_back(1.00);
  EtaBins.push_back(1.25);
  EtaBins.push_back(1.50);
  EtaBins.push_back(1.75);
  EtaBins.push_back(2.00);
  EtaBins.push_back(2.25);

  for (int ipt = 0; ipt < PtBins.size(); ++ipt) {
    vector<TH1F*> tmpPhotonResponseFunctions;
    for (int ieta = 0; ieta < EtaBins.size(); ++ieta) {      
      tmpPhotonResponseFunctions.push_back( new TH1F( Form("PhotonResponseFunction_PtBin%d_EtaBin%d",ipt,ieta), "; (PhotonPt - GenPt)/GenPt; Number of Events", 100, -0.5,0.5 ) );
    }
    photonResponseFunctions.push_back(tmpPhotonResponseFunctions);
  }


  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  cmsana::PhotonTree photree;
  photree.LoadTree(filename.c_str());
  photree.InitTree();

  for (int n=0;n<photree.tree_->GetEntries();n++) { 
    photree.tree_->GetEntry(n);

    Int_t PtBin = -1;
    Int_t EtaBin = -1;

    if (!(photree.fPhoPt > 25 && fabs(photree.fPhoEta) < 2.5)) continue;

    if (photree.fPhoPt > PtBins[0] && photree.fPhoPt <= PtBins[1]) PtBin = 0;
    if (photree.fPhoPt > PtBins[1] && photree.fPhoPt <= PtBins[2]) PtBin = 1;
    if (photree.fPhoPt > PtBins[2] && photree.fPhoPt <= PtBins[3]) PtBin = 2;
    if (photree.fPhoPt > PtBins[3] && photree.fPhoPt <= PtBins[4]) PtBin = 3;
    if (photree.fPhoPt > PtBins[4] && photree.fPhoPt <= PtBins[5]) PtBin = 4;
    if (photree.fPhoPt > PtBins[5] ) PtBin = 5;

    if (fabs(photree.fPhoEta) > EtaBins[0] && fabs(photree.fPhoEta) <= EtaBins[1]) EtaBin = 0;
    if (fabs(photree.fPhoEta) > EtaBins[1] && fabs(photree.fPhoEta) <= EtaBins[2]) EtaBin = 1;
    if (fabs(photree.fPhoEta) > EtaBins[2] && fabs(photree.fPhoEta) <= EtaBins[3]) EtaBin = 2;
    if (fabs(photree.fPhoEta) > EtaBins[3] && fabs(photree.fPhoEta) <= EtaBins[4]) EtaBin = 3;
    if (fabs(photree.fPhoEta) > EtaBins[4] && fabs(photree.fPhoEta) <= EtaBins[5]) EtaBin = 4;
    if (fabs(photree.fPhoEta) > EtaBins[5] && fabs(photree.fPhoEta) <= EtaBins[6]) EtaBin = 5;
    if (fabs(photree.fPhoEta) > EtaBins[6] && fabs(photree.fPhoEta) <= EtaBins[7]) EtaBin = 6;
    if (fabs(photree.fPhoEta) > EtaBins[7] && fabs(photree.fPhoEta) <= EtaBins[8]) EtaBin = 7;
    if (fabs(photree.fPhoEta) > EtaBins[8] && fabs(photree.fPhoEta) <= EtaBins[9]) EtaBin = 8;
    if (fabs(photree.fPhoEta) > EtaBins[9]) EtaBin = 9;

    assert(PtBin >= 0 && PtBin < PtBins.size());
    assert(EtaBin >= 0 && EtaBin < EtaBins.size());

    photonResponseFunctions[PtBin][EtaBin]->Fill( (photree.fPhoPt - photree.fPhoGenPt) / photree.fPhoGenPt);
  }
  

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open("PhotonResolution.root", "UPDATE");
  for (int ipt = 0; ipt < PtBins.size(); ++ipt) {
    for (int ieta = 0; ieta < EtaBins.size(); ++ieta) {      
      file->WriteTObject(photonResponseFunctions[ipt][ieta], photonResponseFunctions[ipt][ieta]->GetName(), "WriteDelete");
    }
  }
 
  file->Close();
  delete file;       

}

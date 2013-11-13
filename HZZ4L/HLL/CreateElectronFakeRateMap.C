//root -l CMSAna/HZZ4L/HLL/CreateEfficiencyMap.C+'("HZZEfficiencyMap_HZZ125.root","HZZ125",0)'


//================================================================================================
//
// Create Efficiency Map
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2F.h>                   // 1D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <TGraphAsymmErrors.h>      // class for parsing strings

// data structs
#include "CMSAna/HZZ4L/interface/HZZEfficiencyMap.hh"
#include "CMSAna/Utils/EfficiencyUtils.hh"


#include "CMSAna/ObjectStudies/interface/ElectronTree.h"
#include "CMSAna/ObjectStudies/interface/MuonTree.h"

// RooFit headers
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooGaussian.h"

#endif
 

//=== FUNCTION DECLARATIONS ======================================================================================


void initialize2DArray(vector<vector<double> > &array, UInt_t NPtBins, UInt_t NEtaBins) {

  array.resize(NPtBins+2);
  for (uint i=0; i < NPtBins+2; ++i) {
    array[i].resize(NEtaBins+2);
    for (uint j=0; j < NEtaBins+2; ++j) {
      array[i][j]= 0;
    }
  }
}


void initializeHistograms(vector<vector<TH1F*> > &hists, string name, UInt_t NPtBins, UInt_t NEtaBins, UInt_t nbins, Double_t xmin, Double_t xmax) {

  hists.resize(NPtBins+2);
  for (uint i=0; i < NPtBins+2; ++i) {
    hists[i].resize(NEtaBins+2);
    for (uint j=0; j < NEtaBins+2; ++j) {      
      hists[i][j]= new TH1F(Form("%s_PtBin%d_EtaBin%d",name.c_str(),i,j), "; (Reco Pt - Gen Pt)/Gen Pt; Number of Events", nbins, xmin, xmax);
    }
  }
}
  

UInt_t FindBin( double value, double bins[], UInt_t nbins) {

  UInt_t nbinboundaries = nbins+1;
  UInt_t bin = 0;
  for (uint i=0; i < nbinboundaries; ++i) {
    if (i < nbinboundaries-1) {
      if (value >= bins[i] && value < bins[i+1]) {
        bin = i+1;
        break;
      }
    } else if (i == nbinboundaries-1) {
      if (value >= bins[i]) {
        bin = nbinboundaries;
        break;
      }
    }    
  }
  return bin;
}
  

void computeEfficiencyPtEta(vector<vector<double> > &numerator, 
                            vector<vector<double> > &denominator,
                            vector<vector<double> > &eff,
                            int leptonType
  ) {

  for (uint i=0; i < numerator.size(); ++i) {
    for (uint j=0; j < numerator[i].size(); ++j) {
      
      if ( denominator[i][j] > 0) {
        eff[i][j] = (numerator[i][j] / denominator[i][j]);
      } else {
        eff[i][j] = 0;
      }
    }
  }
}


//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, const char* tname)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = TFile::Open(infname,"read");
  assert(inf);

  TTree* t = (TTree*)inf->Get(tname);
  
  if (!t) {
    TDirectory *dir = (TDirectory*)inf->FindObjectAny("HwwNtuplerMod");
    if (!dir) {
      cout << "Cannot get Directory HwwNtuplerMod from file " << infname << endl;
      assert(dir);
    }
    t = (TTree*)dir->Get(tname);
  }

  if (!t) {
    cout << "Cannot get Tree with name " << tname << " from file " << infname << endl;
  }
  assert(t);


  if (verbose) {
    cout << "---\tRecovered tree " << t->GetName()
	 << " with "<< t->GetEntries() << " entries" << endl;
  }
  
  return t;
}




//=== MAIN MACRO =================================================================================================

void doCreateElectronFakeRateMap(const string filename, const string Label = "ZZ", Int_t Option = 11) 
{  
  gBenchmark->Start("HZZTemplate");
  string label = Label;
  if (Label != "") label = "_" + Label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================

  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  TH1F *histDenominatorPt = new TH1F ("histDenominatorPt",";Electron p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorPt = new TH1F ("histNumeratorPt",";Electron p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorEta = new TH1F ("histDenominatorEta",";Electron Eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histNumeratorEta = new TH1F ("histNumeratorEta",";Electron Eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histDenominatorPhi = new TH1F ("histDenominatorPhi",";Electron Phi; Number of Events", 50, 0 , 3.2);
  TH1F *histNumeratorPhi = new TH1F ("histNumeratorPhi",";Electron Phi; Number of Events", 50, 0 , 3.2);
  TH1F *histDenominatorRho = new TH1F ("histDenominatorRho",";Electron Rho; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorRho = new TH1F ("histNumeratorRho",";Electron Rho; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpv = new TH1F ("histDenominatorNpv",";Electron Npv; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpv = new TH1F ("histNumeratorNpv",";Electron Npv; Number of Events", 50, 0 , 100);

  TH2F *histDenominatorPtEta = new TH2F ("histDenominatorPtEta",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", 50, 0 , 200, 50, 0, 3.0);
  TH2F *histNumeratorPtEta = new TH2F ("histNumeratorPtEta",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", 50, 0 , 200, 50, 0, 3.0);

  TH1F *histDenominatorPt_lq = new TH1F ("histDenominatorPt_lq",";Electron p_{T} [GeV/c^{2}]; Number of Events", 25, 0 , 100);
  TH1F *histNumeratorPt_lq = new TH1F ("histNumeratorPt_lq",";Electron p_{T} [GeV/c^{2}]; Number of Events", 25, 0 , 100);
  TH1F *histDenominatorPt_b = new TH1F ("histDenominatorPt_b",";Electron p_{T} [GeV/c^{2}]; Number of Events", 25, 0 , 100);
  TH1F *histNumeratorPt_b = new TH1F ("histNumeratorPt_b",";Electron p_{T} [GeV/c^{2}]; Number of Events", 25, 0 , 100);
  TH1F *histDenominatorPt_g = new TH1F ("histDenominatorPt_g",";Electron p_{T} [GeV/c^{2}]; Number of Events", 25, 0 , 100);
  TH1F *histNumeratorPt_g = new TH1F ("histNumeratorPt_g",";Electron p_{T} [GeV/c^{2}]; Number of Events", 25, 0 , 100);
  TH2F *histDenominatorPtEta_lq = new TH2F ("histDenominatorPtEta_lq",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", 50, 0 , 200, 50, 0, 3.0);
  TH2F *histNumeratorPtEta_lq = new TH2F ("histNumeratorPtEta_lq",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", 50, 0 , 200, 50, 0, 3.0);
  TH2F *histDenominatorPtEta_b = new TH2F ("histDenominatorPtEta_b",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", 50, 0 , 200, 50, 0, 3.0);
  TH2F *histNumeratorPtEta_b = new TH2F ("histNumeratorPtEta_b",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", 50, 0 , 200, 50, 0, 3.0);
  TH2F *histDenominatorPtEta_g = new TH2F ("histDenominatorPtEta_g",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", 50, 0 , 200, 50, 0, 3.0);
  TH2F *histNumeratorPtEta_g = new TH2F ("histNumeratorPtEta_g",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", 50, 0 , 200, 50, 0, 3.0);


  //********************************************************
  // Create Arrays to store the map
  //********************************************************
   const UInt_t NPtBins_Response = 4; 
   const UInt_t NEtaBins_Response = 2;
   double ptBins_Response[NPtBins_Response+1] = { 5, 7, 10, 15, 25 };
   double etaBins_Response[NEtaBins_Response+1] = { 0.0, 1.5, 2.5 };

   const UInt_t NPtBins = 7; 
   const UInt_t NEtaBins = 3;
   double ptBins[NPtBins+1] = { 5, 7, 10, 15, 20, 25, 30, 35 };
   double etaBins[NEtaBins+1] = { 0.0, 1.0, 2.0, 2.5 };

   vector<vector<double> > NDenominator_Electrons_PtEta;
   vector<vector<double> > NNumerator_Electrons_PtEta;
   vector<vector<double> > Efficiency_Electrons_PtEta;
   vector<vector<TH1F*> > PtResolution_PtEta_Electrons;

   initialize2DArray(NDenominator_Electrons_PtEta, NPtBins, NEtaBins);
   initialize2DArray(NNumerator_Electrons_PtEta, NPtBins, NEtaBins);
   initialize2DArray(Efficiency_Electrons_PtEta, NPtBins, NEtaBins);
   initializeHistograms(PtResolution_PtEta_Electrons, "LeptonPtResolution_Electrons", NPtBins_Response, NEtaBins_Response, 100, -1.0, 0.25);  

   
   //--------------------------------------------------------------------------------------------------------------
   // Read efficiency map ntuple
   //==============================================================================================================  
   cmsana::ElectronTree eleTree;
   eleTree.LoadTree(filename.c_str());
   eleTree.InitTree(cmsana::ElectronTree::kEleTreeLight);

   cout << "Total : " << eleTree.tree_->GetEntries() << "\n";
  for(UInt_t i=0; i<eleTree.tree_->GetEntries(); i++) {       	
    eleTree.tree_->GetEntry(i);
    if (i % 1000000 == 0) cout << "Ele " << i << endl;

    Int_t tmpPtBin = FindBin( eleTree.fEleGenPt , ptBins, NPtBins);
    Int_t tmpEtaBin = FindBin( fabs(eleTree.fEleGenEta) , etaBins, NEtaBins);
    Int_t tmpPtBin_Response = FindBin( eleTree.fEleGenPt , ptBins_Response, NPtBins_Response);
    Int_t tmpEtaBin_Response = FindBin( fabs(eleTree.fEleGenEta) , etaBins_Response, NEtaBins_Response);

    //selection cuts
    if (!(eleTree.fEleGenPt > 5 && fabs(eleTree.fEleGenEta) < 2.5)) continue;

    //**** PT - ETA ****
    NDenominator_Electrons_PtEta[tmpPtBin][tmpEtaBin] += 1.0;
    histDenominatorPtEta->Fill(eleTree.fEleGenPt,fabs(eleTree.fEleGenEta));
    histDenominatorPt->Fill(eleTree.fEleGenPt);
    histDenominatorEta->Fill(eleTree.fEleGenEta);
    histDenominatorRho->Fill(eleTree.fRho);
    histDenominatorNpv->Fill(eleTree.fNVertices);

    if (abs(eleTree.fPdgId) == 5) {
      histDenominatorPtEta_b->Fill(eleTree.fEleGenPt,fabs(eleTree.fEleGenEta));
      histDenominatorPt_b->Fill(eleTree.fEleGenPt);
    } else if (abs(eleTree.fPdgId) == 21) {
      histDenominatorPtEta_g->Fill(eleTree.fEleGenPt,fabs(eleTree.fEleGenEta));
      histDenominatorPt_g->Fill(eleTree.fEleGenPt);
    } else {
      histDenominatorPtEta_lq->Fill(eleTree.fEleGenPt,fabs(eleTree.fEleGenEta));
      histDenominatorPt_lq->Fill(eleTree.fEleGenPt);
    }

    if(eleTree.fElePt > 0) {
      NNumerator_Electrons_PtEta[tmpPtBin][tmpEtaBin] += 1.0;      
      histNumeratorPtEta->Fill(eleTree.fEleGenPt,fabs(eleTree.fEleGenEta));
      histNumeratorPt->Fill(eleTree.fEleGenPt);
      histNumeratorEta->Fill(eleTree.fEleGenEta);
      histNumeratorRho->Fill(eleTree.fRho);
      histNumeratorNpv->Fill(eleTree.fNVertices);

      if (abs(eleTree.fPdgId) == 5) {
        histNumeratorPtEta_b->Fill(eleTree.fEleGenPt,fabs(eleTree.fEleGenEta));
        histNumeratorPt_b->Fill(eleTree.fEleGenPt);
      } else if (abs(eleTree.fPdgId) == 21) {
        histNumeratorPtEta_g->Fill(eleTree.fEleGenPt,fabs(eleTree.fEleGenEta));        
        histNumeratorPt_g->Fill(eleTree.fEleGenPt);
      } else {
        histNumeratorPtEta_lq->Fill(eleTree.fEleGenPt,fabs(eleTree.fEleGenEta));
        histNumeratorPt_lq->Fill(eleTree.fEleGenPt);
      }

      //fill response function
      PtResolution_PtEta_Electrons[tmpPtBin_Response][tmpEtaBin_Response]->Fill( (eleTree.fElePt - eleTree.fEleGenPt)/eleTree.fEleGenPt , 1.0);
      
    }
    
  }

  //--------------------------------------------------------------------------------------------------------------
  // Make Efficiency Plots
  //==============================================================================================================

  TGraphAsymmErrors *efficiency_pt = cmsana::createEfficiencyGraph(histNumeratorPt, histDenominatorPt, "Efficiency_Pt" , vector<double>() ,  -99, -99, 0, 1, false);
  TGraphAsymmErrors *efficiency_eta = cmsana::createEfficiencyGraph(histNumeratorEta, histDenominatorEta, "Efficiency_Eta" , vector<double>() ,  -99, -99, 0, 1, false);
  TGraphAsymmErrors *efficiency_rho = cmsana::createEfficiencyGraph(histNumeratorRho, histDenominatorRho, "Efficiency_Rho" , vector<double>() ,  -99, -99, 0, 1, false);
  TGraphAsymmErrors *efficiency_npv = cmsana::createEfficiencyGraph(histNumeratorNpv, histDenominatorNpv, "Efficiency_Npv" , vector<double>() ,  -99, -99, 0, 1, false);
  TH2F *efficiency_pteta = cmsana::createEfficiencyHist2D(histNumeratorPtEta, histDenominatorPtEta, "Efficiency_PtEta" , vector<double>() ,vector<double>());  


  TGraphAsymmErrors *efficiency_lq_pt = cmsana::createEfficiencyGraph(histNumeratorPt_lq, histDenominatorPt_lq, "Efficiency_lq_Pt" , vector<double>() ,  -99, -99, 0, 1, false);
  TGraphAsymmErrors *efficiency_b_pt = cmsana::createEfficiencyGraph(histNumeratorPt_b, histDenominatorPt_b, "Efficiency_b_Pt" , vector<double>() ,  -99, -99, 0, 1, false);
  TGraphAsymmErrors *efficiency_g_pt = cmsana::createEfficiencyGraph(histNumeratorPt_g, histDenominatorPt_g, "Efficiency_g_Pt" , vector<double>() ,  -99, -99, 0, 1, false);
  TH2F *efficiency_lq_pteta = cmsana::createEfficiencyHist2D(histNumeratorPtEta_lq, histDenominatorPtEta_lq, "Efficiency_lq_PtEta" , vector<double>() ,vector<double>());  
  TH2F *efficiency_b_pteta = cmsana::createEfficiencyHist2D(histNumeratorPtEta_b, histDenominatorPtEta_b, "Efficiency_b_PtEta" , vector<double>() ,vector<double>());  
  TH2F *efficiency_g_pteta = cmsana::createEfficiencyHist2D(histNumeratorPtEta_g, histDenominatorPtEta_g, "Efficiency_g_PtEta" , vector<double>() ,vector<double>());  

  //--------------------------------------------------------------------------------------------------------------
  // Draw
  //==============================================================================================================
  TCanvas *cv =0;

  cv = new TCanvas("cv","cv",800,600);
  efficiency_pt->Draw("AP");
  //efficiency_pt->SetTitle("");
  efficiency_pt->GetYaxis()->SetRangeUser(0.0,0.1);
  cv->SaveAs("Efficiency_Pt.gif");

  cv = new TCanvas("cv","cv",800,600);
  efficiency_eta->Draw("AP");
  //efficiency_eta->SetTitle("");
  efficiency_eta->GetYaxis()->SetRangeUser(0.0,0.1);
  cv->SaveAs("Efficiency_Eta.gif");

  cv = new TCanvas("cv","cv",800,600);
  efficiency_rho->Draw("AP");
  //efficiency_rho->SetTitle("");
  efficiency_rho->GetYaxis()->SetRangeUser(0.0,0.1);
  cv->SaveAs("Efficiency_Rho.gif");

  cv = new TCanvas("cv","cv",800,600);
  efficiency_npv->Draw("AP");
  //efficiency_npv->SetTitle("");
  efficiency_npv->GetYaxis()->SetRangeUser(0.0,0.1);
  cv->SaveAs("Efficiency_Npv.gif");


  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("FakeRate"+label+".root").c_str(), "UPDATE");
  file->cd();
  file->WriteTObject(efficiency_pt, "Efficiency_Pt", "WriteDelete");
  file->WriteTObject(efficiency_eta, "Efficiency_Eta", "WriteDelete");
  file->WriteTObject(efficiency_rho, "Efficiency_Rho", "WriteDelete");
  file->WriteTObject(efficiency_npv, "Efficiency_NPV", "WriteDelete");
  file->WriteTObject(efficiency_pteta, "Efficiency_PtEta", "WriteDelete");

  file->WriteTObject(efficiency_lq_pt, "Efficiency_lq_Pt", "WriteDelete");
  file->WriteTObject(efficiency_b_pt, "Efficiency_b_Pt", "WriteDelete");
  file->WriteTObject(efficiency_g_pt, "Efficiency_g_Pt", "WriteDelete");
  file->WriteTObject(efficiency_b_pteta, "Efficiency_b_PtEta", "WriteDelete");
  file->WriteTObject(efficiency_lq_pteta, "Efficiency_lq_PtEta", "WriteDelete");
  file->WriteTObject(efficiency_g_pteta, "Efficiency_g_PtEta", "WriteDelete");

  for (uint i=0; i < NPtBins_Response+2; ++i) {
    for (uint j=0; j < NEtaBins_Response+2; ++j) {
      file->WriteTObject(PtResolution_PtEta_Electrons[i][j], PtResolution_PtEta_Electrons[i][j]->GetName(), "WriteDelete");
    }
  }

  file->Close();
  delete file;       


   computeEfficiencyPtEta(NNumerator_Electrons_PtEta, NDenominator_Electrons_PtEta, Efficiency_Electrons_PtEta,11);


  //********************************************************
  // Produce output lookup table
  //******************************************************** 
  ofstream outf_e("FakeElectronEfficiencyMap.h");

  outf_e << "UInt_t FindElectronEfficiencyBin( double value, double bins[], UInt_t nbins) {" << endl;
  outf_e << "  UInt_t nbinboundaries = nbins+1;" << endl;
  outf_e << "  UInt_t bin = 0;" << endl;
  outf_e << "  for (uint i=0; i < nbinboundaries; ++i) {" << endl;
  outf_e << "    if (i < nbinboundaries-1) {" << endl;
  outf_e << "      if (value >= bins[i] && value < bins[i+1]) {" << endl;
  outf_e << "        bin = i+1;" << endl;
  outf_e << "        break;" << endl;
  outf_e << "      }" << endl;
  outf_e << "    } else if (i == nbinboundaries-1) {" << endl;
  outf_e << "      if (value >= bins[i]) {" << endl;
  outf_e << "        bin = nbinboundaries;" << endl;
  outf_e << "        break;" << endl;
  outf_e << "      }" << endl;
  outf_e << "    }    " << endl;
  outf_e << "  }" << endl;
  outf_e << "  return bin;" << endl;
  outf_e << "}" << endl;

  outf_e << endl;
  outf_e << endl;

  outf_e << "Double_t GetElectronEfficiencyPtEta(Double_t Pt, Double_t Eta) {" << endl;

  outf_e << endl;
  outf_e << "  Double_t ptBins[" << NPtBins+1 << "] = {";
  for (uint i=0; i < NPtBins+1; ++i) {
    outf_e << ptBins[i];
    if (i < NPtBins) {
      outf_e << ",";
    }
  }
  outf_e << "};\n";

  outf_e << "  Double_t etaBins[" << NEtaBins+1 << "] = {";
  for (uint i=0; i < NEtaBins+1; ++i) {
    outf_e << etaBins[i];
    if (i < NEtaBins) {
      outf_e << ",";
    }
  }
  outf_e << "};\n";


  outf_e << endl;
  outf_e << endl;

  outf_e << "  Double_t Efficiency[" << NPtBins+2 << "][" << NEtaBins+2 << "] = {";
  outf_e << endl;

  for (uint i=0; i < NPtBins+2; ++i) {
    outf_e << "    {";
    for (uint j=0; j < NEtaBins+2; ++j) {
      outf_e << Efficiency_Electrons_PtEta[i][j];
      if (j< NEtaBins+1) {
        outf_e << ",";
      }
    }
    if (i< NPtBins+1) {
      outf_e << "    },";
    } else {
      outf_e << "}";
    }
    outf_e << endl;
  }
  
  outf_e << "  };" << endl;

  outf_e << endl;
  outf_e << endl;

  outf_e << "  Int_t tmpPtBin = FindElectronEfficiencyBin( Pt , ptBins, " << NPtBins << ");" << endl;
  outf_e << "  Int_t tmpEtaBin = FindElectronEfficiencyBin( Eta , etaBins, " << NEtaBins << ");" << endl;
  outf_e << "  return Efficiency[tmpPtBin][tmpEtaBin];" << endl;
  outf_e << "}" << endl;


  outf_e.close();


  gBenchmark->Show("WWTemplate");       
} 


void FitLeptonResponseModels(const string Label = "ZZ", Int_t Option = 0, Int_t PtBin = -1, Int_t EtaBin = -1) {

  string label = Label;
  if (Label != "") label = "_" + Label;
  TFile *fileInput = new TFile(("FakeRate" + label + ".root").c_str(), "READ");

  //********************************************************
  // Bins
  //********************************************************
  const UInt_t NPtBins = 4;
  const UInt_t NEtaBins = 2;
  double ptBins[NPtBins+1] = { 5, 7, 10, 15, 25 };
  double etaBins[NEtaBins+1] = { 0.0, 1.5, 2.5 };

  
  //********************************************************
  // Output
  //********************************************************
  TFile *fileOutput = new TFile(("FakeElectonPtResolutionModel" + label + ".root").c_str(), "UPDATE");
  
  TH2F* GaussParamArray_Electrons_mean = (TH2F*)fileOutput->Get("GaussParamArray_Electrons_mean");
  TH2F* GaussParamArray_Electrons_sigma = (TH2F*)fileOutput->Get("GaussParamArray_Electrons_sigma");


  if (!GaussParamArray_Electrons_mean) {
    GaussParamArray_Electrons_mean = new TH2F( "GaussParamArray_Electrons_mean", "", NPtBins, 0, NPtBins, NEtaBins, 0, NEtaBins);
    for (uint i=0; i < NPtBins+2; ++i) {
      for (uint j=0; j < NEtaBins+2; ++j) {
        GaussParamArray_Electrons_mean->SetBinContent(i,j,0.0);
      }
    }
  }
  if (!GaussParamArray_Electrons_sigma) {
    GaussParamArray_Electrons_sigma = new TH2F( "GaussParamArray_Electrons_sigma", "", NPtBins, 0, NPtBins, NEtaBins, 0, NEtaBins);
    for (uint i=0; i < NPtBins+2; ++i) {
      for (uint j=0; j < NEtaBins+2; ++j) {
        GaussParamArray_Electrons_sigma->SetBinContent(i,j,0.0);
      }
    }
  }


  //********************************************************
  // Fit for resolution function 
  //******************************************************** 
  for (uint i=0; i < NPtBins+2; ++i) {
    for (uint j=0; j < NEtaBins+2; ++j) {

      if (PtBin >= 0 && EtaBin >= 0) {
        if (!(i==PtBin && j==EtaBin)) continue;
      }
      
      TH1F* hist = (TH1F*)fileInput->Get(Form("LeptonPtResolution_Electrons_PtBin%d_EtaBin%d",i,j));
      assert(hist);

      RooRealVar leptonPtRes("leptonPtRes","leptonPtRes",-1.0,0.25);
      leptonPtRes.setRange("range",-1.00,0.25);
      leptonPtRes.setBins(100);
      RooDataHist *data=0;
      data = new RooDataHist("data","data",RooArgSet(leptonPtRes),hist);

      RooRealVar     *mean  = new RooRealVar("mean","mean",-0.25,-10,10);
      RooRealVar     *sigma = new RooRealVar("sigma","sigma",0.05,0.00001,0.5);

      RooGaussian *model = new RooGaussian("LeptonPtResModel","LeptonPtResModel",leptonPtRes,*mean,*sigma);

      RooFitResult *fitResult=0;
      fitResult = model->fitTo(*data,
                               RooFit::Binned(true),
                               RooFit::Strategy(1),
                               RooFit::Save());

      cout << "Fitted parameters\n";
      cout << mean->getVal() << endl;
      cout << sigma->getVal() << endl;
     
      if (i >= 1 && ( j >= 1 && j <= NEtaBins )) {
        GaussParamArray_Electrons_mean->SetBinContent(i,j,mean->getVal());
        GaussParamArray_Electrons_sigma->SetBinContent(i,j,sigma->getVal());
      } else {
        GaussParamArray_Electrons_mean->SetBinContent(i,j,0);
        GaussParamArray_Electrons_sigma->SetBinContent(i,j,0);        
      }



      //Save Workspace
      RooWorkspace *w = new RooWorkspace(Form("LeptonPtResolutionModel_Electrons_PtBin%d_EtaBin%d",i,j),Form("LeptonPtResolutionModel_Electrons_PtBin%d_EtaBin%d",i,j));
      w->import(*model);
      w->import(*data);
      //w->Print();

      //Make Plot
      RooPlot *frame = leptonPtRes.frame(RooFit::Bins(100));
      data->plotOn(frame,RooFit::MarkerStyle(kFullCircle),RooFit::MarkerSize(0.8),RooFit::DrawOption("ZP"));    
      model->plotOn(frame);

      TCanvas *cv = new TCanvas("cv","cv",800,600);
      
      frame->Draw();
      cv->SaveAs(Form("LeptonPtResolution_Electrons%s_PtBin%d_EtaBin%d.gif",label.c_str(),i,j)); 

      fileOutput->WriteTObject(model, Form("LeptonPtResolutionModel_Electrons_PtBin%d_EtaBin%d",i,j), "WriteDelete");
      fileOutput->WriteTObject(cv, Form("LeptonPtResolutionFit_Electrons_PtBin%d_EtaBin%d",i,j), "WriteDelete");
      fileOutput->WriteTObject(w, w->GetName(), "WriteDelete");

      fileOutput->WriteTObject(GaussParamArray_Electrons_mean, "GaussParamArray_Electrons_mean", "WriteDelete");
      fileOutput->WriteTObject(GaussParamArray_Electrons_sigma, "GaussParamArray_Electrons_sigma", "WriteDelete");


    }
  }

  //********************************************************
  // Produce output lookup table
  //******************************************************** 
  ofstream outf_e("FakeElectronResponseMap.h");

  outf_e << "UInt_t FindFakeElectronResponseBin( double value, double bins[], UInt_t nbins) {" << endl;
  outf_e << "  UInt_t nbinboundaries = nbins+1;" << endl;
  outf_e << "  UInt_t bin = 0;" << endl;
  outf_e << "  for (uint i=0; i < nbinboundaries; ++i) {" << endl;
  outf_e << "    if (i < nbinboundaries-1) {" << endl;
  outf_e << "      if (value >= bins[i] && value < bins[i+1]) {" << endl;
  outf_e << "        bin = i+1;" << endl;
  outf_e << "        break;" << endl;
  outf_e << "      }" << endl;
  outf_e << "    } else if (i == nbinboundaries-1) {" << endl;
  outf_e << "      if (value >= bins[i]) {" << endl;
  outf_e << "        bin = nbinboundaries;" << endl;
  outf_e << "        break;" << endl;
  outf_e << "      }" << endl;
  outf_e << "    }    " << endl;
  outf_e << "  }" << endl;
  outf_e << "  return bin;" << endl;
  outf_e << "}" << endl;

  outf_e << endl;
  outf_e << endl;

  outf_e << "Double_t GetElectronResponseMeanPtEta(Double_t Pt, Double_t Eta) {" << endl;

  outf_e << endl;
  outf_e << "  Double_t ptBins[" << NPtBins+1 << "] = {";
  for (uint i=0; i < NPtBins+1; ++i) {
    outf_e << ptBins[i];
    if (i < NPtBins) {
      outf_e << ",";
    }
  }
  outf_e << "};\n";

  outf_e << "  Double_t etaBins[" << NEtaBins+1 << "] = {";
  for (uint i=0; i < NEtaBins+1; ++i) {
    outf_e << etaBins[i];
    if (i < NEtaBins) {
      outf_e << ",";
    }
  }
  outf_e << "};\n";


  outf_e << endl;
  outf_e << endl;

  outf_e << "  Double_t ResponseMean[" << NPtBins+2 << "][" << NEtaBins+2 << "] = {";
  outf_e << endl;

  for (uint i=0; i < NPtBins+2; ++i) {
    outf_e << "    {";
    for (uint j=0; j < NEtaBins+2; ++j) {
      outf_e << GaussParamArray_Electrons_mean->GetBinContent(i,j);
      if (j< NEtaBins+1) {
        outf_e << ",";
      }
    }
    if (i< NPtBins+1) {
      outf_e << "    },";
    } else {
      outf_e << "}";
    }
    outf_e << endl;
  }
  
  outf_e << "  };" << endl;

  outf_e << endl;
  outf_e << endl;

  outf_e << "  Int_t tmpPtBin = FindFakeElectronResponseBin( Pt , ptBins, " << NPtBins << ");" << endl;
  outf_e << "  Int_t tmpEtaBin = FindFakeElectronResponseBin( Eta , etaBins, " << NEtaBins << ");" << endl;
  outf_e << "  return ResponseMean[tmpPtBin][tmpEtaBin];" << endl;
  outf_e << "}" << endl;

  outf_e << endl;
  outf_e << endl;

  outf_e << "Double_t GetElectronResponseSigmaPtEta(Double_t Pt, Double_t Eta) {" << endl;

  outf_e << endl;
  outf_e << "  Double_t ptBins[" << NPtBins+1 << "] = {";
  for (uint i=0; i < NPtBins+1; ++i) {
    outf_e << ptBins[i];
    if (i < NPtBins) {
      outf_e << ",";
    }
  }
  outf_e << "};\n";

  outf_e << "  Double_t etaBins[" << NEtaBins+1 << "] = {";
  for (uint i=0; i < NEtaBins+1; ++i) {
    outf_e << etaBins[i];
    if (i < NEtaBins) {
      outf_e << ",";
    }
  }
  outf_e << "};\n";


  outf_e << endl;
  outf_e << endl;

  outf_e << "  Double_t ResponseSigma[" << NPtBins+2 << "][" << NEtaBins+2 << "] = {";
  outf_e << endl;

  for (uint i=0; i < NPtBins+2; ++i) {
    outf_e << "    {";
    for (uint j=0; j < NEtaBins+2; ++j) {
      outf_e << GaussParamArray_Electrons_sigma->GetBinContent(i,j);
      if (j< NEtaBins+1) {
        outf_e << ",";
      }
    }
    if (i< NPtBins+1) {
      outf_e << "    },";
    } else {
      outf_e << "}";
    }
    outf_e << endl;
  }
  
  outf_e << "  };" << endl;

  outf_e << endl;
  outf_e << endl;

  outf_e << "  Int_t tmpPtBin = FindFakeElectronResponseBin( Pt , ptBins, " << NPtBins << ");" << endl;
  outf_e << "  Int_t tmpEtaBin = FindFakeElectronResponseBin( Eta , etaBins, " << NEtaBins << ");" << endl;
  outf_e << "  return ResponseSigma[tmpPtBin][tmpEtaBin];" << endl;
  outf_e << "}" << endl;


  outf_e.close();


  fileInput->Close();
  delete fileInput;
  fileOutput->Close();


  gBenchmark->Show("WWTemplate");       
}


void CreateElectronFakeRateMap(const string filename, const string Label = "ZZ", Int_t Option = 0,
                      Int_t PtBin = -1, Int_t EtaBin = -1) {
 
  if (Option == 0) {
    doCreateElectronFakeRateMap(filename, Label);
  }
  if (Option == 1) {
    FitLeptonResponseModels(Label, Option,PtBin, EtaBin);
  }

} 

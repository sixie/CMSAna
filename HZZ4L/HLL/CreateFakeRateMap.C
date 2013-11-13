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

#endif
 

//=== FUNCTION DECLARATIONS ======================================================================================

void initialize3DArray(double ***array, UInt_t NPtBins, UInt_t NEtaBins, UInt_t NPhiBins) {

  array = new double**[NPtBins+2];
  for (uint i=0; i < NPtBins+2; ++i) {
    array[i] = new double*[NEtaBins+2];
    for (uint j=0; j < NEtaBins+2; ++j) {
      array[i][j] = new double[NPhiBins+2];
      for (uint k=0; k < NPhiBins+2; ++k) {
        array[i][j][k] = 0;
      }      
    }
  }
}

void initialize3DArray(vector<vector<vector<double> > > &array, UInt_t NPtBins, UInt_t NEtaBins, UInt_t NPhiBins) {

  array.resize(NPtBins+2);
  for (uint i=0; i < NPtBins+2; ++i) {
    array[i].resize(NEtaBins+2);
    for (uint j=0; j < NEtaBins+2; ++j) {
      array[i][j].resize(NPhiBins+2);
      for (uint k=0; k < NPhiBins+2; ++k) {
        array[i][j][k] = 0;
      }
    }
  }
}

void initialize2DArray(vector<vector<double> > &array, UInt_t NPtBins, UInt_t NEtaBins) {

  array.resize(NPtBins+2);
  for (uint i=0; i < NPtBins+2; ++i) {
    array[i].resize(NEtaBins+2);
    for (uint j=0; j < NEtaBins+2; ++j) {
      array[i][j]= 0;
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
  
void computeEfficiencyPtEtaPhi(vector<vector<vector<double> > > &numerator, 
                       vector<vector<vector<double> > > &denominator,
                       vector<vector<vector<double> > > &eff
  ) {

  for (uint i=0; i < numerator.size(); ++i) {
    for (uint j=0; j < numerator[i].size(); ++j) {
      for (uint k=0; k < numerator[i][j].size(); ++k) {
        if (denominator[i][j][k] > 0) {         
          eff[i][j][k] = numerator[i][j][k] / denominator[i][j][k];        
        } else {
          eff[i][j][k] = 0;
        }
      }
    }
  }
}

void computeEfficiencyPtEta(vector<vector<double> > &numerator, 
                            vector<vector<double> > &denominator,
                            vector<vector<double> > &eff,
                            int leptonType
  ) {

  for (uint i=0; i < numerator.size(); ++i) {
    for (uint j=0; j < numerator[i].size(); ++j) {
      
      double manualCorrFactor = 1.0;
      if (leptonType == 11) {
        if (i < 4) manualCorrFactor = 1.25;
        if (i >= 4 && i < 6) manualCorrFactor = 1.03;
        if (i >= 6 &&  i < 9) manualCorrFactor = 1.02;
        if (i >= 9 &&  i < 11) manualCorrFactor = 1.01;
        if (i >= 11 &&  i < 15) manualCorrFactor = 0.975;
        if (i >= 15 ) manualCorrFactor = 0.97;
        
      } else if (leptonType == 13) {
        if (i < 2) manualCorrFactor *= 1.15;
        if (i >= 2 && i < 3) manualCorrFactor *= 1.05;
        if (i >= 3 && i < 6) manualCorrFactor *= 1.015;
        if (i >= 6 &&  i < 11) manualCorrFactor *= 1.015;
        if (i >= 11 &&  i < 15) manualCorrFactor *= 1.00;
        if (i >= 15 ) manualCorrFactor *= 0.985;
      }
      
      if ( denominator[i][j] > 0) {
        eff[i][j] = manualCorrFactor * (numerator[i][j] / denominator[i][j]);
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

void CreateFakeRateMap(const string filename, const string Label = "ZZ", Int_t Option = 11) 
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
   const UInt_t NPtBins = 7; 
   const UInt_t NEtaBins = 3;
   double ptBins[NPtBins+1] = { 5, 7, 10, 15, 20, 25, 30, 35 };
   double etaBins[NEtaBins+1] = { 0.0, 1.0, 2.0, 2.5 };

   vector<vector<double> > NDenominator_Electrons_PtEta;
   vector<vector<double> > NNumerator_Electrons_PtEta;
   vector<vector<double> > Efficiency_Electrons_PtEta;
   vector<vector<double> > NDenominator_Muons_PtEta;
   vector<vector<double> > NNumerator_Muons_PtEta;
   vector<vector<double> > Efficiency_Muons_PtEta;
   
   initialize2DArray(NDenominator_Electrons_PtEta, NPtBins, NEtaBins);
   initialize2DArray(NNumerator_Electrons_PtEta, NPtBins, NEtaBins);
   initialize2DArray(Efficiency_Electrons_PtEta, NPtBins, NEtaBins);
   initialize2DArray(NDenominator_Muons_PtEta, NPtBins, NEtaBins);
   initialize2DArray(NNumerator_Muons_PtEta, NPtBins, NEtaBins);
   initialize2DArray(Efficiency_Muons_PtEta, NPtBins, NEtaBins);
  
   
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

  file->Close();
  delete file;       




   computeEfficiencyPtEta(NNumerator_Electrons_PtEta, NDenominator_Electrons_PtEta, Efficiency_Electrons_PtEta,11);
//   computeEfficiencyPtEta(NNumerator_Muons_PtEta, NDenominator_Muons_PtEta, Efficiency_Muons_PtEta,13);


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




//   ofstream outf_m("MuonEfficiencyMap.h");

//   outf_m << "UInt_t FindMuonEfficiencyBin( double value, double bins[], UInt_t nbins) {" << endl;
//   outf_m << "  UInt_t nbinboundaries = nbins+1;" << endl;
//   outf_m << "  UInt_t bin = 0;" << endl;
//   outf_m << "  for (uint i=0; i < nbinboundaries; ++i) {" << endl;
//   outf_m << "    if (i < nbinboundaries-1) {" << endl;
//   outf_m << "      if (value >= bins[i] && value < bins[i+1]) {" << endl;
//   outf_m << "        bin = i+1;" << endl;
//   outf_m << "        break;" << endl;
//   outf_m << "      }" << endl;
//   outf_m << "    } else if (i == nbinboundaries-1) {" << endl;
//   outf_m << "      if (value >= bins[i]) {" << endl;
//   outf_m << "        bin = nbinboundaries;" << endl;
//   outf_m << "        break;" << endl;
//   outf_m << "      }" << endl;
//   outf_m << "    }    " << endl;
//   outf_m << "  }" << endl;
//   outf_m << "  return bin;" << endl;
//   outf_m << "}" << endl;

//   outf_m << endl;
//   outf_m << endl;

//   outf_m << "Double_t GetMuonEfficiencyPtEtaPhi(Double_t Pt, Double_t Eta, Double_t Phi) {" << endl;

//   outf_m << endl;
//   outf_m << "  Double_t ptBins[" << NPtBins+1 << "] = {";
//   for (uint i=0; i < NPtBins+1; ++i) {
//     outf_m << ptBins[i];
//     if (i < NPtBins) {
//       outf_m << ",";
//     }
//   }
//   outf_m << "};\n";

//   outf_m << "  Double_t etaBins[" << NEtaBins+1 << "] = {";
//   for (uint i=0; i < NEtaBins+1; ++i) {
//     outf_m << etaBins[i];
//     if (i < NEtaBins) {
//       outf_m << ",";
//     }
//   }
//   outf_m << "};\n";

//   outf_m << "  Double_t phiBins[" << NPhiBins+1 << "] = {";
//   for (uint i=0; i < NPhiBins+1; ++i) {
//     outf_m << phiBins[i];
//     if (i < NPhiBins) {
//       outf_m << ",";
//     }
//   }
//   outf_m << "};\n";

//   outf_m << endl;
//   outf_m << endl;

//   outf_m << "  Double_t Efficiency[" << NPtBins+2 << "][" << NEtaBins+2 << "][" << NPhiBins+2 << "]  = {";
//   outf_m << endl;

//   for (uint i=0; i < NPtBins+2; ++i) {
//     outf_m << "    {" << endl;
//     for (uint j=0; j < NEtaBins+2; ++j) {
//       outf_m << "      {";
//       for (uint k=0; k < NPhiBins+2; ++k) {
//         outf_m << Efficiency_Muons_PtEtaPhi[i][j][k];
//         if (k< NPhiBins+1) {
//           outf_m << ",";
//         }        
//       }
//       if (j< NEtaBins+1) {
//         outf_m << "},";
//       } else {
//         outf_m << "}";
//       }
//       outf_m << endl;
//     }
//     if (i< NPtBins+1) {
//       outf_m << "    },";
//     } else {
//       outf_m << "}";
//     }
//     outf_m << endl;
//   }

//   outf_m << "};" << endl;

//   outf_m << endl;
//   outf_m << endl;

//   outf_m << "  Int_t tmpPtBin = FindMuonEfficiencyBin( Pt , ptBins, " << NPtBins << ");" << endl;
//   outf_m << "  Int_t tmpEtaBin = FindMuonEfficiencyBin( Eta , etaBins, " << NEtaBins << ");" << endl;
//   outf_m << "  Int_t tmpPhiBin = FindMuonEfficiencyBin( Phi , phiBins, " << NPhiBins << ");" << endl;
//   outf_m << "  return Efficiency[tmpPtBin][tmpEtaBin][tmpPhiBin];" << endl;
//   outf_m << "}" << endl;


//   outf_m << endl;
//   outf_m << endl;
//   outf_m << endl;
//   outf_m << endl;


//   outf_m << "Double_t GetMuonEfficiencyPtEta(Double_t Pt, Double_t Eta) {" << endl;

//   outf_m << endl;
//   outf_m << "  Double_t ptBins[" << NPtBins+1 << "] = {";
//   for (uint i=0; i < NPtBins+1; ++i) {
//     outf_m << ptBins[i];
//     if (i < NPtBins) {
//       outf_m << ",";
//     }
//   }
//   outf_m << "};\n";

//   outf_m << "  Double_t etaBins[" << NEtaBins+1 << "] = {";
//   for (uint i=0; i < NEtaBins+1; ++i) {
//     outf_m << etaBins[i];
//     if (i < NEtaBins) {
//       outf_m << ",";
//     }
//   }
//   outf_m << "};\n";


//   outf_m << endl;
//   outf_m << endl;

//   outf_m << "  Double_t Efficiency[" << NPtBins+2 << "][" << NEtaBins+2 << "] = {";
//   outf_m << endl;

//   for (uint i=0; i < NPtBins+2; ++i) {
//     outf_m << "    {";
//     for (uint j=0; j < NEtaBins+2; ++j) {
//       outf_m << Efficiency_Muons_PtEta[i][j];
//       if (j< NEtaBins+1) {
//         outf_m << ",";
//       }
//     }
//     if (i< NPtBins+1) {
//       outf_m << "    },";
//     } else {
//       outf_m << "}";
//     }
//     outf_m << endl;
//   }
  
//   outf_m << "  };" << endl;

//   outf_m << endl;
//   outf_m << endl;

//   outf_m << "  Int_t tmpPtBin = FindMuonEfficiencyBin( Pt , ptBins, " << NPtBins << ");" << endl;
//   outf_m << "  Int_t tmpEtaBin = FindMuonEfficiencyBin( Eta , etaBins, " << NEtaBins << ");" << endl;
//   outf_m << "  return Efficiency[tmpPtBin][tmpEtaBin];" << endl;
//   outf_m << "}" << endl;


//   outf_m.close();




  gBenchmark->Show("WWTemplate");       
} 



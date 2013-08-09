#ifndef __CINT__
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>
#endif

#include <math.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include "CMSAna/HHToBBGG/interface/HHToBBGGEventTree.hh"
#include <vector>
#include <map>
#include "CMSAna/Utils/CommonTools.hh"

Int_t bkgColors[10] = { kRed, kBlue,  kGreen+2, kMagenta, kCyan, kGreen+3, kBlack, kRed-6, kBlue-6, kMagenta+3 };
string bkgLegendLabels[10] = { "HH#rightarrow bb#gamma#gamma", 
                             "ttH,H#rightarrow#gamma#gamma",
                             "ZH#rightarrow bb #gamma#gamma",
                             "bbH (ggH)",
                             "ttbar",
                             "#gamma#gamma + bjets",
			      "#gamma#gamma + 2 mistag",
			      "bb + 2 fake #gamma#gamma",
			      "cc + 2 fake #gamma#gamma",
			      "jjjj"};
string LatexLabels[10] = { "$HH\\rightarrow bb\\gamma\\gamma$", 
                             "$ttH,H\\rightarrow\\gamma\\gamma$",
                             "$ZH\\rightarrow bb \\gamma\\gamma$",
                             "bbH (ggH)",
                             "ttbar",
                             "$bbgg$",
			      "$jjgg$",
			      "$bbjj$",
			      "$ccjj$",
			      "$jjjj$"};

//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
TH1F* NormalizeHist(TH1F *originalHist) {
  TH1F* hist = (TH1F*)originalHist->Clone((string(originalHist->GetName())+"_normalized").c_str());
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return hist;
}


//------------------------------------------------------------------------------
// PlotHiggsRes_LP
//------------------------------------------------------------------------------
void DiHiggsBackgroundCuts( string InputFilename    = "/afs/cern.ch/work/v/vlambert/public/HHbbggBackground/new_normalizedNtuples/signal_bkgd.root", string outputfile = "out.txt") {
  
  double intLumi = 3000; //in units of fb^-1
  ofstream outfile;
  outfile.open(outputfile.c_str());
  
  outfile << "\\begin{multicols}{2}"<<endl;
  outfile << "\\begin{center}"<<endl;
  outfile << "\\begin{tabular}{| C{3.0cm}|| C{3.0cm}|}" <<endl;
  outfile << "\\hline" << endl;
  outfile << "\\multicolumn{2}{|c|}{\\textbf{Stage 0}}\\\\ \\hline" <<endl;
  outfile << "Sample Type & Expected Events \\\\ \\hline"<< endl;
  
  
  //*******************************************************************************************
  //Define Histograms
  //*******************************************************************************************
  //stage 0
  vector<TH1F*> dRgg_ps;
  vector<TH1F*> minDRgb_ps;
  vector<TH1F*> diphotonMass_s0;
  vector<TH1F*> dibjetMass_s0;

  //stage 1
  vector<TH1F*> dRbb_s0;
  vector<TH1F*> minDRgb_s0;
  vector<TH1F*> diphotonMass_s1;
  vector<TH1F*> dibjetMass_s1;

  //stage 2
  vector<TH1F*> SumbbggPt_ps;
  vector<TH1F*> bbggPt_ps;
  vector<TH1F*> diphotonPt_ps;
  vector<TH1F*> dibjetPt_ps;
  vector<TH1F*> diphotonMass_s2;
  vector<TH1F*> dibjetMass_s2;

  // yields
  vector<TH1F*> bbggMass_s0;
  vector<TH1F*> bbggMass_s1;
  vector<TH1F*> bbggMass_s1t;
  vector<TH1F*> bbggMass_s2;
  vector<TH1F*> bbggMass_s2t;

  vector<TH1F*> bbggMass_s0_120;
  vector<TH1F*> bbggMass_s1_120;
  vector<TH1F*> bbggMass_s1t_120;
  vector<TH1F*> bbggMass_s2_120;
  vector<TH1F*> bbggMass_s2t_120;

  vector<TH1F*> bbggMass_s0_122;
  vector<TH1F*> bbggMass_s1_122;
  vector<TH1F*> bbggMass_s1t_122;
  vector<TH1F*> bbggMass_s2_122;
  vector<TH1F*> bbggMass_s2t_122;

    
  for (UInt_t i = 1; i <= 10; ++i) {

    dRgg_ps.push_back( new TH1F( Form("dRgg_%d",i), ";#Delta R(#gamma,#gamma);Number of Events", 50, 0, 5));
    dRgg_ps[i-1]->SetFillColor(bkgColors[i-1]);
    dRgg_ps[i-1]->SetLineColor(bkgColors[i-1]);
    dRgg_ps[i-1]->SetStats(false);

    dRbb_s0.push_back( new TH1F( Form("dRbb_%d",i), ";#Delta R(b,b);Number of Events", 25, 0, 5));
    dRbb_s0[i-1]->SetFillColor(bkgColors[i-1]);
    dRbb_s0[i-1]->SetLineColor(bkgColors[i-1]);
    dRbb_s0[i-1]->SetStats(false);
    
    minDRgb_ps.push_back( new TH1F( Form("minDRgb_PS_%d",i), ";min #Delta R(#gamma,b);Number of Events", 40, 0, 4));
    minDRgb_ps[i-1]->SetFillColor(bkgColors[i-1]);
    minDRgb_ps[i-1]->SetLineColor(bkgColors[i-1]);
    minDRgb_ps[i-1]->SetStats(false);

    minDRgb_s0.push_back( new TH1F( Form("minDRgbS0_%d",i), ";min #Delta R(#gamma,b);Number of Events", 40, 0, 4));
    minDRgb_s0[i-1]->SetFillColor(bkgColors[i-1]);
    minDRgb_s0[i-1]->SetLineColor(bkgColors[i-1]);
    minDRgb_s0[i-1]->SetStats(false);



    diphotonPt_ps.push_back( new TH1F( Form("diphotonPt_%d",i), "; p_{T #gamma#gamma}[GeV/c];Number of Events", 30, 50, 200));
    diphotonPt_ps[i-1]->SetFillColor(bkgColors[i-1]);
    diphotonPt_ps[i-1]->SetLineColor(bkgColors[i-1]);
    diphotonPt_ps[i-1]->SetStats(false);
    
    dibjetPt_ps.push_back( new TH1F( Form("dibjetPt_%d",i), "; p_{T bb}[GeV/c];Number of Events", 25, 0, 400));
    dibjetPt_ps[i-1]->SetFillColor(bkgColors[i-1]);
    dibjetPt_ps[i-1]->SetLineColor(bkgColors[i-1]);
    dibjetPt_ps[i-1]->SetStats(false);

    SumbbggPt_ps.push_back( new TH1F( Form("bb+ggPt_%d",i), "; p_{T bb}+p_{T#gamma#gamma} [GeV/c];Number of Events", 50, 25, 1000));
    SumbbggPt_ps[i-1]->SetFillColor(bkgColors[i-1]);
    SumbbggPt_ps[i-1]->SetLineColor(bkgColors[i-1]);
    SumbbggPt_ps[i-1]->SetStats(false); 

    bbggPt_ps.push_back( new TH1F( Form("bbggPt_%d",i), "; p_{T bb+#gamma#gamma} [GeV/c];Number of Events", 50, 0, 200));
    bbggPt_ps[i-1]->SetFillColor(bkgColors[i-1]);
    bbggPt_ps[i-1]->SetLineColor(bkgColors[i-1]);
    bbggPt_ps[i-1]->SetStats(false); 

    diphotonMass_s0.push_back( new TH1F( Form("diphotonMass0_%d",i), ";M_{#gamma#gamma} [GeV/cˆ{2}];Number of Events", 15, 100, 150));
    diphotonMass_s0[i-1]->SetFillColor(bkgColors[i-1]);
    diphotonMass_s0[i-1]->SetLineColor(bkgColors[i-1]);
    diphotonMass_s0[i-1]->SetStats(false);

   dibjetMass_s0.push_back( new TH1F( Form("dibjetMass0_%d",i), ";M_{bb} [GeV/c^{2}];Number of Events",20 , 70, 200));
    dibjetMass_s0[i-1]->SetFillColor(bkgColors[i-1]);
    dibjetMass_s0[i-1]->SetLineColor(bkgColors[i-1]);
    dibjetMass_s0[i-1]->SetStats(false);



    diphotonMass_s1.push_back( new TH1F( Form("diphotonMass1_%d",i), ";M_{#gamma#gamma} [GeV/cˆ{2}];Number of Events", 15, 100, 150));
    diphotonMass_s1[i-1]->SetFillColor(bkgColors[i-1]);
    diphotonMass_s1[i-1]->SetLineColor(bkgColors[i-1]);
    diphotonMass_s1[i-1]->SetStats(false);

   dibjetMass_s1.push_back( new TH1F( Form("dibjetMass1_%d",i), ";M_{bb} [GeV/c^{2}];Number of Events", 20, 70, 200));
    dibjetMass_s1[i-1]->SetFillColor(bkgColors[i-1]);
    dibjetMass_s1[i-1]->SetLineColor(bkgColors[i-1]);
    dibjetMass_s1[i-1]->SetStats(false);



    diphotonMass_s2.push_back( new TH1F( Form("diphotonMass2_%d",i), ";M_{#gamma#gamma} [GeV/cˆ{2}];Number of Events", 15, 100, 150));
    diphotonMass_s2[i-1]->SetFillColor(bkgColors[i-1]);
    diphotonMass_s2[i-1]->SetLineColor(bkgColors[i-1]);
    diphotonMass_s2[i-1]->SetStats(false);

   dibjetMass_s2.push_back( new TH1F( Form("dibjetMass2_%d",i), ";M_{bb} [GeV/c^{2}];Number of Events", 20, 70, 200));
    dibjetMass_s2[i-1]->SetFillColor(bkgColors[i-1]);
    dibjetMass_s2[i-1]->SetLineColor(bkgColors[i-1]);
    dibjetMass_s2[i-1]->SetStats(false);
    
    

    
    bbggMass_s0.push_back( new TH1F( Form("bbggMass0_%d",i), "; m_{bb#gamma#gamma} [GeV/c^{2}];Number of Events", 100, 0, 1000));
    bbggMass_s0[i-1]->SetFillColor(bkgColors[i-1]);
    bbggMass_s0[i-1]->SetLineColor(bkgColors[i-1]);
    bbggMass_s0[i-1]->SetStats(false);
    
    bbggMass_s1.push_back( new TH1F( Form("bbggMass1_%d",i), "; m_{bb#gamma#gamma} [GeV/c^{2}];Number of Events", 100, 0, 1000));
    bbggMass_s1[i-1]->SetFillColor(bkgColors[i-1]);
    bbggMass_s1[i-1]->SetLineColor(bkgColors[i-1]);
    bbggMass_s1[i-1]->SetStats(false);

    bbggMass_s1t.push_back( new TH1F( Form("bbggMass1t_%d",i), "; m_{bb#gamma#gamma} [GeV/c^{2}];Number of Events", 100, 0, 1000));
    bbggMass_s1t[i-1]->SetFillColor(bkgColors[i-1]);
    bbggMass_s1t[i-1]->SetLineColor(bkgColors[i-1]);
    bbggMass_s1t[i-1]->SetStats(false);
    
    bbggMass_s2.push_back( new TH1F( Form("bbggMass2_%d",i), "; m_{bb#gamma#gamma} [GeV/c^{2}];Number of Events", 100, 0, 1000));
    bbggMass_s2[i-1]->SetFillColor(bkgColors[i-1]);
    bbggMass_s2[i-1]->SetLineColor(bkgColors[i-1]);
    bbggMass_s2[i-1]->SetStats(false);

    bbggMass_s2t.push_back( new TH1F( Form("bbggMass2t_%d",i), "; m_{bb#gamma#gamma} [GeV/c^{2}];Number of Events", 100, 0, 1000));
    bbggMass_s2t[i-1]->SetFillColor(bkgColors[i-1]);
    bbggMass_s2t[i-1]->SetLineColor(bkgColors[i-1]);
    bbggMass_s2t[i-1]->SetStats(false);



    bbggMass_s0_120.push_back( new TH1F( Form("bbggMass0_120_%d",i), "; m_{bb#gamma#gamma} [GeV/c^{2}];Number of Events", 100, 0, 1000));
    bbggMass_s0_120[i-1]->SetFillColor(bkgColors[i-1]);
    bbggMass_s0_120[i-1]->SetLineColor(bkgColors[i-1]);
    bbggMass_s0_120[i-1]->SetStats(false);
    
    bbggMass_s1_120.push_back( new TH1F( Form("bbggMass1_120_%d",i), "; m_{bb#gamma#gamma} [GeV/c^{2}];Number of Events", 100, 0, 1000));
    bbggMass_s1_120[i-1]->SetFillColor(bkgColors[i-1]);
    bbggMass_s1_120[i-1]->SetLineColor(bkgColors[i-1]);
    bbggMass_s1_120[i-1]->SetStats(false);

    bbggMass_s1t_120.push_back( new TH1F( Form("bbggMass1t_120_%d",i), "; m_{bb#gamma#gamma} [GeV/c^{2}];Number of Events", 100, 0, 1000));
    bbggMass_s1t_120[i-1]->SetFillColor(bkgColors[i-1]);
    bbggMass_s1t_120[i-1]->SetLineColor(bkgColors[i-1]);
    bbggMass_s1t_120[i-1]->SetStats(false);
    
    bbggMass_s2_120.push_back( new TH1F( Form("bbggMass2_120_%d",i), "; m_{bb#gamma#gamma} [GeV/c^{2}];Number of Events", 100, 0, 1000));
    bbggMass_s2_120[i-1]->SetFillColor(bkgColors[i-1]);
    bbggMass_s2_120[i-1]->SetLineColor(bkgColors[i-1]);
    bbggMass_s2_120[i-1]->SetStats(false);

    bbggMass_s2t_120.push_back( new TH1F( Form("bbggMass2t_120_%d",i), "; m_{bb#gamma#gamma} [GeV/c^{2}];Number of Events", 100, 0, 1000));
    bbggMass_s2t_120[i-1]->SetFillColor(bkgColors[i-1]);
    bbggMass_s2t_120[i-1]->SetLineColor(bkgColors[i-1]);
    bbggMass_s2t_120[i-1]->SetStats(false);


    bbggMass_s0_122.push_back( new TH1F( Form("bbggMass0_122_%d",i), "; m_{bb#gamma#gamma} [GeV/c^{2}];Number of Events", 100, 0, 1000));
    bbggMass_s0_122[i-1]->SetFillColor(bkgColors[i-1]);
    bbggMass_s0_122[i-1]->SetLineColor(bkgColors[i-1]);
    bbggMass_s0_122[i-1]->SetStats(false);
    
    bbggMass_s1_122.push_back( new TH1F( Form("bbggMass1_122_%d",i), "; m_{bb#gamma#gamma} [GeV/c^{2}];Number of Events", 100, 0, 1000));
    bbggMass_s1_122[i-1]->SetFillColor(bkgColors[i-1]);
    bbggMass_s1_122[i-1]->SetLineColor(bkgColors[i-1]);
    bbggMass_s1_122[i-1]->SetStats(false);

    bbggMass_s1t_122.push_back( new TH1F( Form("bbggMass1_122t_%d",i), "; m_{bb#gamma#gamma} [GeV/c^{2}];Number of Events", 100, 0, 1000));
    bbggMass_s1t_122[i-1]->SetFillColor(bkgColors[i-1]);
    bbggMass_s1t_122[i-1]->SetLineColor(bkgColors[i-1]);
    bbggMass_s1t_122[i-1]->SetStats(false);
    
    bbggMass_s2_122.push_back( new TH1F( Form("bbggMass2_122_%d",i), "; m_{bb#gamma#gamma} [GeV/c^{2}];Number of Events", 100, 0, 1000));
    bbggMass_s2_122[i-1]->SetFillColor(bkgColors[i-1]);
    bbggMass_s2_122[i-1]->SetLineColor(bkgColors[i-1]);
    bbggMass_s2_122[i-1]->SetStats(false);

    bbggMass_s2t_122.push_back( new TH1F( Form("bbggMass2t_122_%d",i), "; m_{bb#gamma#gamma} [GeV/c^{2}];Number of Events", 100, 0, 1000));
    bbggMass_s2t_122[i-1]->SetFillColor(bkgColors[i-1]);
    bbggMass_s2t_122[i-1]->SetLineColor(bkgColors[i-1]);
    bbggMass_s2t_122[i-1]->SetStats(false);


  }
  
  
  //*******************************************************************************************
  //Do analysis
  //*******************************************************************************************
  
  cmsana::HHToBBGGEventTree event;
  event.LoadTree(InputFilename.c_str());
  event.InitTree();
  
  //Double_t deltaR_bbgg;
  //Double_t deltaPhi_bbgg;
  
  Int_t Preselection = 0;
  Int_t Stage0 = 0;
  Int_t Stage1 = 0;
  Int_t Stage1tight = 0;
  Int_t Stage2 = 0;
  Int_t Stage2tight = 0;
  Int_t tightwindow = 0;
  Int_t optwindow = 0;

  for (int n=0;n<event.tree_->GetEntries();n++) { 
    event.tree_->GetEntry(n);
    double weight = event.weight*intLumi;

    Preselection = 0;
    Stage0 = 0;
    Stage1 = 0;
    Stage1tight = 0;
    Stage2 = 0;
    Stage2tight = 0;
    tightwindow = 0;
    optwindow = 0;
    
    
    //make sure the sample type makes sense
    assert(event.sampletype >= 1 && event.sampletype <= 13);
    
    if (event.bjet1.Pt() > 25.0 && event.bjet2.Pt() > 25.0
	&& event.pho1.Pt() > 25.0 && event.pho2.Pt() > 25.0
	&& fmax(event.pho1.Pt(),event.pho2.Pt()) > 40.0
	&& event.nlep <=0 && event.ncentraljets <= 3
	&& event.diphoton.M() > 100 && event.diphoton.M() < 150
	&& event.dibjet.M() > 70 && event.dibjet.M() < 200
	) {
      Preselection = 1;
    }
    
    if (2.0 > event.DRgg && event.minDRgb > 1.0) {
      Stage0 = 1;
    }
    
    if (2.0 > event.DRgg && 2.0 > event.DRbb && event.minDRgb > 1.5) {
      Stage1 = 1;
    }

    if (1.6 > event.DRgg && 1.6 > event.DRbb && event.minDRgb > 1.5) {
      Stage1tight = 1;
    }
    
    if (event.bbgg.Pt() <110 && event.bbgg.Pt() > 10
	&& (event.diphoton.Pt() + event.dibjet.Pt()) > 260
	) {
      Stage2 = 1;
    }

    if (event.bbgg.Pt() <85 && event.bbgg.Pt() > 10
	&& (event.diphoton.Pt() + event.dibjet.Pt()) > 340
	) {
      Stage2tight = 1;
    }
   
    if (event.diphoton.M() > 120 && event.diphoton.M() < 130
	&& event.dibjet.M() > 105 && event.dibjet.M() < 145) {
      tightwindow = 1;      
    }

    if (event.diphoton.M() > 122 && event.diphoton.M() < 128
	&& event.dibjet.M() > 110 && event.dibjet.M() < 140) {
      optwindow = 1;      
    }
	
      
    //Post Selection
    if (Preselection == 1) {
      dRgg_ps[event.sampletype-1]->Fill( event.DRgg , weight);
      minDRgb_ps[event.sampletype-1]->Fill( event.minDRgb , weight);
      SumbbggPt_ps[event.sampletype-1]->Fill( event.dibjet.Pt() + event.diphoton.Pt() , weight);
      bbggPt_ps[event.sampletype-1]->Fill( event.bbgg.Pt() , weight);
      diphotonPt_ps[event.sampletype-1]->Fill( event.diphoton.Pt() , weight);
      dibjetPt_ps[event.sampletype-1]->Fill( event.dibjet.Pt() , weight);
    }
    
    // Stage 0
    if (  Preselection == 1 && Stage0 == 1) {
      diphotonMass_s0[event.sampletype-1]->Fill( event.diphoton.M() , weight);
      dibjetMass_s0[event.sampletype-1]->Fill( event.dibjet.M() , weight);
      dRbb_s0[event.sampletype-1]->Fill( event.DRbb , weight);
      minDRgb_s0[event.sampletype-1]->Fill( event.minDRgb , weight);
      bbggMass_s0[event.sampletype-1]->Fill( event.bbgg.M() , weight);
      if (tightwindow == 1) {
	bbggMass_s0_120[event.sampletype-1]->Fill( event.bbgg.M() , weight);
      }
      if (optwindow == 1) {
	bbggMass_s0_122[event.sampletype-1]->Fill( event.bbgg.M() , weight);
      }
    }
  
    // Stage 1
    if (  Preselection == 1 && Stage1 == 1) {
      diphotonMass_s1[event.sampletype-1]->Fill( event.diphoton.M() , weight);
      dibjetMass_s1[event.sampletype-1]->Fill( event.dibjet.M() , weight);
      bbggMass_s1[event.sampletype-1]->Fill( event.bbgg.M() , weight);
      if (tightwindow == 1) {
	bbggMass_s1_120[event.sampletype-1]->Fill( event.bbgg.M() , weight);
      }
      if (optwindow == 1) {
	bbggMass_s1_122[event.sampletype-1]->Fill( event.bbgg.M() , weight);
      }
    }
    
    // Stage 1t
    if (  Preselection == 1 && Stage1tight == 1) {
      bbggMass_s1t[event.sampletype-1]->Fill( event.bbgg.M() , weight);
      if (tightwindow == 1) {
	bbggMass_s1t_120[event.sampletype-1]->Fill( event.bbgg.M() , weight);
      }
      if (optwindow == 1) {
	bbggMass_s1t_122[event.sampletype-1]->Fill( event.bbgg.M() , weight);
      }
    }

    // Stage 2
    if (  Preselection == 1 && Stage2 == 1) {
      diphotonMass_s2[event.sampletype-1]->Fill( event.diphoton.M() , weight);
      dibjetMass_s2[event.sampletype-1]->Fill( event.dibjet.M() , weight);
      bbggMass_s2[event.sampletype-1]->Fill( event.bbgg.M() , weight);
      if (tightwindow == 1) {
	bbggMass_s2_120[event.sampletype-1]->Fill( event.bbgg.M() , weight);
      }
      if (optwindow == 1) {
	bbggMass_s2_122[event.sampletype-1]->Fill( event.bbgg.M() , weight);
      }
    }
    
    // Stage 2t
    if (  Preselection == 1 && Stage2tight == 1) {  
      bbggMass_s2t[event.sampletype-1]->Fill( event.bbgg.M() , weight);
      if (tightwindow == 1) {
	bbggMass_s2t_120[event.sampletype-1]->Fill( event.bbgg.M() , weight);
      }
      if (optwindow == 1) {
	bbggMass_s2t_122[event.sampletype-1]->Fill( event.bbgg.M() , weight);
      }
    }
  }

  //Draw Plots
  TCanvas *cv = 0;
  TLegend *legend = 0;
  bool firstdrawn = false;



  //*******************************************************************************************
  //Diphoton mass normalized Stage 0
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.6,0.6,0.80,0.80);
  legend->SetTextSize(0.05);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = diphotonMass_s0.size()-1; i >= 0; i--) {
    if (!(i==5)) continue;
    if (diphotonMass_s0[i]->Integral() > 0) {
      TH1F* normHist = NormalizeHist(diphotonMass_s0[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  legend->Draw();
  cv->SaveAs("MassGG_s0.pdf");


  //*******************************************************************************************
  //Dibjet mass normalized Stage 0
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.6,0.6,0.8,0.8);
  legend->SetTextSize(0.05);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = dibjetMass_s0.size()-1; i >= 0; i--) {
    if (!(i==5)) continue;
    if (dibjetMass_s0[i]->Integral() > 0) {
      TH1F* normHist = NormalizeHist(dibjetMass_s0[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  legend->Draw();
  cv->SaveAs("MassBB_s0.pdf");

  //*******************************************************************************************
  //Diphoton mass normalized Stage 1
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.6,0.6,0.8,0.8);
  legend->SetTextSize(0.05);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = diphotonMass_s1.size()-1; i >= 0; i--) {
    if (!(i==5)) continue;
    if (diphotonMass_s1[i]->Integral() > 0) {
      TH1F* normHist = NormalizeHist(diphotonMass_s1[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  legend->Draw();
  cv->SaveAs("MassGG_s1.pdf");


  //*******************************************************************************************
  //Dibjet mass normalized Stage 1
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.6,0.6,0.8,0.8);
  legend->SetTextSize(0.05);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = dibjetMass_s1.size()-1; i >= 0; i--) {
    if (!(i==5)) continue;
    if (dibjetMass_s1[i]->Integral() > 0) {
      TH1F* normHist = NormalizeHist(dibjetMass_s1[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  legend->Draw();
  cv->SaveAs("MassBB_s1.pdf");

  //*******************************************************************************************
  //Diphoton mass normalized Stage 2
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.6,0.6,0.8,0.8);
  legend->SetTextSize(0.05);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = diphotonMass_s2.size()-1; i >= 0; i--) {
    if (!(i==5)) continue;
    if (diphotonMass_s2[i]->Integral() > 0) {
      TH1F* normHist = NormalizeHist(diphotonMass_s2[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  legend->Draw();
  cv->SaveAs("MassGG_s2.pdf");


  //*******************************************************************************************
  //Dibjet mass normalized Stage 2
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.6,0.6,0.8,0.8);
  legend->SetTextSize(0.05);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = dibjetMass_s2.size()-1; i >= 0; i--) {
    if (!(i==5)) continue;
    if (dibjetMass_s2[i]->Integral() > 0) {
      TH1F* normHist = NormalizeHist(dibjetMass_s2[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  legend->Draw();
  cv->SaveAs("MassBB_s2.pdf");
 


  //*******************************************************************************************
  //dRgg Post Selection normalized
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.60,0.90,0.90);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = dRgg_ps.size()-1; i >= 0; i--) {
    if (dRgg_ps[i]->Integral() > 0) {
      if (!(i==0 || i==1 || i==5 || i==6)) continue;
      TH1F* normHist = NormalizeHist(dRgg_ps[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  legend->Draw();
  cv->SaveAs("dRgg_ps_normalized.pdf");


  //*******************************************************************************************
  //dRbb normalized
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.60,0.90,0.90);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = dRbb_s0.size()-1; i >= 0; i--) {
    if (dRbb_s0[i]->Integral() > 0) {
      if (!(i==0 || i==1 || i==5 || i==6)) continue;
      TH1F* normHist = NormalizeHist(dRbb_s0[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  legend->Draw();
  cv->SaveAs("dRbb_s0_normalized.pdf");

  //*******************************************************************************************
  //minDRgb normalized
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.60,0.90,0.90);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = minDRgb_s0.size()-1; i >= 0; i--) {
    if (minDRgb_s0[i]->Integral() > 0) {
      if (!(i==0 || i==1 || i==5 || i==6)) continue;
      TH1F* normHist = NormalizeHist(minDRgb_s0[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  legend->Draw();
  cv->SaveAs("minDRgb_s0_normalized.pdf");

  //*******************************************************************************************
  //minDRgb normalized
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.60,0.90,0.90);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = minDRgb_ps.size()-1; i >= 0; i--) {
    if (minDRgb_ps[i]->Integral() > 0) {
      if (!(i==0 || i==1 || i==5 || i==6)) continue;
      TH1F* normHist = NormalizeHist(minDRgb_ps[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  legend->Draw();
  cv->SaveAs("minDRgb_ps_normalized.pdf");


  //*******************************************************************************************
  //diphoton pt normalized Post Selection
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.6,0.90,0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = diphotonPt_ps.size()-1; i >= 0; i--) {
    if (diphotonPt_ps[i]->Integral() > 0) {
      if (!(i==0 || i==1 || i==5 || i==6)) continue;
      TH1F* normHist = NormalizeHist(diphotonPt_ps[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  legend->Draw();
  cv->SaveAs("diphotonPt_ps_normalized.pdf");


  //*******************************************************************************************
  //dibjet pt normalized Post Selection
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.6,0.90,0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = dibjetPt_ps.size()-1; i >= 0; i--) {
    if (dibjetPt_ps[i]->Integral() > 0) {
      if (!(i==0 || i==1 || i==5 || i==6)) continue;
      TH1F* normHist = NormalizeHist(dibjetPt_ps[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  legend->Draw();
  cv->SaveAs("dibjetPt_ps_normalized.pdf");


  //*******************************************************************************************
  //bbgg pt normalized
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.60,0.90,0.90);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = bbggPt_ps.size()-1; i >= 0; i--) {
    if (bbggPt_ps[i]->Integral() > 0) {
      if (!(i==0 || i==1 || i==5 || i==6)) continue;
      TH1F* normHist = NormalizeHist(bbggPt_ps[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  legend->Draw();
  cv->SaveAs("bbggPt_ps_normalized.pdf");

  //*******************************************************************************************
  //diphoton + dibjet Pt  normalized
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.60,0.90,0.90);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i =SumbbggPt_ps.size()-1; i >= 0; i--) {
    if (SumbbggPt_ps[i]->Integral() > 0) {
      if (!(i==0 || i==1 || i==5 || i==6)) continue;
      TH1F* normHist = NormalizeHist(SumbbggPt_ps[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  legend->Draw();
  cv->SaveAs("SumbbggPt_ps_normalized.pdf");
  

  //*******************************************************************************************
  // M bbgg normalized Stage 0 Wide Window
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.54,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = bbggMass_s0.size()-1; i >= 0; i--) {
    if (bbggMass_s0[i]->Integral() > 0) {
      outfile << Form("%s & %.4f \\\\ \\hline", LatexLabels[i].c_str(), bbggMass_s0[i]->Integral())<<endl;

      TH1F* normHist = NormalizeHist(bbggMass_s0[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  Float_t bkgd1 = 0;
  for (Int_t j = bbggMass_s0.size()-1; j>=1; j--) {
    bkgd1+=bbggMass_s0[j]->Integral();
  }
  outfile << Form("Sig/Bkgd Ratio & %.3f \\\\ \\hline", (bbggMass_s0[0]->Integral())/bkgd1)<<endl;
  outfile << Form("Significance & %.3f \\\\ \\hline", (bbggMass_s0[0]->Integral())/sqrt(bkgd1))<<endl;
  
  legend->Draw();
  cv->SaveAs("bbggMass_s0_normalized.pdf");
  outfile << "\\multicolumn{2}{|c|}{\\textbf{Stage 1}}\\\\ \\hline" <<endl;
  outfile << "Sample Type & Expected Events \\\\ \\hline"<< endl;

  
  //*******************************************************************************************
  // M bbgg normalized Stage 1 Wide Window
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.54,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = bbggMass_s1.size()-1; i >= 0; i--) {
    if (bbggMass_s1[i]->Integral() > 0) {
      outfile << Form("%s & %.4f \\\\ \\hline", LatexLabels[i].c_str(), bbggMass_s1[i]->Integral())<<endl;

      TH1F* normHist = NormalizeHist(bbggMass_s1[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  Float_t bkgd2 = 0;
  for (Int_t j = bbggMass_s1.size()-1; j>=1; j--) {
    bkgd2+=bbggMass_s1[j]->Integral();
  }
  outfile << Form("Sig/Bkgd Ratio & %.3f \\\\ \\hline", (bbggMass_s1[0]->Integral())/bkgd2)<<endl;
  outfile << Form("Significance & %.3f \\\\ \\hline", (bbggMass_s1[0]->Integral())/sqrt(bkgd2))<<endl;
  
  legend->Draw();
  cv->SaveAs("bbggMass_s1_normalized.pdf");


  outfile << "\\multicolumn{2}{|c|}{\\textbf{Stage 1 Tight}}\\\\ \\hline" <<endl;
  outfile << "Sample Type & Expected Events \\\\ \\hline"<< endl;



  //*******************************************************************************************
  // M bbgg normalized Stage 1 Tight Wide Window
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.54,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = bbggMass_s1t.size()-1; i >= 0; i--) {
    if (bbggMass_s1t[i]->Integral() > 0) {
      outfile << Form("%s & %.4f \\\\ \\hline", LatexLabels[i].c_str(), bbggMass_s1t[i]->Integral())<<endl;

      TH1F* normHist = NormalizeHist(bbggMass_s1t[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  Float_t bkgd2t = 0;
  for (Int_t j = bbggMass_s1t.size()-1; j>=1; j--) {
    bkgd2t+=bbggMass_s1t[j]->Integral();
  }
  outfile << Form("Sig/Bkgd Ratio & %.3f \\\\ \\hline", (bbggMass_s1t[0]->Integral())/bkgd2t)<<endl;
  outfile << Form("Significance & %.3f \\\\ \\hline", (bbggMass_s1t[0]->Integral())/sqrt(bkgd2t))<<endl;
  
  legend->Draw();
  cv->SaveAs("bbggMass_s1_normalized.pdf");

  outfile << "\\end{tabular}"<<endl;
  outfile << "\\end{center}" <<endl;

  outfile << "\\begin{center}"<<endl;
  outfile << "\\begin{tabular}{| C{3.0cm}|| C{3.0cm}|}" <<endl;
  outfile << "\\hline" << endl;
  outfile << "\\multicolumn{2}{|c|}{\\textbf{Stage 2}}\\\\ \\hline" <<endl;
  outfile << "Sample Type & Expected Events \\\\ \\hline"<< endl;


  //*******************************************************************************************
  // M bbgg normalized Stage 2 Wide Window
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.54,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = bbggMass_s2.size()-1; i >= 0; i--) {
    if (bbggMass_s2[i]->Integral() > 0) {
      outfile << Form("%s & %.4f \\\\ \\hline", LatexLabels[i].c_str(), bbggMass_s2[i]->Integral())<<endl;

      TH1F* normHist = NormalizeHist(bbggMass_s2[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  Float_t bkgd3 = 0;
  for (Int_t j = bbggMass_s2.size()-1; j>=1; j--) {
    bkgd3+=bbggMass_s2[j]->Integral();
  }
  outfile << Form("Sig/Bkgd Ratio & %.3f \\\\ \\hline", (bbggMass_s2[0]->Integral())/bkgd3)<<endl;
  outfile << Form("Significance & %.3f \\\\ \\hline", (bbggMass_s2[0]->Integral())/sqrt(bkgd3))<<endl;
  
  legend->Draw();
  cv->SaveAs("bbggMass_s2_normalized.pdf");
  outfile << "\\multicolumn{2}{|c|}{\\textbf{Stage 2 Tight}}\\\\ \\hline" <<endl;
  outfile << "Sample Type & Expected Events \\\\ \\hline"<< endl;

  //*******************************************************************************************
  // M bbgg normalized Stage 2 Tight Wide Window
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.54,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = bbggMass_s2t.size()-1; i >= 0; i--) {
    if (bbggMass_s2t[i]->Integral() > 0) {
      outfile << Form("%s & %.4f \\\\ \\hline", LatexLabels[i].c_str(), bbggMass_s2t[i]->Integral())<<endl;

      TH1F* normHist = NormalizeHist(bbggMass_s2t[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  Float_t bkgd4 = 0;
  for (Int_t j = bbggMass_s2t.size()-1; j>=1; j--) {
    bkgd4+=bbggMass_s2t[j]->Integral();
  }
  outfile << Form("Sig/Bkgd Ratio & %.3f \\\\ \\hline", (bbggMass_s2t[0]->Integral())/bkgd4)<<endl;
  outfile << Form("Significance & %.3f \\\\ \\hline", (bbggMass_s2t[0]->Integral())/sqrt(bkgd4))<<endl;
  
  legend->Draw();
  cv->SaveAs("bbggMass_s2t_normalized.pdf");

  outfile << "\\end{tabular}"<<endl;
  outfile << "\\end{center}" <<endl;
  outfile << "\\end{multicols}"<<endl;
  outfile << "                "<<endl;
  outfile << "\\begin{multicols}{2}"<<endl;
  outfile << "\\begin{center}"<<endl;
  outfile << "\\begin{tabular}{| C{3.0cm}|| C{3.0cm}|}" <<endl;
  outfile << "\\hline" << endl;
  outfile << "\\multicolumn{2}{|c|}{\\textbf{Stage 0}}\\\\ \\hline" <<endl;
  outfile << "Sample Type & Expected Events \\\\ \\hline"<< endl;

  //*******************************************************************************************
  // M bbgg normalized Stage 0 Tighter Window
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.54,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = bbggMass_s0_120.size()-1; i >= 0; i--) {
    if (bbggMass_s0_120[i]->Integral() > 0) {
      outfile << Form("%s & %.4f \\\\ \\hline", LatexLabels[i].c_str(), bbggMass_s0_120[i]->Integral())<<endl;

      TH1F* normHist = NormalizeHist(bbggMass_s0_120[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  Float_t bkgd1_120 = 0;
  for (Int_t j = bbggMass_s0_120.size()-1; j>=1; j--) {
    bkgd1_120+=bbggMass_s0_120[j]->Integral();
  }
  outfile << Form("Sig/Bkgd Ratio & %.3f \\\\ \\hline", (bbggMass_s0_120[0]->Integral())/bkgd1_120)<<endl;
  outfile << Form("Significance & %.3f \\\\ \\hline", (bbggMass_s0_120[0]->Integral())/sqrt(bkgd1_120))<<endl;
  
  legend->Draw();
  cv->SaveAs("bbggMass_s0__120normalized.pdf");
  outfile << "\\multicolumn{2}{|c|}{\\textbf{Stage 1}}\\\\ \\hline" <<endl;
  outfile << "Sample Type & Expected Events \\\\ \\hline"<< endl;

  
  //*******************************************************************************************
  // M bbgg normalized Stage 1 Tighter Window
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.54,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = bbggMass_s1_120.size()-1; i >= 0; i--) {
    if (bbggMass_s1_120[i]->Integral() > 0) {
      outfile << Form("%s & %.4f \\\\ \\hline", LatexLabels[i].c_str(), bbggMass_s1_120[i]->Integral())<<endl;

      TH1F* normHist = NormalizeHist(bbggMass_s1_120[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  Float_t bkgd2_120 = 0;
  for (Int_t j = bbggMass_s1_120.size()-1; j>=1; j--) {
    bkgd2_120+=bbggMass_s1_120[j]->Integral();
  }
  outfile << Form("Sig/Bkgd Ratio & %.3f \\\\ \\hline", (bbggMass_s1_120[0]->Integral())/bkgd2_120)<<endl;
  outfile << Form("Significance & %.3f \\\\ \\hline", (bbggMass_s1_120[0]->Integral())/sqrt(bkgd2_120))<<endl;
  
  legend->Draw();
  cv->SaveAs("bbggMass_s1_120_normalized.pdf");

  outfile << "\\end{tabular}"<<endl;
  outfile << "\\end{center}" <<endl;

  outfile << "\\begin{center}"<<endl;
  outfile << "\\begin{tabular}{| C{3.0cm}|| C{3.0cm}|}" <<endl;
  outfile << "\\hline" << endl;
  outfile << "\\multicolumn{2}{|c|}{\\textbf{Stage 1 Tight}}\\\\ \\hline" <<endl;
  outfile << "Sample Type & Expected Events \\\\ \\hline"<< endl;

  //*******************************************************************************************
  // M bbgg normalized Stage 1 Tighter Window
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.54,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = bbggMass_s1t_120.size()-1; i >= 0; i--) {
    if (bbggMass_s1t_120[i]->Integral() > 0) {
      outfile << Form("%s & %.4f \\\\ \\hline", LatexLabels[i].c_str(), bbggMass_s1t_120[i]->Integral())<<endl;

      TH1F* normHist = NormalizeHist(bbggMass_s1t_120[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  Float_t bkgd2t_120 = 0;
  for (Int_t j = bbggMass_s1t_120.size()-1; j>=1; j--) {
    bkgd2t_120+=bbggMass_s1t_120[j]->Integral();
  }
  outfile << Form("Sig/Bkgd Ratio & %.3f \\\\ \\hline", (bbggMass_s1t_120[0]->Integral())/bkgd2t_120)<<endl;
  outfile << Form("Significance & %.3f \\\\ \\hline", (bbggMass_s1t_120[0]->Integral())/sqrt(bkgd2t_120))<<endl;
  
  legend->Draw();
  cv->SaveAs("bbggMass_s1t_120_normalized.pdf");


  outfile << "\\multicolumn{2}{|c|}{\\textbf{Stage 2}}\\\\ \\hline" <<endl;
  outfile << "Sample Type & Expected Events \\\\ \\hline"<< endl;


  //*******************************************************************************************
  // M bbgg normalized Stage 2 Tighter Window
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.54,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = bbggMass_s2_120.size()-1; i >= 0; i--) {
    if (bbggMass_s2_120[i]->Integral() > 0) {
      outfile << Form("%s & %.4f \\\\ \\hline", LatexLabels[i].c_str(), bbggMass_s2_120[i]->Integral())<<endl;

      TH1F* normHist = NormalizeHist(bbggMass_s2_120[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  Float_t bkgd3_120 = 0;
  for (Int_t j = bbggMass_s2_120.size()-1; j>=1; j--) {
    bkgd3_120+=bbggMass_s2_120[j]->Integral();
  }
  outfile << Form("Sig/Bkgd Ratio & %.3f \\\\ \\hline", (bbggMass_s2_120[0]->Integral())/bkgd3_120)<<endl;
  outfile << Form("Significance & %.3f \\\\ \\hline", (bbggMass_s2_120[0]->Integral())/sqrt(bkgd3_120))<<endl;
  
  legend->Draw();
  cv->SaveAs("bbggMass_s2_120_normalized.pdf");
  outfile << "\\multicolumn{2}{|c|}{\\textbf{Stage 2 Tight}}\\\\ \\hline" <<endl;
  outfile << "Sample Type & Expected Events \\\\ \\hline"<< endl;

  //*******************************************************************************************
  // M bbgg normalized Stage 2 Tight Tighter Window
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.54,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = bbggMass_s2t_120.size()-1; i >= 0; i--) {
    if (bbggMass_s2t_120[i]->Integral() > 0) {
      outfile << Form("%s & %.4f \\\\ \\hline", LatexLabels[i].c_str(), bbggMass_s2t_120[i]->Integral())<<endl;

      TH1F* normHist = NormalizeHist(bbggMass_s2t_120[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  Float_t bkgd4_120 = 0;
  for (Int_t j = bbggMass_s2t_120.size()-1; j>=1; j--) {
    bkgd4_120+=bbggMass_s2t_120[j]->Integral();
  }
  outfile << Form("Sig/Bkgd Ratio & %.3f \\\\ \\hline", (bbggMass_s2t_120[0]->Integral())/bkgd4_120)<<endl;
  outfile << Form("Significance & %.3f \\\\ \\hline", (bbggMass_s2t_120[0]->Integral())/sqrt(bkgd4_120))<<endl;
  
  legend->Draw();
  cv->SaveAs("bbggMass_s2t_120_normalized.pdf");

  outfile << "\\end{tabular}"<<endl;
  outfile << "\\end{center}" <<endl;
  outfile << "\\end{multicols}"<<endl;
  outfile << "                "<<endl;
  outfile << "\\begin{multicols}{2}"<<endl;
  outfile << "\\begin{center}"<<endl;
  outfile << "\\begin{tabular}{| C{3.0cm}|| C{3.0cm}|}" <<endl;
  outfile << "\\hline" << endl;
  outfile << "\\multicolumn{2}{|c|}{\\textbf{Stage 0}}\\\\ \\hline" <<endl;
  outfile << "Sample Type & Expected Events \\\\ \\hline"<< endl;


  //*******************************************************************************************
  // M bbgg normalized Stage 0 Optimistic Window
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.54,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = bbggMass_s0_122.size()-1; i >= 0; i--) {
    if (bbggMass_s0_122[i]->Integral() > 0) {
      outfile << Form("%s & %.4f \\\\ \\hline", LatexLabels[i].c_str(), bbggMass_s0_122[i]->Integral())<<endl;

      TH1F* normHist = NormalizeHist(bbggMass_s0_122[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  Float_t bkgd1_122 = 0;
  for (Int_t j = bbggMass_s0_122.size()-1; j>=1; j--) {
    bkgd1_122+=bbggMass_s0_122[j]->Integral();
  }
  outfile << Form("Sig/Bkgd Ratio & %.3f \\\\ \\hline", (bbggMass_s0_122[0]->Integral())/bkgd1_122)<<endl;
  outfile << Form("Significance & %.3f \\\\ \\hline", (bbggMass_s0_122[0]->Integral())/sqrt(bkgd1_122))<<endl;
  
  legend->Draw();
  cv->SaveAs("bbggMass_s0__122normalized.pdf");
  outfile << "\\multicolumn{2}{|c|}{\\textbf{Stage 1}}\\\\ \\hline" <<endl;
  outfile << "Sample Type & Expected Events \\\\ \\hline"<< endl;

  
  //*******************************************************************************************
  // M bbgg normalized Stage 1 Optimistic Window
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.54,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = bbggMass_s1_122.size()-1; i >= 0; i--) {
    if (bbggMass_s1_122[i]->Integral() > 0) {
      outfile << Form("%s & %.4f \\\\ \\hline", LatexLabels[i].c_str(), bbggMass_s1_122[i]->Integral())<<endl;

      TH1F* normHist = NormalizeHist(bbggMass_s1_122[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  Float_t bkgd2_122 = 0;
  for (Int_t j = bbggMass_s1_122.size()-1; j>=1; j--) {
    bkgd2_122+=bbggMass_s1_122[j]->Integral();
  }
  outfile << Form("Sig/Bkgd Ratio & %.3f \\\\ \\hline", (bbggMass_s1_122[0]->Integral())/bkgd2_122)<<endl;
  outfile << Form("Significance & %.3f \\\\ \\hline", (bbggMass_s1_122[0]->Integral())/sqrt(bkgd2_122))<<endl;
  
  legend->Draw();
  cv->SaveAs("bbggMass_s1_122_normalized.pdf");


  outfile << "\\multicolumn{2}{|c|}{\\textbf{Stage 1 Tight}}\\\\ \\hline" <<endl;
  outfile << "Sample Type & Expected Events \\\\ \\hline"<< endl;


  //*******************************************************************************************
  // M bbgg normalized Stage 1 Tight Optimistic Window
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.54,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = bbggMass_s1t_122.size()-1; i >= 0; i--) {
    if (bbggMass_s1t_122[i]->Integral() > 0) {
      outfile << Form("%s & %.4f \\\\ \\hline", LatexLabels[i].c_str(), bbggMass_s1t_122[i]->Integral())<<endl;

      TH1F* normHist = NormalizeHist(bbggMass_s1t_122[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  Float_t bkgd2t_122 = 0;
  for (Int_t j = bbggMass_s1t_122.size()-1; j>=1; j--) {
    bkgd2t_122+=bbggMass_s1t_122[j]->Integral();
  }
  outfile << Form("Sig/Bkgd Ratio & %.3f \\\\ \\hline", (bbggMass_s1t_122[0]->Integral())/bkgd2t_122)<<endl;
  outfile << Form("Significance & %.3f \\\\ \\hline", (bbggMass_s1t_122[0]->Integral())/sqrt(bkgd2t_122))<<endl;
  
  legend->Draw();
  cv->SaveAs("bbggMass_s1t_122_normalized.pdf");

  outfile << "\\end{tabular}"<<endl;
  outfile << "\\end{center}" <<endl;

  outfile << "\\begin{center}"<<endl;
  outfile << "\\begin{tabular}{| C{3.0cm}|| C{3.0cm}|}" <<endl;
  outfile << "\\hline" << endl;
  outfile << "\\multicolumn{2}{|c|}{\\textbf{Stage 2}}\\\\ \\hline" <<endl;
  outfile << "Sample Type & Expected Events \\\\ \\hline"<< endl;


  //*******************************************************************************************
  // M bbgg normalized Stage 2 Optimistic Window
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.54,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = bbggMass_s2_122.size()-1; i >= 0; i--) {
    if (bbggMass_s2_122[i]->Integral() > 0) {
      outfile << Form("%s & %.4f \\\\ \\hline", LatexLabels[i].c_str(), bbggMass_s2_122[i]->Integral())<<endl;

      TH1F* normHist = NormalizeHist(bbggMass_s2_122[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  Float_t bkgd3_122 = 0;
  for (Int_t j = bbggMass_s2_122.size()-1; j>=1; j--) {
    bkgd3_122+=bbggMass_s2_122[j]->Integral();
  }
  outfile << Form("Sig/Bkgd Ratio & %.3f \\\\ \\hline", (bbggMass_s2_122[0]->Integral())/bkgd3_122)<<endl;
  outfile << Form("Significance & %.3f \\\\ \\hline", (bbggMass_s2_122[0]->Integral())/sqrt(bkgd3_122))<<endl;
  
  legend->Draw();
  cv->SaveAs("bbggMass_s2_122_normalized.pdf");
  outfile << "\\multicolumn{2}{|c|}{\\textbf{Stage 2 Tight}}\\\\ \\hline" <<endl;
  outfile << "Sample Type & Expected Events \\\\ \\hline"<< endl;

  //*******************************************************************************************
  // M bbgg normalized Stage 2 Tight Optimistic Window
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.65,0.54,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = bbggMass_s2t_122.size()-1; i >= 0; i--) {
    if (bbggMass_s2t_122[i]->Integral() > 0) {
      outfile << Form("%s & %.4f \\\\ \\hline", LatexLabels[i].c_str(), bbggMass_s2t_122[i]->Integral())<<endl;

      TH1F* normHist = NormalizeHist(bbggMass_s2t_122[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        normHist->SetMaximum(normHist->GetMaximum()*1.5);
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  Float_t bkgd4_122 = 0;
  for (Int_t j = bbggMass_s2t_122.size()-1; j>=1; j--) {
    bkgd4_122+=bbggMass_s2t_122[j]->Integral();
  }
  outfile << Form("Sig/Bkgd Ratio & %.3f \\\\ \\hline", (bbggMass_s2t_122[0]->Integral())/bkgd4_122)<<endl;
  outfile << Form("Significance & %.3f \\\\ \\hline", (bbggMass_s2t_122[0]->Integral())/sqrt(bkgd4_122))<<endl;
  
  legend->Draw();
  cv->SaveAs("bbggMass_s2t_122_normalized.pdf");

  outfile << "\\end{tabular}"<<endl;
  outfile << "\\end{center}" <<endl;
  outfile << "\\end{multicols}"<<endl;




  outfile.close();
}

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

// Int_t bkgColors[10] = { kRed, kBlue,  kGreen+2, kMagenta, kCyan, kGreen+3, kBlack, kRed-6, kBlue-6, kMagenta+3 };
// string bkgLegendLabels[10] = { "HH#rightarrow bb#gamma#gamma", 
//                              "ttH,H#rightarrow#gamma#gamma",
//                              "ZH#rightarrow bb #gamma#gamma",
//                              "bbH (ggH)",
//                              "ttbar",
//                              "#gamma#gamma bb",
//                              "#gamma#gamma jj",
//                              "bbjj",
//                              "ccjj",
//                              "jjjj" };

const Int_t NBkgComponents = 7;
Int_t bkgColors[NBkgComponents] = { kRed, kGreen+2, kBlue, kMagenta+2, kCyan+2, kRed+2, kBlack };
string bkgLegendLabels[NBkgComponents] = { "HH#rightarrow bb#gamma#gamma", 
                                           "ZH",
                                           "ttH & bbH ",
                                           "#gamma#gamma bb",
                                           "#gamma#gamma jj",
                                           "4-jets",
                                           "t #bar{t}"};
Int_t bkgSampleTypes[NBkgComponents][3] = { {1,-1,-1} , 
                                            {3,-1,-1} , 
                                            {2, 4,-1} , 
                                            {6,-1,-1} , 
                                            {7,-1,-1} , 
                                            {8, 9,10} , 
                                            {5,-1,-1}   };

Int_t BkgComponentIndex( Int_t sampletype ) {

  for (uint i=0; i<NBkgComponents; ++i) {
    for (uint j=0; j<3; ++j) {
      if (sampletype == bkgSampleTypes[i][j]) {
        return i;        
      }
    }   
  }
  cout << "sample type " << sampletype << " does not correspond to any known bkg component. Quitting...\n";
  assert(false);
  return -1;  
}

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
void CutBasedAnalysis ( string InputFilename    = "/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/HHToBBGG/normalizedNtuples/HHToBBGGNtuple.main.root")
{

  double intLumi = 3000; //in units of fb^-1

 
  //*******************************************************************************************
  //Define Histograms
  //*******************************************************************************************
  vector<TH1F*> diphotonMass;
  vector<TH1F*> dibjetMass;
  vector<TH1F*> nlep;
  vector<TH1F*> njets;
  vector<TH1F*> ncentraljets;
  vector<TH1F*> dRgg;
  vector<TH1F*> dRbb;
  vector<TH1F*> minDRgb;

  vector<TH1F*> diphotonMassAfterScheme1Selection;
  vector<TH1F*> dibjetMassAfterScheme1Selection;
  vector<TH1F*> diphotonPt;
  vector<TH1F*> dibjetPt;
  vector<TH1F*> bbggMass;

  for (UInt_t i = 1; i <= NBkgComponents; ++i) {
    diphotonMass.push_back( new TH1F( Form("diphotonMass_%d",i), ";M_{#gamma#gamma} [GeV/c^{2}];Number of Events", 50, 100, 150));
    diphotonMass[i-1]->SetFillColor(bkgColors[i-1]);
    diphotonMass[i-1]->SetLineColor(bkgColors[i-1]);
    diphotonMass[i-1]->SetStats(false);

    dibjetMass.push_back( new TH1F( Form("dibjetMass_%d",i), ";M_{bb} [GeV/c^{2}];Number of Events", 26, 70, 200));
    dibjetMass[i-1]->SetFillColor(bkgColors[i-1]);
    dibjetMass[i-1]->SetLineColor(bkgColors[i-1]);
    dibjetMass[i-1]->SetStats(false);

    nlep.push_back( new TH1F( Form("nlep_%d",i), ";Number of Leptons;Number of Events", 6, -0.5, 5.5));
    nlep[i-1]->SetFillColor(bkgColors[i-1]);
    nlep[i-1]->SetLineColor(bkgColors[i-1]);
    nlep[i-1]->SetStats(false);

    njets.push_back( new TH1F( Form("njets_%d",i), ";Number of Jets;Number of Events", 15, -0.5, 14.5));
    njets[i-1]->SetFillColor(bkgColors[i-1]);
    njets[i-1]->SetLineColor(bkgColors[i-1]);
    njets[i-1]->SetStats(false);

    ncentraljets.push_back( new TH1F( Form("ncentraljets_%d",i), ";Number of Jets with |#eta|<2.5 ;Number of Events", 15, -0.5, 14.5));
    ncentraljets[i-1]->SetFillColor(bkgColors[i-1]);
    ncentraljets[i-1]->SetLineColor(bkgColors[i-1]);
    ncentraljets[i-1]->SetStats(false);

    dRgg.push_back( new TH1F( Form("dRgg_%d",i), ";#Delta R(#gamma,#gamma);Number of Events", 25, 0, 5));
    dRgg[i-1]->SetFillColor(bkgColors[i-1]);
    dRgg[i-1]->SetLineColor(bkgColors[i-1]);
    dRgg[i-1]->SetStats(false);

    dRbb.push_back( new TH1F( Form("dRgg_%d",i), ";#Delta R(b,b);Number of Events", 25, 0, 5));
    dRbb[i-1]->SetFillColor(bkgColors[i-1]);
    dRbb[i-1]->SetLineColor(bkgColors[i-1]);
    dRbb[i-1]->SetStats(false);

    minDRgb.push_back( new TH1F( Form("dRgg_%d",i), ";min #Delta R(#gamma,b);Number of Events", 25, 0, 5));
    minDRgb[i-1]->SetFillColor(bkgColors[i-1]);
    minDRgb[i-1]->SetLineColor(bkgColors[i-1]);
    minDRgb[i-1]->SetStats(false);


    diphotonPt.push_back( new TH1F( Form("diphotonPt_%d",i), "; p_{T #gamma#gamma};Number of Events", 25, 0, 300));
    diphotonPt[i-1]->SetFillColor(bkgColors[i-1]);
    diphotonPt[i-1]->SetLineColor(bkgColors[i-1]);
    diphotonPt[i-1]->SetStats(false);

    dibjetPt.push_back( new TH1F( Form("dibjetPt_%d",i), "; p_{T bb};Number of Events", 25, 0, 300));
    dibjetPt[i-1]->SetFillColor(bkgColors[i-1]);
    dibjetPt[i-1]->SetLineColor(bkgColors[i-1]);
    dibjetPt[i-1]->SetStats(false);


    diphotonMassAfterScheme1Selection.push_back( new TH1F( Form("diphotonMassAfterScheme1Selection_%d",i), ";M_{#gamma#gamma} [GeV/c^{2}];Number of Events", 25, 100, 150));
    diphotonMassAfterScheme1Selection[i-1]->SetFillColor(bkgColors[i-1]);
    diphotonMassAfterScheme1Selection[i-1]->SetLineColor(bkgColors[i-1]);
    diphotonMassAfterScheme1Selection[i-1]->SetStats(false);

    dibjetMassAfterScheme1Selection.push_back( new TH1F( Form("dibjetMassAfterScheme1Selection_%d",i), ";M_{bb} [GeV/c^{2}];Number of Events", 13, 70, 200));
    dibjetMassAfterScheme1Selection[i-1]->SetFillColor(bkgColors[i-1]);
    dibjetMassAfterScheme1Selection[i-1]->SetLineColor(bkgColors[i-1]);
    dibjetMassAfterScheme1Selection[i-1]->SetStats(false);

    bbggMass.push_back( new TH1F( Form("bbggMass_%d",i), "; m_{bb#gamma#gamma} [GeV/c^{2}];Number of Events", 50, 0, 1000));
    bbggMass[i-1]->SetFillColor(bkgColors[i-1]);
    bbggMass[i-1]->SetLineColor(bkgColors[i-1]);
    bbggMass[i-1]->SetStats(false);

 }

  //*******************************************************************************************
  //Define CutFlow Event Counts
  //*******************************************************************************************
  vector<double> NEventsPassStage1;
  vector<double> NEventsPassStage2;
  vector<double> NEventsPassStage3;
  for (UInt_t i = 1; i <= NBkgComponents; ++i) {
    NEventsPassStage1.push_back(0);
    NEventsPassStage2.push_back(0);
    NEventsPassStage3.push_back(0);
  }

  //*******************************************************************************************
  //Do analysis
  //*******************************************************************************************
                
  cmsana::HHToBBGGEventTree event;
  event.LoadTree(InputFilename.c_str());
  event.InitTree();

  cout << "Total Events: " << event.tree_->GetEntries() << "\n";
  for (int n=0;n<event.tree_->GetEntries();n++) { 
    
    event.tree_->GetEntry(n);
    if (n % 100000 == 0) cout << "Processing Event " << n << "\n";
        
    double weight = event.weight*intLumi;

    //make sure the sample type makes sense
    assert(event.sampletype >= 1 && event.sampletype <= 10);
    
     if (event.bjet1.Pt() > 30.0 && event.bjet2.Pt() > 30.0
         && event.pho1.Pt() > 25.0 && event.pho2.Pt() > 25.0
         && fmax(event.pho1.Pt(),event.pho2.Pt()) > 40.0
         && fabs(event.pho1.Eta()) < 2.5 && fabs(event.pho2.Eta()) < 2.5
         && fabs(event.bjet1.Eta()) < 2.4 && fabs(event.bjet2.Eta()) < 2.4
       ) {

      if (event.diphoton.M() > 100 && event.diphoton.M() < 150
           && event.dibjet.M() > 70 && event.dibjet.M() < 200 ) {

        diphotonMass[BkgComponentIndex(event.sampletype)]->Fill( event.diphoton.M() , weight);
        dibjetMass[BkgComponentIndex(event.sampletype)]->Fill( event.dibjet.M() , weight);
        nlep[BkgComponentIndex(event.sampletype)]->Fill( event.nlep , weight);
        njets[BkgComponentIndex(event.sampletype)]->Fill( event.njets , weight);
        ncentraljets[BkgComponentIndex(event.sampletype)]->Fill( event.ncentraljets , weight);

        dRgg[BkgComponentIndex(event.sampletype)]->Fill( event.DRgg , weight);
        dRbb[BkgComponentIndex(event.sampletype)]->Fill( event.DRbb , weight);
        minDRgb[BkgComponentIndex(event.sampletype)]->Fill( event.minDRgb , weight);        
        diphotonPt[BkgComponentIndex(event.sampletype)]->Fill( event.diphoton.Pt(), weight);
        dibjetPt[BkgComponentIndex(event.sampletype)]->Fill( event.dibjet.Pt(), weight);

        //*********************************************
        //mass distributions after Scheme 1 Selection
        //*********************************************
        if (event.nlep <= 0 && event.ncentraljets < 4
            && event.DRgg < 2.0 && event.DRbb < 2.0 && event.minDRgb > 1.5) {
          diphotonMassAfterScheme1Selection[BkgComponentIndex(event.sampletype)]->Fill( event.diphoton.M() , weight);
          if (event.diphoton.M() > 120 && event.diphoton.M() < 130 ) {
            dibjetMassAfterScheme1Selection[BkgComponentIndex(event.sampletype)]->Fill( event.dibjet.M() , weight);
          }
        }
        
        //*******************************************************************
        //distributions after Scheme 1 Selection and Mass Window Cuts
        //*******************************************************************
        if (event.nlep <= 0 && event.ncentraljets < 4
            && event.DRgg < 2.0 && event.DRbb < 2.0 && event.minDRgb > 1.5
            && event.diphoton.M() > 120 && event.diphoton.M() < 130
            && event.dibjet.M() > 105 && event.dibjet.M() < 145
          ) {         
          bbggMass[BkgComponentIndex(event.sampletype)]->Fill( event.bbgg.M(), weight);           
        }

      } //fit mass window requirements






      //*********************************************
      //Event Counts for Cut Based Analysis
      //*********************************************
      if (event.diphoton.M() > 100 && event.diphoton.M() < 150
          && event.dibjet.M() > 70 && event.dibjet.M() < 200 ) {
        NEventsPassStage1[BkgComponentIndex(event.sampletype)] += weight;
        
        if (event.nlep <= 0 && event.ncentraljets < 4 && 
            event.DRgg < 2.0 && event.DRbb < 2.0 && event.minDRgb > 1.5) {
          NEventsPassStage2[BkgComponentIndex(event.sampletype)] += weight;
          if (event.diphoton.M() > 120 && event.diphoton.M() < 130
              && event.dibjet.M() > 105 && event.dibjet.M() < 145) {
            NEventsPassStage3[BkgComponentIndex(event.sampletype)] += weight;
          }
        }
      }
      

     } //object selection requirements
  }
  
  
  //Draw Plots
  TCanvas *cv = 0;
  TLegend *legend = 0;
  bool firstdrawn = false;


  //*******************************************************************************************
  //Diphoton Mass
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.70,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackDiphotonMass = new THStack();
  cout << "after object selection\n";
  for (Int_t i = diphotonMass.size()-1; i >= 0; i--) {
    cout << "sample " << i << " : " << diphotonMass[i]->GetSumOfWeights() << "\n";
    if (diphotonMass[i]->Integral() > 0) {
      stackDiphotonMass->Add(diphotonMass[i]);
    }
  }
  for (Int_t i = 0 ; i < diphotonMass.size(); ++i) {
    legend->AddEntry(diphotonMass[i],bkgLegendLabels[i].c_str(), "F");
  }


  stackDiphotonMass->Draw();
  stackDiphotonMass->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackDiphotonMass->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackDiphotonMass->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackDiphotonMass->GetHists()->At(0)))->GetYaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs("MassGG.gif");
  cv->SaveAs("MassGG.pdf");


  //*******************************************************************************************
  //Dibjet Mass
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.70,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackDibjetMass = new THStack();
  for (Int_t i = dibjetMass.size()-1; i >= 0; i--) {
    if (dibjetMass[i]->Integral() > 0) {
      stackDibjetMass->Add(dibjetMass[i]);
    }
  }
  for (Int_t i = 0 ; i < dibjetMass.size(); ++i) {
    legend->AddEntry(dibjetMass[i],bkgLegendLabels[i].c_str(), "F");
  }

  stackDibjetMass->Draw();
  stackDibjetMass->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackDibjetMass->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackDibjetMass->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackDibjetMass->GetHists()->At(0)))->GetYaxis()->GetTitle());

  legend->Draw();
  cv->SaveAs("MassBB.gif");
  cv->SaveAs("MassBB.pdf");

  //*******************************************************************************************
  //Dibjet mass normalized
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.55,0.54,0.70,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = dibjetMass.size()-1; i >= 0; i--) {
    if (dibjetMass[i]->Integral() > 0) {
      
      TH1F* normHist = NormalizeHist(dibjetMass[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(2);
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
  cv->SaveAs("MassBB_normalized.gif");


  //*******************************************************************************************
  //NLeptons
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.55,0.54,0.70,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = nlep.size()-1; i >= 0; i--) {
    if (nlep[i]->Integral() > 0) {
      
      TH1F* normHist = NormalizeHist(nlep[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      normHist->GetYaxis()->SetTitle("Fraction of Events");
      if (!firstdrawn) {
        normHist->Draw("hist");
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  legend->Draw();
  cv->SaveAs("NLeptons.gif");
  cv->SaveAs("NLeptons.pdf");

  //*******************************************************************************************
  //NJets
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.55,0.54,0.70,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = njets.size()-1; i >= 0; i--) {
    if (njets[i]->Integral() > 0) {
      
      TH1F* normHist = NormalizeHist(njets[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(2);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      if (!firstdrawn) {
        normHist->Draw("hist");
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  legend->Draw();
  cv->SaveAs("Njets.gif");


  //*******************************************************************************************
  //Ncentraljets
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.55,0.54,0.70,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = ncentraljets.size()-1; i >= 0; i--) {
    if (ncentraljets[i]->Integral() > 0) {
      
      TH1F* normHist = NormalizeHist(ncentraljets[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(3);
      legend->AddEntry(normHist,bkgLegendLabels[i].c_str(), "L");
      normHist->GetYaxis()->SetTitle("Fraction of Events");
      if (!firstdrawn) {
        normHist->Draw("hist");
        firstdrawn = true;
      } else {
        normHist->Draw("hist,same");
      }

    }
  }

  legend->Draw();
  cv->SaveAs("Ncentraljets.gif");
  cv->SaveAs("Ncentraljets.pdf");


  //*******************************************************************************************
  //dRgg
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.15,0.54,0.55,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackDRgg = new THStack();
  cout << "After Mass windows + jet counting\n";
  for (Int_t i = dRgg.size()-1; i >= 0; i--) {
    cout << "sample " << i << " : " << dRgg[i]->GetSumOfWeights() << "\n";
    if (dRgg[i]->Integral() > 0) {
      stackDRgg->Add(dRgg[i]);
    }
  }
  for (Int_t i = 0 ; i < dRgg.size(); ++i) {
    legend->AddEntry(dRgg[i],bkgLegendLabels[i].c_str(), "F");
  }

  stackDRgg->Draw();
  stackDRgg->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackDRgg->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackDRgg->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackDRgg->GetHists()->At(0)))->GetYaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs("dRgg.gif");
  cv->SaveAs("dRgg.pdf");

  //*******************************************************************************************
  //dRbb
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.15,0.54,0.55,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackDRbb = new THStack();
  for (Int_t i = dRbb.size()-1; i >= 0; i--) {
    if (dRbb[i]->Integral() > 0) {
      stackDRbb->Add(dRbb[i]);
    }
  }
  for (Int_t i = 0 ; i < dRbb.size(); ++i) {
    legend->AddEntry(dRbb[i],bkgLegendLabels[i].c_str(), "F");
  }

  stackDRbb->Draw();
  stackDRbb->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackDRbb->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackDRbb->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackDRbb->GetHists()->At(0)))->GetYaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs("dRbb.gif");
  cv->SaveAs("dRbb.pdf");

  //*******************************************************************************************
  //minDRgb
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackminDRgb = new THStack();
  for (Int_t i = dRgg.size()-1; i >= 0; i--) {
    if (minDRgb[i]->Integral() > 0) {
      stackminDRgb->Add(minDRgb[i]);
    }
  }
  for (Int_t i = 0 ; i < minDRgb.size(); ++i) {
    legend->AddEntry(minDRgb[i],bkgLegendLabels[i].c_str(), "F");
  }
  stackminDRgb->Draw();
  stackminDRgb->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackminDRgb->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackminDRgb->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackminDRgb->GetHists()->At(0)))->GetYaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs("minDRgb.gif");
  cv->SaveAs("minDRgb.pdf");


  //*******************************************************************************************
  //dRgg normalized
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.55,0.54,0.70,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = ncentraljets.size()-1; i >= 0; i--) {
    if (ncentraljets[i]->Integral() > 0) {
      
      TH1F* normHist = NormalizeHist(dRgg[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(2);
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
  cv->SaveAs("dRgg_normalized.gif");

  //*******************************************************************************************
  //dRbb normalized
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.55,0.54,0.70,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = ncentraljets.size()-1; i >= 0; i--) {
    if (ncentraljets[i]->Integral() > 0) {
      
      TH1F* normHist = NormalizeHist(dRbb[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(2);
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
  cv->SaveAs("dRbb_normalized.gif");

  //*******************************************************************************************
  //minDRgb normalized
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.55,0.54,0.70,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = ncentraljets.size()-1; i >= 0; i--) {
    if (ncentraljets[i]->Integral() > 0) {
      
      TH1F* normHist = NormalizeHist(minDRgb[i]);
      normHist->SetFillStyle(0);
      normHist->SetLineWidth(2);
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
  cv->SaveAs("minDRgb_normalized.gif");

  //*******************************************************************************************
  //diphoton pt
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackdiphotonPt = new THStack();
  cout << "after DR cuts\n";
  for (Int_t i = diphotonPt.size()-1; i >= 0; i--) {
    cout << "sample " << i << " : " << diphotonPt[i]->GetSumOfWeights() << "\n";
    if (diphotonPt[i]->Integral() > 0) {
      stackdiphotonPt->Add(diphotonPt[i]);
    }
  }
  for (Int_t i = 0 ; i < diphotonPt.size(); ++i) {
    legend->AddEntry(diphotonPt[i],bkgLegendLabels[i].c_str(), "F");
  }
  stackdiphotonPt->Draw();
  stackdiphotonPt->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackdiphotonPt->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackdiphotonPt->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackdiphotonPt->GetHists()->At(0)))->GetYaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs("diphotonPt.gif");

  //*******************************************************************************************
  //dibjetPt
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackdibjetPt = new THStack();
  for (Int_t i = dibjetPt.size()-1; i >= 0; i--) {
    if (dibjetPt[i]->Integral() > 0) {
      stackdibjetPt->Add(dibjetPt[i]);
    }
  }
  for (Int_t i = 0 ; i < dibjetPt.size(); ++i) {
    legend->AddEntry(dibjetPt[i],bkgLegendLabels[i].c_str(), "F");
  }
  stackdibjetPt->Draw();
  stackdibjetPt->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackdibjetPt->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackdibjetPt->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackdibjetPt->GetHists()->At(0)))->GetYaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs("dibjetPt.gif");




  //*******************************************************************************************
  //Diphoton Mass After Scheme 1 Selection
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.70,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackDiphotonMassAfterScheme1Selection = new THStack();
  cout << "after object selection\n";
  for (Int_t i = diphotonMassAfterScheme1Selection.size()-1; i >= 0; i--) {
    cout << "sample " << i << " : " << diphotonMassAfterScheme1Selection[i]->GetSumOfWeights() << "\n";
    if (diphotonMassAfterScheme1Selection[i]->Integral() > 0) {
      stackDiphotonMassAfterScheme1Selection->Add(diphotonMassAfterScheme1Selection[i]);
      legend->AddEntry(diphotonMassAfterScheme1Selection[i],bkgLegendLabels[i].c_str(), "F");
    }
  }
  stackDiphotonMassAfterScheme1Selection->Draw();
  stackDiphotonMassAfterScheme1Selection->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackDiphotonMassAfterScheme1Selection->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackDiphotonMassAfterScheme1Selection->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackDiphotonMassAfterScheme1Selection->GetHists()->At(0)))->GetYaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs("MassGGAfterScheme1Selection.gif");
  cv->SaveAs("MassGGAfterScheme1Selection.pdf");


  //*******************************************************************************************
  //Dibjet Mass After Scheme 1 Selection
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.70,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackDibjetMassAfterScheme1Selection = new THStack();
  for (Int_t i = dibjetMassAfterScheme1Selection.size()-1; i >= 0; i--) {
    if (dibjetMassAfterScheme1Selection[i]->Integral() > 0) {
      stackDibjetMassAfterScheme1Selection->Add(dibjetMassAfterScheme1Selection[i]);
      legend->AddEntry(dibjetMassAfterScheme1Selection[i],bkgLegendLabels[i].c_str(), "F");
    }
  }

  stackDibjetMassAfterScheme1Selection->Draw();
  stackDibjetMassAfterScheme1Selection->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackDibjetMassAfterScheme1Selection->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackDibjetMassAfterScheme1Selection->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackDibjetMassAfterScheme1Selection->GetHists()->At(0)))->GetYaxis()->GetTitle());

  legend->Draw();
  cv->SaveAs("MassBBAfterScheme1Selection.gif");
  cv->SaveAs("MassBBAfterScheme1Selection.pdf");




  //*******************************************************************************************
  //bbgg mass
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackbbggMass = new THStack();
  for (Int_t i = bbggMass.size()-1; i >= 0; i--) {
    if (bbggMass[i]->Integral() > 0) {
      stackbbggMass->Add(bbggMass[i]);
      legend->AddEntry(bbggMass[i],bkgLegendLabels[i].c_str(), "F");
    }
  }
  stackbbggMass->Draw();
  stackbbggMass->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackbbggMass->GetHists()->At(0)))->GetXaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs("bbggMass_AfterSelStep2.gif");




  //*******************************************************************************************
  //Show event count summary
  //*******************************************************************************************
  cout << "\n";
  cout << "Event Counts\n";
  for (UInt_t i = 0; i < NBkgComponents; ++i) {
    cout << bkgLegendLabels[i] << " : "
         << NEventsPassStage1[i] << " " 
         << NEventsPassStage2[i] << " " 
         << NEventsPassStage3[i] << " " 
         << "\n";
  }

}

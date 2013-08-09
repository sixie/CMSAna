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

Int_t bkgColors[10] = { kRed, kBlue,  kGreen+2, kMagenta, kCyan, kBlue-9, kBlack, kRed-6, kBlue-6, kMagenta+3 };
string bkgLegendLabels[10] = { "HH#rightarrow bb#gamma#gamma", 
                             "ttH,H#rightarrow#gamma#gamma",
                             "ZH#rightarrow bb #gamma#gamma",
                             "bbH (ggH)",
                             "ttbar",
                             "#gamma#gamma bb",
                             "#gamma#gamma jj",
                             "bbjj",
                             "ccjj",
                             "jjjj" };

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

  vector<TH1F*> diphotonPt;
  vector<TH1F*> dibjetPt;
  vector<TH1F*> bbggMass;

  for (UInt_t i = 1; i <= 10; ++i) {
    diphotonMass.push_back( new TH1F( Form("diphotonMass_%d",i), ";M_{#gamma#gamma} [GeV/c^{2}];Number of Events", 100, 100, 150));
    diphotonMass[i-1]->SetFillColor(bkgColors[i-1]);
    diphotonMass[i-1]->SetLineColor(bkgColors[i-1]);
    diphotonMass[i-1]->SetStats(false);

    dibjetMass.push_back( new TH1F( Form("dibjetMass_%d",i), ";M_{bb} [GeV/c^{2}];Number of Events", 100, 0, 300));
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

    ncentraljets.push_back( new TH1F( Form("ncentraljets_%d",i), ";Number of Central Jets;Number of Events", 15, -0.5, 14.5));
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

    diphotonPt.push_back( new TH1F( Form("diphotonPt_%d",i), "; p_{T #gamma#gamma};Number of Events", 25, 50, 350));
    diphotonPt[i-1]->SetFillColor(bkgColors[i-1]);
    diphotonPt[i-1]->SetLineColor(bkgColors[i-1]);
    diphotonPt[i-1]->SetStats(false);

    dibjetPt.push_back( new TH1F( Form("dibjetPt_%d",i), "; p_{T bb};Number of Events", 25, 0, 400));
    dibjetPt[i-1]->SetFillColor(bkgColors[i-1]);
    dibjetPt[i-1]->SetLineColor(bkgColors[i-1]);
    dibjetPt[i-1]->SetStats(false);

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
  for (UInt_t i = 1; i <= 10; ++i) {
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

  for (int n=0;n<event.tree_->GetEntries();n++) { 
    event.tree_->GetEntry(n);

    double weight = event.weight*intLumi;

    //make sure the sample type makes sense
     assert(event.sampletype >= 1 && event.sampletype <= 10);
    
//      if (event.sampletype != 1) continue;

     
     if (event.bjet1.Pt() > 30.0 && event.bjet2.Pt() > 30.0
         && event.pho1.Pt() > 25.0 && event.pho2.Pt() > 25.0
         && fmax(event.pho1.Pt(),event.pho2.Pt()) > 40.0
         && fabs(event.pho1.Eta()) < 2.5 && fabs(event.pho2.Eta()) < 2.5
         && fabs(event.bjet1.Eta()) < 2.4 && fabs(event.bjet2.Eta()) < 2.4
       ) {

       diphotonMass[event.sampletype-1]->Fill( event.diphoton.M() , weight);
       dibjetMass[event.sampletype-1]->Fill( event.dibjet.M() , weight);
       nlep[event.sampletype-1]->Fill( event.nlep , weight);
       njets[event.sampletype-1]->Fill( event.njets , weight);
       ncentraljets[event.sampletype-1]->Fill( event.ncentraljets , weight);

       if (event.diphoton.M() > 100 && event.diphoton.M() < 150
           && event.dibjet.M() > 60 && event.dibjet.M() < 200 ) {
         NEventsPassStage1[event.sampletype-1] += weight;
       }

       if (event.nlep <= 0 && event.ncentraljets <= 4
           && event.diphoton.M() > 120 && event.diphoton.M() < 130
           && event.dibjet.M() > 105 && event.dibjet.M() < 145
         ) {
         
         dRgg[event.sampletype-1]->Fill( event.DRgg , weight);
         dRbb[event.sampletype-1]->Fill( event.DRbb , weight);
         minDRgb[event.sampletype-1]->Fill( event.minDRgb , weight);
         NEventsPassStage2[event.sampletype-1] += weight;

         if (event.DRgg < 2.0 && event.minDRgb > 1.0) {
           diphotonPt[event.sampletype-1]->Fill( event.diphoton.Pt(), weight);
           dibjetPt[event.sampletype-1]->Fill( event.dibjet.Pt(), weight);
           bbggMass[event.sampletype-1]->Fill( event.bbgg.M(), weight);
           NEventsPassStage3[event.sampletype-1] += weight;
         }

       }

     }
  }
  
  
  //Draw Plots
  TCanvas *cv = 0;
  TLegend *legend = 0;
  bool firstdrawn = false;


  //*******************************************************************************************
  //Diphoton Mass
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackDiphotonMass = new THStack();
  cout << "after object selection\n";
  for (Int_t i = diphotonMass.size()-1; i >= 0; i--) {
    cout << "sample " << i << " : " << diphotonMass[i]->GetSumOfWeights() << "\n";
    if (diphotonMass[i]->Integral() > 0) {
      stackDiphotonMass->Add(diphotonMass[i]);
      legend->AddEntry(diphotonMass[i],bkgLegendLabels[i].c_str(), "F");
    }
  }
  stackDiphotonMass->Draw();
  stackDiphotonMass->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackDiphotonMass->GetHists()->At(0)))->GetXaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs("MassGG.gif");


  //*******************************************************************************************
  //Dibjet Mass
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.15,0.54,0.40,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackDibjetMass = new THStack();
  for (Int_t i = dibjetMass.size()-1; i >= 0; i--) {
    if (dibjetMass[i]->Integral() > 0) {
      stackDibjetMass->Add(dibjetMass[i]);
      legend->AddEntry(dibjetMass[i],bkgLegendLabels[i].c_str(), "F");
    }
  }

  stackDibjetMass->Draw();
  stackDibjetMass->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackDibjetMass->GetHists()->At(0)))->GetXaxis()->GetTitle());

  legend->Draw();
  cv->SaveAs("MassBB.gif");

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
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = nlep.size()-1; i >= 0; i--) {
    if (nlep[i]->Integral() > 0) {
      
      TH1F* normHist = NormalizeHist(nlep[i]);
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
  cv->SaveAs("NLeptons.gif");

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
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  
  firstdrawn = false;
  for (Int_t i = ncentraljets.size()-1; i >= 0; i--) {
    if (ncentraljets[i]->Integral() > 0) {
      
      TH1F* normHist = NormalizeHist(ncentraljets[i]);
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
  cv->SaveAs("Ncentraljets.gif");


  //*******************************************************************************************
  //dRgg
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackDRgg = new THStack();
  cout << "After Mass windows + jet counting\n";
  for (Int_t i = dRgg.size()-1; i >= 0; i--) {
    cout << "sample " << i << " : " << dRgg[i]->GetSumOfWeights() << "\n";
    if (dRgg[i]->Integral() > 0) {
      stackDRgg->Add(dRgg[i]);
      legend->AddEntry(dRgg[i],bkgLegendLabels[i].c_str(), "F");
    }
  }
  stackDRgg->Draw();
  stackDRgg->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackDRgg->GetHists()->At(0)))->GetXaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs("dRgg_AfterSelStep1.gif");

  //*******************************************************************************************
  //dRbb
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.54,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackDRbb = new THStack();
  for (Int_t i = dRbb.size()-1; i >= 0; i--) {
    if (dRbb[i]->Integral() > 0) {
      stackDRbb->Add(dRbb[i]);
      legend->AddEntry(dRbb[i],bkgLegendLabels[i].c_str(), "F");
    }
  }
  stackDRbb->Draw();
  stackDRbb->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackDRbb->GetHists()->At(0)))->GetXaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs("dRbb_AfterSelStep1.gif");

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
      legend->AddEntry(minDRgb[i],bkgLegendLabels[i].c_str(), "F");
    }
  }
  stackminDRgb->Draw();
  stackminDRgb->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackminDRgb->GetHists()->At(0)))->GetXaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs("minDRgb_AfterSelStep1.gif");


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
      legend->AddEntry(diphotonPt[i],bkgLegendLabels[i].c_str(), "F");
    }
  }
  stackdiphotonPt->Draw();
  stackdiphotonPt->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackdiphotonPt->GetHists()->At(0)))->GetXaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs("diphotonPt_AfterSelStep2.gif");

  //*******************************************************************************************
  //diphoton pt
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
      legend->AddEntry(dibjetPt[i],bkgLegendLabels[i].c_str(), "F");
    }
  }
  stackdibjetPt->Draw();
  stackdibjetPt->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackdibjetPt->GetHists()->At(0)))->GetXaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs("dibjetPt_AfterSelStep2.gif");

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
  for (UInt_t i = 0; i < 10; ++i) {
    cout << bkgLegendLabels[i] << " : "
         << NEventsPassStage1[i] << " " 
         << NEventsPassStage2[i] << " " 
         << NEventsPassStage3[i] << " " 
         << "\n";
  }

}

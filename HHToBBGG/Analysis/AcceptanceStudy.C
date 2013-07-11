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

TH1F* MakeCutEfficiencyPlot(TH1F *hist, bool cutBelow = true) {

  TH1F *eff = (TH1F*)hist->Clone((string(hist->GetName())+"_cuteff").c_str());
  Double_t integral = 0;
  eff->SetTitle("");
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    integral += hist->GetBinContent(b);
  }

  for (UInt_t b=0; int(b)<eff->GetXaxis()->GetNbins()+2; ++b) {
    double tmp = 0;
    if (cutBelow) {
      //integrate from bin b to the last bin
      for (UInt_t c=b; int(c)<hist->GetXaxis()->GetNbins()+2; ++c) {
        tmp += hist->GetBinContent(c);
      }
    } else {
      //integrate from bin 0 up to bin b
      for (UInt_t c=0; int(c)<=b; ++c) {
        tmp += hist->GetBinContent(c);
      }
    }

    //cout << "bin " << b << " : " << tmp << " " << integral << " \n";
    eff->SetBinContent(b, tmp / integral);
    //ignore uncertainties for now
  }

  eff->GetYaxis()->SetTitle("Cut Efficiency");
  return eff;
}


//------------------------------------------------------------------------------
// PlotHiggsRes_LP
//------------------------------------------------------------------------------
void AcceptanceStudy ( string InputFilename    = "/afs/cern.ch/work/s/sixie/public/HHToBBGG/normalizedNtuples/HHToBBGGNtuple.HHtoBBGG-14tev-START53_V7A.normalized.root")
{

  const double intLumi = 3000;

  //*******************************************************************************************
  //Define Histograms
  //*******************************************************************************************
  double TotalEvents = 0;
  double EventsPassPhotonPt_25_25 = 0;
  double EventsPassPhotonPt_40_25 = 0;
  double EventsPassPhotonPt_40_20 = 0;
  double EventsPassPhotonPt_40_15 = 0;
  double EventsPassPhotonPt_50_20 = 0;
  double EventsPassPhotonPt_60_15 = 0;
  double EventsPassPhotonPt_40_25_BothEta1p5 = 0;
  double EventsPassPhotonPt_40_25_BothEta2p0 = 0;
  double EventsPassPhotonPt_40_25_BothEta2p5 = 0;
  double EventsPassPhotonPt_40_25_BothEta3p0 = 0;
  double EventsPassPhotonPt_40_25_BothEta3p5 = 0;
  double EventsPassPhotonPt_40_25_BothEta4p0 = 0;

  double EventsPassBJetPt_25_25 = 0;
  double EventsPassBJetPt_30_30 = 0;
  double EventsPassBJetPt_30_25 = 0;
  double EventsPassBJetPt_30_20 = 0;
  double EventsPassBJetPt_45_25 = 0;
  double EventsPassBJetPt_60_20 = 0;
  double EventsPassBJetPt_30_30_BothEta2p4 = 0;
  double EventsPassBJetPt_30_30_BothEta3p0 = 0;
  double EventsPassBJetPt_30_30_BothEta3p5 = 0;
  double EventsPassBJetPt_30_30_BothEta4p0 = 0;


  TH1F *hist_photon1Pt =  new TH1F( "photon1Pt", "; p_{T #gamma 1} [GeV/c];Number of Events", 100, 0, 400);
  TH1F *hist_photon2Pt =  new TH1F( "photon2Pt", "; p_{T #gamma 2} [GeV/c];Number of Events", 100, 0, 200);

  TH1F *hist_photonEta =  new TH1F( "photonEta", "; #eta_{#gamma } ;Number of Events", 100, 0, 5);
  TH1F *hist_photonMaxEta =  new TH1F( "photonMaxEta", "; max( #eta_{#gamma 1}, #eta_{#gamma 2}) ;Number of Events", 100, 0, 5);
  TH1F *hist_photonMaxEta_AfterPtCuts =  new TH1F( "photonMaxEta_AfterPtCuts", "; max( #eta_{#gamma 1}, #eta_{#gamma 2}) ;Number of Events", 100, 0, 5);

  TH1F *hist_bjet1Pt =  new TH1F( "bjet1Pt", "; p_{T b1} [GeV/c];Number of Events", 100, 0, 400);
  TH1F *hist_bjet2Pt =  new TH1F( "bjet2Pt", "; p_{T b2} [GeV/c];Number of Events", 100, 0, 200);

  TH1F *hist_bjetEta =  new TH1F( "bjetEta", "; #eta_{b} ;Number of Events", 100, 0, 5);
  TH1F *hist_bjetMaxEta =  new TH1F( "bjetMaxEta", "; max( #eta_{b1}, #eta_{b2}) ;Number of Events", 100, 0, 5);

  TH1F *hist_bjet2Pt_AfterBaselinePhotonCuts =  new TH1F( "bjet2Pt_AfterBaselinePhotonCuts", "; p_{T b2} [GeV/c];Number of Events", 100, 0, 200);
  TH1F *hist_bjetMaxEta_AfterBaselinePhotonCuts =  new TH1F( "bjetMaxEta_AfterBaselinePhotonCuts", "; max( #eta_{b1}, #eta_{b2}) ;Number of Events", 100, 0, 5);
  TH1F *hist_photon2Pt_AfterBaselineBJetCuts =  new TH1F( "photon2Pt_AfterBaselineBJetCuts", "; p_{T #gamma 2} [GeV/c];Number of Events", 100, 0, 200);
  TH1F *hist_photonMaxEta_AfterBaselineBJetCuts =  new TH1F( "photonMaxEta_AfterBaselineBJetCuts", "; max( #eta_{#gamma 1}, #eta_{#gamma 2}) ;Number of Events", 100, 0, 5);



  double EventsPass_PhotonPt4025Eta1p5_BJetPt3030Eta2p4 = 0;
  double EventsPass_PhotonPt4025Eta2p0_BJetPt3030Eta2p4 = 0;
  double EventsPass_PhotonPt4025Eta2p5_BJetPt3030Eta2p4 = 0;
  double EventsPass_PhotonPt4025Eta2p5_BJetPt3030Eta3p0 = 0;
  double EventsPass_PhotonPt4025Eta3p0_BJetPt3030Eta2p4 = 0;
  double EventsPass_PhotonPt4025Eta3p0_BJetPt3030Eta3p0 = 0;

  double EventsPass_PhotonPt6015Eta2p5_BJetPt6020Eta2p4 = 0;
  double EventsPass_PhotonPt6015Eta2p5_BJetPt6020Eta3p0 = 0;
  double EventsPass_PhotonPt6015Eta3p0_BJetPt6020Eta2p4 = 0;
  double EventsPass_PhotonPt6015Eta3p0_BJetPt6020Eta3p0 = 0;


  //*******************************************************************************************
  //Do analysis
  //*******************************************************************************************
                
  cmsana::HHToBBGGEventTree event;
  event.LoadTree(InputFilename.c_str());
  event.InitTree();

  for (int n=0;n<event.tree_->GetEntries();n++) { 
    event.tree_->GetEntry(n);
    double weight = event.weight * intLumi;
    TotalEvents += weight;

    hist_photon1Pt->Fill( fmin( fmax(event.genpho1.Pt(), event.genpho2.Pt()), 399.99) );
    hist_photon2Pt->Fill( fmin( fmin(event.genpho1.Pt(), event.genpho2.Pt()), 199.99) );
    hist_bjet1Pt->Fill( fmin( fmax(event.genbjet1.Pt(), event.genbjet2.Pt()), 399.99) );
    hist_bjet2Pt->Fill( fmin( fmin(event.genbjet1.Pt(), event.genbjet2.Pt()), 199.99) );
    
    double photonMaxEta = fmax( fabs(event.genpho1.Eta()),fabs(event.genpho2.Eta())); 
    hist_photonMaxEta->Fill( photonMaxEta );
    hist_photonEta->Fill( fabs(event.genpho1.Eta()) );
    hist_photonEta->Fill( fabs(event.genpho2.Eta()) );

    double bjetMaxEta = fmax( fabs(event.genbjet1.Eta()),fabs(event.genbjet2.Eta()));    
    if (event.genbjet1.Pt() > 30 && event.genbjet2.Pt() > 30) {
      hist_bjetMaxEta->Fill( bjetMaxEta );    
      hist_bjetEta->Fill( fabs(event.genbjet1.Eta()) );
      hist_bjetEta->Fill( fabs(event.genbjet2.Eta()) );
    }
   
      
    //*******************************************************************************************
    //Study Photon Acceptance
    //*******************************************************************************************
    if ( event.genpho1.Pt() > 25 && event.genpho2.Pt() > 25 ) EventsPassPhotonPt_25_25 += weight;    

    if ( fmax(event.genpho1.Pt(),event.genpho2.Pt()) > 40 && event.genpho1.Pt() > 25 && event.genpho2.Pt() > 25) {
      EventsPassPhotonPt_40_25 += weight;        
      hist_photonMaxEta_AfterPtCuts->Fill( photonMaxEta ); 

      if ( photonMaxEta < 1.5) EventsPassPhotonPt_40_25_BothEta1p5 += weight;
      if ( photonMaxEta < 2.0) EventsPassPhotonPt_40_25_BothEta2p0 += weight;
      if ( photonMaxEta < 2.5) EventsPassPhotonPt_40_25_BothEta2p5 += weight;
      if ( photonMaxEta < 3.0) EventsPassPhotonPt_40_25_BothEta3p0 += weight;
      if ( photonMaxEta < 3.5) EventsPassPhotonPt_40_25_BothEta3p5 += weight;
      if ( photonMaxEta < 4.0) EventsPassPhotonPt_40_25_BothEta4p0 += weight;

    }    

    //Lower photon2 pt cuts
    if ( fmax(event.genpho1.Pt(),event.genpho2.Pt()) > 40 && fmin(event.genpho1.Pt(),event.genpho2.Pt()) > 20 ) 
      EventsPassPhotonPt_40_20 += weight;
    if ( fmax(event.genpho1.Pt(),event.genpho2.Pt()) > 40 && fmin(event.genpho1.Pt(),event.genpho2.Pt()) > 15 ) 
      EventsPassPhotonPt_40_15 += weight;
    if ( fmax(event.genpho1.Pt(),event.genpho2.Pt()) > 50 && fmin(event.genpho1.Pt(),event.genpho2.Pt()) > 20 ) 
      EventsPassPhotonPt_50_20 += weight;
    if ( fmax(event.genpho1.Pt(),event.genpho2.Pt()) > 60 && fmin(event.genpho1.Pt(),event.genpho2.Pt()) > 15 ) 
      EventsPassPhotonPt_60_15 += weight;

    //AfterBaselineBJetCuts
    if (event.genbjet1.Pt() > 30 && event.genbjet2.Pt() > 30 
        && fabs(event.genbjet1.Eta()) < 2.4 && fabs(event.genbjet2.Eta()) < 2.4) {
      hist_photon2Pt_AfterBaselineBJetCuts->Fill( fmin( fmin(event.genpho1.Pt(), event.genpho2.Pt()), 199.99) );
      hist_photonMaxEta_AfterBaselineBJetCuts->Fill( photonMaxEta );
    }


    //*******************************************************************************************
    //Study BJet Acceptance
    //*******************************************************************************************
    if ( event.genbjet1.Pt() > 25 && event.genbjet2.Pt() > 25 ) EventsPassBJetPt_25_25 += weight;    
    if ( event.genbjet1.Pt() > 30 && event.genbjet2.Pt() > 30) {
      EventsPassBJetPt_30_30 += weight;

      if ( bjetMaxEta < 2.4) EventsPassBJetPt_30_30_BothEta2p4 += weight;
      if ( bjetMaxEta < 3.0) EventsPassBJetPt_30_30_BothEta3p0 += weight;
      if ( bjetMaxEta < 3.5) EventsPassBJetPt_30_30_BothEta3p5 += weight;
      if ( bjetMaxEta < 4.0) EventsPassBJetPt_30_30_BothEta4p0 += weight;

    }    

    //Lower photon2 pt cuts
    if ( fmax(event.genbjet1.Pt(),event.genbjet2.Pt()) > 30 && fmin(event.genbjet1.Pt(),event.genbjet2.Pt()) > 25 ) 
      EventsPassBJetPt_30_25 += weight;
    if ( fmax(event.genbjet1.Pt(),event.genbjet2.Pt()) > 30 && fmin(event.genbjet1.Pt(),event.genbjet2.Pt()) > 20 ) 
      EventsPassBJetPt_30_20 += weight;
    if ( fmax(event.genbjet1.Pt(),event.genbjet2.Pt()) > 45 && fmin(event.genbjet1.Pt(),event.genbjet2.Pt()) > 25 ) 
      EventsPassBJetPt_45_25 += weight;
    if ( fmax(event.genbjet1.Pt(),event.genbjet2.Pt()) > 60 && fmin(event.genbjet1.Pt(),event.genbjet2.Pt()) > 20 ) 
      EventsPassBJetPt_60_20 += weight;

    //AfterBaselinePhotonCuts
    if ( fmax(event.genpho1.Pt(),event.genpho2.Pt()) > 40 && event.genpho1.Pt() > 25 && event.genpho2.Pt() > 25
         && fabs(event.genpho1.Eta()) < 2.5 && fabs(event.genpho2.Eta()) < 2.5) {
      hist_bjet2Pt_AfterBaselinePhotonCuts->Fill( fmin( fmin(event.genbjet1.Pt(), event.genbjet2.Pt()), 199.99) );
      hist_bjetMaxEta_AfterBaselinePhotonCuts->Fill( bjetMaxEta );
    }


    //*******************************************************************************************
    //Study Full Acceptance
    //*******************************************************************************************
    if ( fmax( event.genpho1.Pt(), event.genpho2.Pt()) > 40 && 
         fmin( event.genpho1.Pt(), event.genpho2.Pt()) > 25 &&
         fmax( event.genbjet1.Pt(), event.genbjet2.Pt()) > 30 && 
         fmin( event.genbjet1.Pt(), event.genbjet2.Pt()) > 30 
      ) {
      //no Ecal Endcap scenario
      if ( fabs( event.genpho1.Eta() ) < 1.5 && fabs( event.genpho2.Eta() ) < 1.5
           && fabs( event.genbjet1.Eta() ) < 2.4 && fabs( event.genbjet2.Eta() ) < 2.4 ) 
        EventsPass_PhotonPt4025Eta1p5_BJetPt3030Eta2p4 += weight;
      if ( fabs( event.genpho1.Eta() ) < 2.0 && fabs( event.genpho2.Eta() ) < 2.0
           && fabs( event.genbjet1.Eta() ) < 2.4 && fabs( event.genbjet2.Eta() ) < 2.4 ) 
        EventsPass_PhotonPt4025Eta2p0_BJetPt3030Eta2p4 += weight;

      if ( fabs( event.genpho1.Eta() ) < 2.5 && fabs( event.genpho2.Eta() ) < 2.5
           && fabs( event.genbjet1.Eta() ) < 2.4 && fabs( event.genbjet2.Eta() ) < 2.4 ) 
        EventsPass_PhotonPt4025Eta2p5_BJetPt3030Eta2p4 += weight;
      if ( fabs( event.genpho1.Eta() ) < 2.5 && fabs( event.genpho2.Eta() ) < 2.5
           && fabs( event.genbjet1.Eta() ) < 3.0 && fabs( event.genbjet2.Eta() ) < 3.0 ) 
        EventsPass_PhotonPt4025Eta2p5_BJetPt3030Eta3p0 += weight;
      if ( fabs( event.genpho1.Eta() ) < 3.0 && fabs( event.genpho2.Eta() ) < 3.0
           && fabs( event.genbjet1.Eta() ) < 2.4 && fabs( event.genbjet2.Eta() ) < 2.4 )
        EventsPass_PhotonPt4025Eta3p0_BJetPt3030Eta2p4 += weight;
      if ( fabs( event.genpho1.Eta() ) < 3.0 && fabs( event.genpho2.Eta() ) < 3.0
           && fabs( event.genbjet1.Eta() ) < 3.0 && fabs( event.genbjet2.Eta() ) < 3.0 )
        EventsPass_PhotonPt4025Eta3p0_BJetPt3030Eta3p0 += weight;
    }
    
    if ( fmax( event.genpho1.Pt(), event.genpho2.Pt()) > 60 && 
         fmin( event.genpho1.Pt(), event.genpho2.Pt()) > 15 &&
         fmax( event.genbjet1.Pt(), event.genbjet2.Pt()) > 60 && 
         fmin( event.genbjet1.Pt(), event.genbjet2.Pt()) > 20
      ) {
      if ( fabs( event.genpho1.Eta() ) < 2.5 && fabs( event.genpho2.Eta() ) < 2.5
           && fabs( event.genbjet1.Eta() ) < 2.4 && fabs( event.genbjet2.Eta() ) < 2.4 ) 
        EventsPass_PhotonPt6015Eta2p5_BJetPt6020Eta2p4 += weight;
      if ( fabs( event.genpho1.Eta() ) < 2.5 && fabs( event.genpho2.Eta() ) < 2.5
           && fabs( event.genbjet1.Eta() ) < 3.0 && fabs( event.genbjet2.Eta() ) < 3.0 ) 
        EventsPass_PhotonPt6015Eta2p5_BJetPt6020Eta3p0 += weight;
      if ( fabs( event.genpho1.Eta() ) < 3.0 && fabs( event.genpho2.Eta() ) < 3.0
           && fabs( event.genbjet1.Eta() ) < 2.4 && fabs( event.genbjet2.Eta() ) < 2.4 )
        EventsPass_PhotonPt6015Eta3p0_BJetPt6020Eta2p4 += weight;
      if ( fabs( event.genpho1.Eta() ) < 3.0 && fabs( event.genpho2.Eta() ) < 3.0
           && fabs( event.genbjet1.Eta() ) < 3.0 && fabs( event.genbjet2.Eta() ) < 3.0 )
        EventsPass_PhotonPt6015Eta3p0_BJetPt6020Eta3p0 += weight;
    }
    
   


  }
  
  
  //Draw Plots
  TCanvas *cv = 0;
  TLegend *legend = 0;

  //*******************************************************************************************
  //Normalize
  //*******************************************************************************************
  NormalizeHist(hist_photon1Pt);
  NormalizeHist(hist_photon2Pt);
  NormalizeHist(hist_photonEta);
  NormalizeHist(hist_photonMaxEta);
  NormalizeHist(hist_photonMaxEta_AfterPtCuts);

  NormalizeHist(hist_bjet1Pt);
  NormalizeHist(hist_bjet2Pt);
  NormalizeHist(hist_bjetEta);
  NormalizeHist(hist_bjetMaxEta);

  NormalizeHist(hist_photon2Pt_AfterBaselineBJetCuts);
  NormalizeHist(hist_photonMaxEta_AfterBaselineBJetCuts);
  NormalizeHist(hist_bjet2Pt_AfterBaselinePhotonCuts);
  NormalizeHist(hist_bjetMaxEta_AfterBaselinePhotonCuts);

  //*******************************************************************************************
  //Photon Pt
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  hist_photon1Pt->Draw("hist");
  cv->SaveAs("GenPhoton1Pt.gif");

  cv = new TCanvas("cv","cv", 800,600);
  TH1F *hist_photon1PtCutEff = MakeCutEfficiencyPlot(hist_photon1Pt, true);
  hist_photon1PtCutEff->Draw("hist");
  cv->SaveAs("GenPhoton1PtCutEff.gif");
  

  cv = new TCanvas("cv","cv", 800,600);
  hist_photon2Pt->Draw("hist");
  cv->SaveAs("GenPhoton2Pt.gif");

  cv = new TCanvas("cv","cv", 800,600);
  TH1F *hist_photon2PtCutEff = MakeCutEfficiencyPlot(hist_photon2Pt, true);
  hist_photon2PtCutEff->Draw("hist");
  cv->SaveAs("GenPhoton2PtCutEff.gif");

  //*******************************************************************************************
  //Photon Eta
  //*******************************************************************************************

  cv = new TCanvas("cv","cv", 800,600);
  hist_photonEta->Draw("hist");
  cv->SaveAs("GenPhotonEta.gif");

  cv = new TCanvas("cv","cv", 800,600);
  hist_photonMaxEta->Draw("hist");
  cv->SaveAs("GenPhotonMaxEta.gif");

  cv = new TCanvas("cv","cv", 800,600);
  hist_photonMaxEta_AfterPtCuts->Draw("hist");
  cv->SaveAs("GenPhotonMaxEta_AfterPtCuts.gif");

  cv = new TCanvas("cv","cv", 800,600);
  TH1F *hist_photonMaxEtaCutEff = MakeCutEfficiencyPlot(hist_photonMaxEta, false);
  TH1F *hist_photonMaxEtaCutEff_AfterPtCuts = MakeCutEfficiencyPlot(hist_photonMaxEta_AfterPtCuts, false);
  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(hist_photonMaxEtaCutEff, "All Events", "L");
  legend->AddEntry(hist_photonMaxEtaCutEff_AfterPtCuts, "Events Passing 40/25 p_{T} Cut", "L");

  hist_photonMaxEtaCutEff_AfterPtCuts->SetLineColor(kBlue);
  hist_photonMaxEtaCutEff->Draw("hist");
  hist_photonMaxEtaCutEff_AfterPtCuts->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("GenPhotonMaxEtaCutEff.gif");



  //*******************************************************************************************
  //BJet Pt
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  hist_bjet1Pt->Draw("hist");
  cv->SaveAs("GenBjet1Pt.gif");

  cv = new TCanvas("cv","cv", 800,600);
  TH1F *hist_bjet1PtCutEff = MakeCutEfficiencyPlot(hist_bjet1Pt, true);
  hist_bjet1PtCutEff->Draw("hist");
  cv->SaveAs("GenBjet1PtCutEff.gif");
  

  cv = new TCanvas("cv","cv", 800,600);
  hist_bjet2Pt->Draw("hist");
  cv->SaveAs("GenBjet2Pt.gif");

  cv = new TCanvas("cv","cv", 800,600);
  TH1F *hist_bjet2PtCutEff = MakeCutEfficiencyPlot(hist_bjet2Pt, true);
  hist_bjet2PtCutEff->Draw("hist");
  cv->SaveAs("GenBjet2PtCutEff.gif");


  //*******************************************************************************************
  //B-Jet Eta
  //*******************************************************************************************

  cv = new TCanvas("cv","cv", 800,600);
  hist_bjetEta->Draw("hist");
  cv->SaveAs("GenBjetEta.gif");

  cv = new TCanvas("cv","cv", 800,600);
  hist_bjetMaxEta->Draw("hist");
  cv->SaveAs("GenBjetMaxEta.gif");

  cv = new TCanvas("cv","cv", 800,600);
  TH1F *hist_bjetMaxEtaCutEff = MakeCutEfficiencyPlot(hist_bjetMaxEta, false);
  hist_bjetMaxEtaCutEff->Draw("hist");
  legend->Draw();
  cv->SaveAs("GenBjetMaxEtaCutEff.gif");



  //*******************************************************************************************
  //Photon Pt and Eta After Baseline BJet Cuts
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  TH1F *hist_photon2PtCutEff_AfterBaselineBJetCuts = MakeCutEfficiencyPlot(hist_photon2Pt_AfterBaselineBJetCuts, true);

  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(hist_photon2PtCutEff, "All Events", "L");
  legend->AddEntry(hist_photon2PtCutEff_AfterBaselineBJetCuts, "Events Passing b-Jet Acceptance", "L");

  hist_photon2PtCutEff_AfterBaselineBJetCuts->SetLineColor(kBlue);
  hist_photon2PtCutEff->Draw("hist");
  hist_photon2PtCutEff_AfterBaselineBJetCuts->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("GenPhoton2PtCutEff_AfterBaselineBJetCuts.gif");

  cv = new TCanvas("cv","cv", 800,600);
  TH1F *hist_photonMaxEtaCutEff_AfterBaselineBJetCuts = MakeCutEfficiencyPlot(hist_photonMaxEta_AfterBaselineBJetCuts, false);

  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(hist_photonMaxEtaCutEff, "All Events", "L");
  legend->AddEntry(hist_photonMaxEtaCutEff_AfterBaselineBJetCuts, "Events Passing b-Jet Acceptance", "L");

  hist_photonMaxEtaCutEff_AfterBaselineBJetCuts->SetLineColor(kBlue);
  hist_photonMaxEtaCutEff->Draw("hist");
  hist_photonMaxEtaCutEff_AfterBaselineBJetCuts->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("GenPhotonMaxEtaCutEff_AfterBaselineBJetCuts.gif");


  //*******************************************************************************************
  //BJet Pt and Eta After Baseline Photon Cuts
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  TH1F *hist_bjet2PtCutEff_AfterBaselinePhotonCuts = MakeCutEfficiencyPlot(hist_bjet2Pt_AfterBaselinePhotonCuts, true);

  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(hist_bjet2PtCutEff, "All Events", "L");
  legend->AddEntry(hist_bjet2PtCutEff_AfterBaselinePhotonCuts, "Events Passing Photon Acceptance", "L");

  hist_bjet2PtCutEff_AfterBaselinePhotonCuts->SetLineColor(kBlue);
  hist_bjet2PtCutEff->Draw("hist");
  hist_bjet2PtCutEff_AfterBaselinePhotonCuts->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("GenBjet2PtCutEff_AfterBaselinePhotonCuts.gif");

  cv = new TCanvas("cv","cv", 800,600);
  TH1F *hist_bjetMaxEtaCutEff_AfterBaselinePhotonCuts = MakeCutEfficiencyPlot(hist_bjetMaxEta_AfterBaselinePhotonCuts, false);

  legend = new TLegend(0.54,0.24,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(hist_bjetMaxEtaCutEff, "All Events", "L");
  legend->AddEntry(hist_bjetMaxEtaCutEff_AfterBaselinePhotonCuts, "Events Passing Photon Acceptance", "L");

  hist_bjetMaxEtaCutEff_AfterBaselinePhotonCuts->SetLineColor(kBlue);
  hist_bjetMaxEtaCutEff->Draw("hist");
  hist_bjetMaxEtaCutEff_AfterBaselinePhotonCuts->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("GenBJetMaxEtaCutEff_AfterBaselinePhotonCuts.gif");





  //*******************************************************************************************
  //Summarize Yields and Acceptance
  //*******************************************************************************************
  cout << "Total Events: " << TotalEvents << "\n\n";


  cout << "EventsPassPhotonPt_40_25 : " << EventsPassPhotonPt_40_25 << " | " << EventsPassPhotonPt_40_25/TotalEvents << "\n";
  cout << "EventsPassPhotonPt_25_25 : " << EventsPassPhotonPt_25_25 << " | " << EventsPassPhotonPt_25_25/TotalEvents << "\n";
  cout << "EventsPassPhotonPt_40_20 : " << EventsPassPhotonPt_40_20 << " | " << EventsPassPhotonPt_40_20/TotalEvents << "\n";
  cout << "EventsPassPhotonPt_40_15 : " << EventsPassPhotonPt_40_15 << " | " << EventsPassPhotonPt_40_15/TotalEvents << "\n";
  cout << "EventsPassPhotonPt_50_20 : " << EventsPassPhotonPt_50_20 << " | " << EventsPassPhotonPt_50_20/TotalEvents << "\n";
  cout << "EventsPassPhotonPt_60_15 : " << EventsPassPhotonPt_60_15 << " | " << EventsPassPhotonPt_60_15/TotalEvents << "\n";


  cout << "\n\n";
  cout << "EventsPassPhotonPt_40_25_BothEta2p5 : " << EventsPassPhotonPt_40_25_BothEta2p5 << " | " << EventsPassPhotonPt_40_25_BothEta2p5/TotalEvents << " | baseline " << "\n";
  cout << "EventsPassPhotonPt_40_25_BothEta3p0 : " << EventsPassPhotonPt_40_25_BothEta3p0 << " | " << EventsPassPhotonPt_40_25_BothEta3p0/TotalEvents << " | " << EventsPassPhotonPt_40_25_BothEta3p0/EventsPassPhotonPt_40_25_BothEta2p5 << "\n";
  cout << "EventsPassPhotonPt_40_25_BothEta3p5 : " << EventsPassPhotonPt_40_25_BothEta3p5 << " | " << EventsPassPhotonPt_40_25_BothEta3p5/TotalEvents << " | " << EventsPassPhotonPt_40_25_BothEta3p5/EventsPassPhotonPt_40_25_BothEta2p5 << "\n";
  cout << "EventsPassPhotonPt_40_25_BothEta4p0 : " << EventsPassPhotonPt_40_25_BothEta4p0 << " | " << EventsPassPhotonPt_40_25_BothEta4p0/TotalEvents << " | " << EventsPassPhotonPt_40_25_BothEta4p0/EventsPassPhotonPt_40_25_BothEta2p5 << "\n";
  cout << "EventsPassPhotonPt_40_25_BothEta2p0 : " << EventsPassPhotonPt_40_25_BothEta2p0 << " | " << EventsPassPhotonPt_40_25_BothEta2p0/TotalEvents << " | baseline " << "\n";
  cout << "EventsPassPhotonPt_40_25_BothEta1p5 : " << EventsPassPhotonPt_40_25_BothEta1p5 << " | " << EventsPassPhotonPt_40_25_BothEta1p5/TotalEvents << " | baseline " << "\n";


  cout << "\n\n";
  cout << "EventsPassBJetPt_30_30 : " << EventsPassBJetPt_30_30 << " | " << EventsPassBJetPt_30_30/TotalEvents << "\n";
  cout << "EventsPassBJetPt_25_25 : " << EventsPassBJetPt_25_25 << " | " << EventsPassBJetPt_25_25/TotalEvents << "\n";
  cout << "EventsPassBJetPt_30_25 : " << EventsPassBJetPt_30_25 << " | " << EventsPassBJetPt_30_25/TotalEvents << "\n";
  cout << "EventsPassBJetPt_30_20 : " << EventsPassBJetPt_30_20 << " | " << EventsPassBJetPt_30_20/TotalEvents << "\n";
  cout << "EventsPassBJetPt_45_25 : " << EventsPassBJetPt_45_25 << " | " << EventsPassBJetPt_45_25/TotalEvents << "\n";
  cout << "EventsPassBJetPt_60_20 : " << EventsPassBJetPt_60_20 << " | " << EventsPassBJetPt_60_20/TotalEvents << "\n";


  cout << "\n\n";
  cout << "EventsPassBJetPt_30_30_BothEta2p4 : " << EventsPassBJetPt_30_30_BothEta2p4 << " | " << EventsPassBJetPt_30_30_BothEta2p4/TotalEvents << " | baseline " << "\n";
  cout << "EventsPassBJetPt_30_30_BothEta3p0 : " << EventsPassBJetPt_30_30_BothEta3p0 << " | " << EventsPassBJetPt_30_30_BothEta3p0/TotalEvents << " | " << EventsPassBJetPt_30_30_BothEta3p0/EventsPassBJetPt_30_30_BothEta2p4 << "\n";
  cout << "EventsPassBJetPt_30_30_BothEta3p5 : " << EventsPassBJetPt_30_30_BothEta3p5 << " | " << EventsPassBJetPt_30_30_BothEta3p5/TotalEvents << " | " << EventsPassBJetPt_30_30_BothEta3p5/EventsPassBJetPt_30_30_BothEta2p4 << "\n";
  cout << "EventsPassBJetPt_30_30_BothEta4p0 : " << EventsPassBJetPt_30_30_BothEta4p0 << " | " << EventsPassBJetPt_30_30_BothEta4p0/TotalEvents << " | " << EventsPassBJetPt_30_30_BothEta4p0/EventsPassBJetPt_30_30_BothEta2p4 << "\n";


  cout << "\n\n";
  cout << "EventsPass_PhotonPt4025Eta2p5_BJetPt3030Eta2p4 : " << EventsPass_PhotonPt4025Eta2p5_BJetPt3030Eta2p4 << " | " << EventsPass_PhotonPt4025Eta2p5_BJetPt3030Eta2p4/TotalEvents << " | baseline " << "\n";
  cout << "EventsPass_PhotonPt4025Eta2p5_BJetPt3030Eta3p0 : " << EventsPass_PhotonPt4025Eta2p5_BJetPt3030Eta3p0 << " | " << EventsPass_PhotonPt4025Eta2p5_BJetPt3030Eta3p0/TotalEvents << " | " << EventsPass_PhotonPt4025Eta2p5_BJetPt3030Eta3p0/EventsPass_PhotonPt4025Eta2p5_BJetPt3030Eta2p4 << "\n";
  cout << "EventsPass_PhotonPt4025Eta3p0_BJetPt3030Eta2p4 : " << EventsPass_PhotonPt4025Eta3p0_BJetPt3030Eta2p4 << " | " << EventsPass_PhotonPt4025Eta3p0_BJetPt3030Eta2p4/TotalEvents << " | " << EventsPass_PhotonPt4025Eta3p0_BJetPt3030Eta2p4/EventsPass_PhotonPt4025Eta2p5_BJetPt3030Eta2p4 << "\n";
  cout << "EventsPass_PhotonPt4025Eta3p0_BJetPt3030Eta3p0 : " << EventsPass_PhotonPt4025Eta3p0_BJetPt3030Eta3p0 << " | " << EventsPass_PhotonPt4025Eta3p0_BJetPt3030Eta3p0/TotalEvents << " | " << EventsPass_PhotonPt4025Eta3p0_BJetPt3030Eta3p0/EventsPass_PhotonPt4025Eta2p5_BJetPt3030Eta2p4 << "\n";
  cout << "EventsPass_PhotonPt4025Eta2p0_BJetPt3030Eta2p4 : " << EventsPass_PhotonPt4025Eta2p0_BJetPt3030Eta2p4 << " | " << EventsPass_PhotonPt4025Eta2p0_BJetPt3030Eta2p4/TotalEvents << " | baseline " << "\n";
  cout << "EventsPass_PhotonPt4025Eta1p5_BJetPt3030Eta2p4 : " << EventsPass_PhotonPt4025Eta1p5_BJetPt3030Eta2p4 << " | " << EventsPass_PhotonPt4025Eta1p5_BJetPt3030Eta2p4/TotalEvents << " | " << EventsPass_PhotonPt4025Eta1p5_BJetPt3030Eta2p4/EventsPass_PhotonPt4025Eta2p5_BJetPt3030Eta2p4 << "\n";


  cout << "\n";
  cout << "EventsPass_PhotonPt6015Eta2p5_BJetPt6020Eta2p4 : " << EventsPass_PhotonPt6015Eta2p5_BJetPt6020Eta2p4 << " | " << EventsPass_PhotonPt6015Eta2p5_BJetPt6020Eta2p4/TotalEvents << " | baseline " << "\n";
  cout << "EventsPass_PhotonPt6015Eta2p5_BJetPt6020Eta3p0 : " << EventsPass_PhotonPt6015Eta2p5_BJetPt6020Eta3p0 << " | " << EventsPass_PhotonPt6015Eta2p5_BJetPt6020Eta3p0/TotalEvents << " | " << EventsPass_PhotonPt6015Eta2p5_BJetPt6020Eta3p0/EventsPass_PhotonPt6015Eta2p5_BJetPt6020Eta2p4 << "\n";
  cout << "EventsPass_PhotonPt6015Eta3p0_BJetPt6020Eta2p4 : " << EventsPass_PhotonPt6015Eta3p0_BJetPt6020Eta2p4 << " | " << EventsPass_PhotonPt6015Eta3p0_BJetPt6020Eta2p4/TotalEvents << " | " << EventsPass_PhotonPt6015Eta3p0_BJetPt6020Eta2p4/EventsPass_PhotonPt6015Eta2p5_BJetPt6020Eta2p4 << "\n";
  cout << "EventsPass_PhotonPt6015Eta3p0_BJetPt6020Eta3p0 : " << EventsPass_PhotonPt6015Eta3p0_BJetPt6020Eta3p0 << " | " << EventsPass_PhotonPt6015Eta3p0_BJetPt6020Eta3p0/TotalEvents << " | " << EventsPass_PhotonPt6015Eta3p0_BJetPt6020Eta3p0/EventsPass_PhotonPt6015Eta2p5_BJetPt6020Eta2p4 << "\n";



  
}

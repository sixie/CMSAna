#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include "CMSAna/HHToBBGG/interface/HHToBBGGEventTree.hh"

//------------------------------------------------------------------------------
// PlotHiggsRes_LP
//------------------------------------------------------------------------------
void ExampleAnalysis ( string InputFilename    = "HHToBBGG.14TeV.root")
{

  TH1F* diphotonMass = new TH1F( "diphotonMass", ";M_{#gamma#gamma};Number of Events", 100, 100, 150);

  cmsana::HHToBBGGEventTree event;
  event.LoadTree(InputFilename.c_str());
  event.InitTree();

  for (int n=0;n<event.tree_->GetEntries();n++) { 
    event.tree_->GetEntry(n);

    if (event.bjet1.Pt() > 25 && event.bjet2.Pt() > 25
        && event.pho1.Pt() > 25 && event.pho2.Pt() > 25) {

      diphotonMass->Fill( event.diphoton.M() , event.weight);
    }
    
  }

  diphotonMass->Draw();

}
